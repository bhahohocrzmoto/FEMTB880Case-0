# -*- coding: utf-8 -*-
"""IEC 60287-1-1 loss engine for the TB 880 Case #0 benchmark workflow.

I use this module as the single authoritative implementation of the IEC 60287-1-1
loss chain for CIGRE TB 880 Case #0. Both the analytical solver and the ANSYS
Picard FEM loop call into the same API so that conductor, dielectric, and sheath
losses are always computed with one traceable formulation.
"""
# ============================================================
# cable_losses_tb880_case0.py
# ------------------------------------------------------------
# IEC 60287-1-1 style losses for TB 880 Case #0 (touching trefoil)
#
# Default repository interpretation:
# - benchmark-faithful TB 880 Case #0-1 IEC base path
# - solid bonding includes circulating sheath losses
# - solid-bonded eddy sheath losses are neglected by default
#
# Optional comparison mode is retained via sheath_eddy_policy.
# ============================================================

import math
from tb880_case0_data import CASE


def normalize_bonding_mode(bonding):
    """Normalize bonding mode to: solid, single, cross (legacy alias: both->solid)."""
    mode = str(bonding).strip().lower()
    legacy = {"both": "solid"}
    normalized = legacy.get(mode, mode)
    if normalized not in ("solid", "single", "cross"):
        raise ValueError("Unsupported bonding mode '{0}'. Use: solid, single, cross".format(bonding))
    return normalized


def normalize_sheath_eddy_policy(policy):
    """Normalize sheath eddy policy mode."""
    mode = str(policy).strip().lower()
    valid = (
        "auto",
        "exclude",
        "include_for_solid",
        "include_for_all",
    )
    if mode not in valid:
        raise ValueError(
            "Unsupported sheath_eddy_policy '{0}'. Use: {1}".format(policy, ", ".join(valid))
        )
    return mode


class Cable(object):
    def __init__(self, cid,
                 d_cond, d_inner_semicon, d_ins, d_outer_semicon, d_screen_mean, d_screen_out, d_oversheath,
                 A_cond_elec, A_screen_elec, A_cond_FE, A_innerins_FE, A_screen_FE,
                 r_cond_dc_20_ohm_per_m,
                 rho_cond_20, alpha_cond, ks, kp,
                 rho_screen_20, alpha_screen,
                 eps_r, tan_delta,
                 f_hz, U_LL_V, I_rms_A,
                 bonding="solid",
                 sheath_eddy_policy="auto"):
        self.cid = cid
        self.I = float(I_rms_A)
        self.f = float(f_hz)
        self.U_LL_V = float(U_LL_V)

        self.d_cond = d_cond
        self.d_inner_semicon = d_inner_semicon
        self.d_ins = d_ins
        self.d_outer_semicon = d_outer_semicon
        self.d_screen_mean = d_screen_mean
        self.d_screen_out = d_screen_out
        self.d_oversheath = d_oversheath

        self.A_cond_elec = A_cond_elec
        self.A_screen_elec = A_screen_elec
        self.A_cond_FE = A_cond_FE
        self.A_innerins_FE = A_innerins_FE
        self.A_screen_FE = A_screen_FE

        # I preserve the benchmark input R_dc(20 degC) because TB 880 takes this as
        # a primary conductor datum and IEC 60287-1-1 Eq. (1) updates it with temperature.
        self.r_cond_dc_20_ohm_per_m = float(r_cond_dc_20_ohm_per_m)
        self.rho_cond_20 = rho_cond_20
        self.alpha_cond = alpha_cond
        self.ks = ks
        self.kp = kp
        self.rho_screen_20 = rho_screen_20
        self.alpha_screen = alpha_screen
        self.eps_r = eps_r
        self.tan_delta = tan_delta

        self.Rac = None
        self.Rdc = None
        self.Rs = None
        self.lambda1 = None

        self.bonding = normalize_bonding_mode(bonding)
        self.sheath_eddy_policy = normalize_sheath_eddy_policy(sheath_eddy_policy)

    def _resistance_ratio(self, alpha, T_C):
        # I apply IEC 60287-1-1:2023 Eq. (1) in ratio form so that the same linear
        # temperature correction can be reused for both conductor and sheath metals.
        return 1.0 + alpha * (T_C - 20.0)

    def Rdc_at_temp(self, T_cond_C):
        """Apply IEC 60287-1-1 Eq. (1) to the stored benchmark R_dc(20 degC)."""
        return self.r_cond_dc_20_ohm_per_m * self._resistance_ratio(self.alpha_cond, T_cond_C)

    def Rac_at_temp(self, T_cond_C):
        # I first compute the conductor DC resistance at operating temperature using
        # IEC 60287-1-1 Eq. (1), because Eq. (2), Eq. (5), and Eq. (9) all depend on it.
        Rdc = self.Rdc_at_temp(T_cond_C)

        # I evaluate x_s^2 from IEC 60287-1-1 Eq. (5), then select the appropriate
        # skin-effect branch from Eq. (6), Eq. (7), or Eq. (8) based on the x_s range.
        xs2 = (8.0 * math.pi * self.f * 1e-7 * self.ks) / Rdc
        xs = math.sqrt(xs2)
        if xs <= 2.8:
            ys = xs**4 / (192.0 + 0.8 * xs**4)
        elif xs <= 3.8:
            ys = -0.136 - 0.0177 * xs + 0.0563 * xs**2
        else:
            ys = 0.354 * xs - 0.733

        # I evaluate x_p^2 from IEC 60287-1-1 Eq. (9), then form g with the same
        # rational expression used in Eq. (10) for the proximity-effect calculation.
        xp2 = (8.0 * math.pi * self.f * 1e-7 * self.kp) / Rdc
        xp = math.sqrt(xp2)
        g = xp**4 / (192.0 + 0.8 * xp**4)
        # In TB 880 touching trefoil, I take s as the overall cable diameter so that
        # the Eq. (10) geometry factor (d_c / s)^2 reflects touching installation.
        s = self.d_oversheath
        r = (self.d_cond / s) ** 2
        # I apply IEC 60287-1-1 Eq. (10) for three single-core cables in trefoil,
        # which adds the proximity term y_p to the AC resistance formulation.
        yp = g * r * (0.312 * r + 1.18 / (g + 0.27))

        # I complete IEC 60287-1-1 Eq. (2) by multiplying R_dc(theta) by the sum of
        # the skin and proximity correction factors.
        return Rdc * (1.0 + ys + yp)

    def Rs_at_temp(self, T_screen_C):
        # I compute the sheath DC resistance from material resistivity and nominal
        # electrical area because IEC 60287-1-1 models the metallic sheath this way.
        rho_screen = self.rho_screen_20 * self._resistance_ratio(self.alpha_screen, T_screen_C)
        return rho_screen / self.A_screen_elec

    def sheath_resistance_at_temp(self, T_screen_C):
        return self.Rs_at_temp(T_screen_C)

    def sheath_reactance_X(self):
        # I use the IEC 60287-1-1 Section 2.3.1 trefoil reactance expression
        # X = 2 * omega * 1e-7 * ln(2s/d), where s is centre spacing and d is the
        # mean sheath diameter, because circulating current losses depend on X.
        s = self.d_oversheath
        d = self.d_screen_mean
        if (2.0 * s) <= d:
            raise ValueError("Need 2*s > d for ln(2s/d)")
        omega = 2.0 * math.pi * self.f
        return 2.0 * omega * 1e-7 * math.log((2.0 * s) / d)

    def _lambda1_prime(self, Rs, Rac, X):
        # I apply IEC 60287-1-1 Section 2.3.1 for the solid-bonded circulating
        # sheath loss factor lambda1' in touching trefoil.
        return (Rs / Rac) * (1.0 / (1.0 + (Rs / X) ** 2))

    def _lambda1_doubleprime(self, Rs, Rac):
        # I implement IEC 60287-1-1 Section 2.3.5 and Section 2.3.6, where the eddy
        # term depends on m = omega * 1e-7 / R_s and on the geometric factor lambda_0.
        omega = 2.0 * math.pi * self.f
        m = omega * 1e-7 / Rs

        d = self.d_screen_mean
        s = self.d_oversheath

        # I use lambda_0 = 3 * [m^2 / (1 + m^2)] * (d / (2s))^2 from IEC 60287-1-1
        # as the base eddy-current coupling term for single-core cables.
        lambda0 = 3.0 * (m**2 / (1.0 + m**2)) * (d / (2.0 * s)) ** 2
        # I keep delta_1 explicitly because IEC 60287-1-1 applies this correction to
        # account for the non-ideal distribution in the bonded sheath loss formula.
        delta1 = (1.14 * m**2.45 + 0.33) * (d / (2.0 * s)) ** (0.92 * m + 1.66)
        # I set delta_2 to zero because there is no additional correction active for
        # this thin single-layer sheath representation in the benchmark case.
        delta2 = 0.0
        # I use g_s = 1.0 for a single metallic layer, which matches the simplified
        # TB 880 sheath model used in this repository.
        gs = 1.0
        # I set the thickness term to zero because I treat the laminated sheath as a
        # thin screen and therefore neglect the extra thickness correction.
        thickness_term = 0.0

        lambda1pp = (Rs / Rac) * (gs * lambda0 * (1.0 + delta1 + delta2) + thickness_term)

        if self.bonding == "solid":
            X = self.sheath_reactance_X()
            # For solid bonding, IEC 60287-1-1 applies the reduction factor F based
            # on M = R_s / X so that eddy losses are adjusted for the bonding circuit.
            M = Rs / X
            F = (4.0 * M**4 + (2.0 * M)**2) / (4.0 * (M**2 + 1.0)**2)
            lambda1pp *= F

        return lambda1pp

    def _include_eddy_for_current_setup(self):
        # I centralize this policy logic because the repository keeps two use cases:
        # benchmark-faithful TB 880 regression, and comparison studies that include
        # eddy losses to examine their impact on the current rating.
        mode = self.sheath_eddy_policy
        if mode == "exclude":
            return False
        if mode == "include_for_solid":
            return self.bonding == "solid"
        if mode == "include_for_all":
            return True

        # auto -> repository default interpretation from CASE assumptions.
        # I exclude solid-bonded eddy losses by default because the published TB 880
        # Case #0-1 results align with circulating-current losses only.
        if self.bonding == "solid":
            return bool(CASE.assumptions.solid_bonded_include_eddy_default)
        return True

    def dielectric_loss_W_per_m(self):
        # I use IEC 60287-1-1 Eq. (15) for capacitance, with the semicon screen
        # diameters D_i and d_c taken from the stored TB 880 geometry definition.
        Di = self.d_outer_semicon
        d = self.d_inner_semicon
        if Di <= d:
            raise Exception("Need Di > d")
        C = 2.0 * math.pi * CASE.constants.eps0_f_per_m * self.eps_r / math.log(Di / d)

        # I then apply IEC 60287-1-1 Eq. (14), using U_0 = U_LL / sqrt(3), because
        # dielectric losses in a single-core cable depend on phase-to-earth voltage.
        U0 = self.U_LL_V / math.sqrt(3.0)
        omega = 2.0 * math.pi * self.f
        return omega * C * U0**2 * self.tan_delta

    def calculate_losses(self, T_core_C, T_screen_C):
        """Single IEC 60287 loss chain used by the analytical and FEM workflows."""
        # I begin with the conductor resistance chain R_dc -> R_ac because conductor
        # Joule loss W_c = I^2 * R_ac is the dominant source term in the benchmark.
        Rac = self.Rac_at_temp(T_core_C)
        Rdc = self.Rdc_at_temp(T_core_C)
        Rs = self.Rs_at_temp(T_screen_C)

        # I evaluate dielectric loss independently from current because IEC 60287-1-1
        # Eq. (14) depends on voltage, capacitance, and tan(delta), not on load current.
        Wd = self.dielectric_loss_W_per_m()
        Wc = self.I**2 * Rac

        # I next compute sheath losses through lambda1' and lambda1'' because IEC
        # 60287-1-1 defines W_s = lambda1 * W_c for the metallic sheath.
        X = self.sheath_reactance_X()
        lambda1_prime = self._lambda1_prime(Rs, Rac, X) if self.bonding == "solid" else 0.0
        lambda1_doubleprime = self._lambda1_doubleprime(Rs, Rac) if self._include_eddy_for_current_setup() else 0.0

        # I combine the sheath terms according to IEC 60287-1-1 Section 2.3.7:
        # solid bonding uses lambda1' + lambda1'', while single-point or cross
        # bonding keeps only the eddy component because circulating current is absent.
        if self.bonding == "solid":
            lambda1 = lambda1_prime + lambda1_doubleprime
        else:
            lambda1 = lambda1_doubleprime

        Ws = lambda1 * Wc

        # I cache these values so that downstream solvers can inspect the last loss
        # state without recomputing the full IEC chain.
        self.Rac = Rac
        self.Rdc = Rdc
        self.Rs = Rs
        self.lambda1 = lambda1
        self.Wc = Wc
        self.Ws = Ws
        self.Wd = Wd

        return {
            "core": Wc,
            "screen": Ws,
            "dielectric": Wd,
            "lambda1": lambda1,
            "lambda1_prime": lambda1_prime,
            "lambda1_doubleprime": lambda1_doubleprime,
            "Rac": Rac,
            "Rdc": Rdc,
            "Rs": Rs,
        }

    def q_core(self, T_core_C, T_screen_C):
        # I convert conductor loss from W/m to W/m^3 by dividing by the geometric FE
        # copper area, which follows the volumetric source convention described in TB 963.
        losses = self.calculate_losses(T_core_C, T_screen_C)
        return losses["core"] / self.A_cond_FE

    def q_screen(self, T_core_C, T_screen_C):
        # I convert sheath loss from W/m to W/m^3 in the same way so that ANSYS can
        # apply the metallic screen heating as a body load.
        losses = self.calculate_losses(T_core_C, T_screen_C)
        return losses["screen"] / self.A_screen_FE

    def q_innerins(self, T_core_C, T_screen_C):
        # I convert dielectric loss from W/m to W/m^3 using the FE inner-insulation
        # annulus area because TB 963 couples IEC losses back as volumetric heating.
        losses = self.calculate_losses(T_core_C, T_screen_C)
        return losses["dielectric"] / self.A_innerins_FE


def create_case0_cables(I_rms_A=None, bonding=None, sheath_eddy_policy="auto"):
    """Return three symmetric trefoil Cable objects built from the centralized CASE data."""
    if I_rms_A is None:
        I_rms_A = CASE.benchmark.i_final_a
    if bonding is None:
        bonding = CASE.installation.bonding
    bonding = normalize_bonding_mode(bonding)

    d_core = CASE.geometry.d_cond_m
    d_inner_semicon = CASE.geometry.d_inner_semicon_m
    d_ins = CASE.geometry.d_ins_m
    d_outer_semicon = CASE.geometry.d_outer_semicon_m
    d_screen_mean = CASE.geometry.d_screen_mean_m
    d_screen_out = CASE.geometry.d_screen_out_m
    d_oversheath = CASE.geometry.d_oversheath_m

    A_cond_elec = CASE.geometry.a_cond_elec_m2
    A_screen_elec = CASE.geometry.a_screen_elec_m2
    A_cond_FE = CASE.geometry.a_cond_geom_m2
    A_innerins_FE = CASE.geometry.a_innerins_geom_m2
    A_screen_FE = CASE.geometry.a_screen_geom_m2

    r_cond_dc_20_ohm_per_m = CASE.material.r_cond_dc_20_ohm_per_m
    rho_cond_20 = CASE.material.rho_cond_20_ohm_m
    alpha_cond = CASE.material.alpha_cond_20_per_k
    ks = CASE.material.ks
    kp = CASE.material.kp

    rho_screen_20 = CASE.material.rho_screen_20_ohm_m
    alpha_screen = CASE.material.alpha_screen_20_per_k

    eps_r = CASE.material.eps_r
    tan_delta = CASE.material.tan_delta

    f_hz = CASE.installation.frequency_hz
    U_LL_V = CASE.installation.u_ll_v

    # I create three identical Cable objects because the touching trefoil benchmark is
    # symmetric and each phase shares the same centralized geometry and material data.
    cables = {}
    for cid in ["C01", "C02", "C03"]:
        cables[cid] = Cable(
            cid=cid,
            d_cond=d_core,
            d_inner_semicon=d_inner_semicon,
            d_ins=d_ins,
            d_outer_semicon=d_outer_semicon,
            d_screen_mean=d_screen_mean,
            d_screen_out=d_screen_out,
            d_oversheath=d_oversheath,
            A_cond_elec=A_cond_elec,
            A_screen_elec=A_screen_elec,
            A_cond_FE=A_cond_FE,
            A_innerins_FE=A_innerins_FE,
            A_screen_FE=A_screen_FE,
            r_cond_dc_20_ohm_per_m=r_cond_dc_20_ohm_per_m,
            rho_cond_20=rho_cond_20,
            alpha_cond=alpha_cond,
            ks=ks,
            kp=kp,
            rho_screen_20=rho_screen_20,
            alpha_screen=alpha_screen,
            eps_r=eps_r,
            tan_delta=tan_delta,
            f_hz=f_hz,
            U_LL_V=U_LL_V,
            I_rms_A=I_rms_A,
            bonding=bonding,
            sheath_eddy_policy=sheath_eddy_policy,
        )
    return cables
