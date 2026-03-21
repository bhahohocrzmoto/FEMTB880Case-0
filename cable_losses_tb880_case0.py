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
# - central CASE flags include circulating sheath losses by default
# - central CASE flags exclude eddy sheath losses by default
# ============================================================

import math
from tb880_case0_data import CASE



class Cable(object):
    def __init__(self, cid,
                 d_cond, d_inner_semicon, d_ins, d_outer_semicon, d_screen_mean, d_screen_out, d_oversheath,
                 A_cond_elec, A_screen_elec, A_cond_FE_geom, A_innerins_FE_geom, A_screen_FE_geom,
                 r_cond_dc_20_ohm_per_m,
                 rho_cond_20, alpha_cond, ks, kp,
                 rho_screen_20, alpha_screen,
                 eps_r, tan_delta,
                 f_hz, U_LL_V, I_rms_A):
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
        self.A_cond_FE_geom = A_cond_FE_geom
        self.A_innerins_FE_geom = A_innerins_FE_geom
        self.A_screen_FE_geom = A_screen_FE_geom
        self.A_cond_FE_model = A_cond_FE_geom
        self.A_innerins_FE_model = A_innerins_FE_geom
        self.A_screen_FE_model = A_screen_FE_geom

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

    def _get_area(self, region, area_mode):
        """
        Return the cross-sectional area for `region` based on `area_mode`.

        Parameters
        ----------
        region : str
            One of "conductor", "screen", "innerins".
        area_mode : str
            "analytical" -> use the IEC nominal/electrical areas stored in the
            central data file (A_cond_elec, A_screen_elec). These are the areas
            used for IEC 60287 resistance and loss-factor calculations.
            "fem" -> use the FEM model areas (A_cond_FE_model,
            A_screen_FE_model, A_innerins_FE_model). These are initially set from
            geometry-derived reference values, and overwritten at runtime by the
            ANSYS driver with actual FE body areas.

        Returns
        -------
        float
            Cross-sectional area in m^2.
        """
        if area_mode == "analytical":
            mapping = {
                "conductor": self.A_cond_elec,
                "screen": self.A_screen_elec,
                "innerins": None,
            }
        elif area_mode == "fem":
            mapping = {
                "conductor": self.A_cond_FE_model,
                "screen": self.A_screen_FE_model,
                "innerins": self.A_innerins_FE_model,
            }
        else:
            raise ValueError(
                "area_mode must be 'analytical' or 'fem', got '{0}'".format(area_mode)
            )
        area = mapping.get(region)
        if area is None:
            raise ValueError(
                "No {0} area defined for region '{1}'".format(area_mode, region)
            )
        return area

    def Rs_at_temp(self, T_screen_C, area_mode="analytical"):
        # I compute the sheath DC resistance from material resistivity and the
        # selected screen area so analytical and FEM workflows cannot mix modes.
        rho_screen = self.rho_screen_20 * self._resistance_ratio(self.alpha_screen, T_screen_C)
        return rho_screen / self._get_area("screen", area_mode)

    def sheath_resistance_at_temp(self, T_screen_C, area_mode="analytical"):
        return self.Rs_at_temp(T_screen_C, area_mode=area_mode)

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

    def _beta1(self, Rs, area_mode="analytical"):
        # IEC 60287-1-1:2023 Sec 5.3.7.1 defines beta_1 = sqrt(4 * pi * omega /
        # (rho_s * 1e7)), where rho_s is the sheath resistivity in ohm.m recovered
        # from R_s = rho_s / A_s using the selected screen area for the active mode.
        rho_s = Rs * self._get_area("screen", area_mode)
        omega = 2.0 * math.pi * self.f
        return math.sqrt(4.0 * math.pi * omega / (rho_s * 1e7))

    def _gs_factor(self, beta1):
        # IEC 60287-1-1:2023 Sec 5.3.7.1 defines C_gs = 1 + (t_s / D_s)^1.74 *
        # (beta_1 * D_s * 1e-3 - 1.6), with sheath thickness t_s in mm and D_s as
        # the external diameter of the metallic sheath in mm.
        t_s_mm = 0.5 * (self.d_screen_out - self.d_outer_semicon) * 1000.0
        D_s_mm = self.d_screen_out * 1000.0
        ratio = t_s_mm / D_s_mm
        return 1.0 + ratio ** 1.74 * (beta1 * D_s_mm * 1e-3 - 1.6)

    def _eddy_thickness_term(self, beta1):
        # IEC 60287-1-1:2023 Sec 5.3.7.1 adds (beta_1 * t_s)^4 / (12 * 1e12), with
        # sheath thickness t_s taken in mm using the metallic sheath radial build.
        t_s_mm = 0.5 * (self.d_screen_out - self.d_outer_semicon) * 1000.0
        return (beta1 * t_s_mm) ** 4 / (12.0 * 1e12)

    def _lambda1_doubleprime(self, Rs, Rac, apply_F=False):
        # I implement IEC 60287-1-1:2023 Sec 5.3.7.1 for sheath eddy-current
        # losses. The parameter m = omega * 1e-7 / R_s quantifies the sheath
        # penetration depth relative to its thickness. lambda_0 is the base
        # trefoil coupling term, and beta_1 / C_gs / thickness use the
        # metallic sheath geometry from the stored cable diameters.
        omega = 2.0 * math.pi * self.f
        m = omega * 1e-7 / Rs

        d = self.d_screen_mean
        s = self.d_oversheath

        # I use lambda_0 = 3 * [m^2 / (1 + m^2)] * (d / (2s))^2 from
        # IEC 60287-1-1:2023 Sec 5.3.7.1 as the base eddy-current coupling
        # term for three single-core cables in trefoil formation.
        # The coefficient 3 is C_C0 for trefoil in IEC 60287-1-1 Sec 5.3.7.1.
        lambda0 = 3.0 * (m**2 / (1.0 + m**2)) * (d / (2.0 * s)) ** 2

        # I evaluate Delta_1 from IEC 60287-1-1:2023 Sec 5.3.7.1 as the
        # first correction to the base eddy-current coupling term.
        delta1 = (1.14 * m**2.45 + 0.33) * (d / (2.0 * s)) ** (0.92 * m + 1.66)

        # IEC 60287-1-1:2023 Sec 5.3.7.1: Delta_2 = 0 for three single-core
        # cables in trefoil formation. Non-zero Delta_2 applies only to flat
        # formations, namely centre cable, outer leading, and outer lagging.
        delta2 = 0.0

        # IEC 60287-1-1:2023 Sec 5.3.7.1 defines beta_1 from the sheath
        # resistivity and angular frequency.
        beta1 = self._beta1(Rs)

        # IEC 60287-1-1:2023 Sec 5.3.7.1 defines C_gs from beta_1, sheath
        # thickness t_s, and sheath outer diameter D_s.
        gs = self._gs_factor(beta1)

        # IEC 60287-1-1:2023 Sec 5.3.7.1 adds (beta_1 * t_s)^4 / (12 * 1e12)
        # using the metallic sheath thickness t_s in mm.
        thickness_term = self._eddy_thickness_term(beta1)

        # I assemble the base eddy-current loss factor before the F correction.
        lambda1pp = (Rs / Rac) * (gs * lambda0 * (1.0 + delta1 + delta2) + thickness_term)

        # I apply the F factor from IEC 60287-1-1:2023 Section 2.3.5 when the
        # caller requests it. TB 880 Guidance Point 31 extends F to all
        # conductor designs, not only Milliken. The physical reason is that
        # circulating sheath currents reduce the magnetic field that drives
        # eddy currents, so lambda1'' is smaller when both loss mechanisms
        # coexist. F is always between 0 and 1.
        #
        # The caller derives apply_F from three central flags:
        # include_sheath_circulating_losses,
        # include_sheath_eddy_losses, and
        # include_F_factor_for_eddy_reduction.
        # If any of these is False, apply_F will be False.
        F = 1.0
        if apply_F:
            X = self.sheath_reactance_X()
            # I use M = R_s / X for trefoil formation because M = N in trefoil.
            M = Rs / X
            F = (4.0 * M**4 + (2.0 * M)**2) / (4.0 * (M**2 + 1.0)**2)
            lambda1pp = lambda1pp * F

        return {
            "lambda1_doubleprime": lambda1pp,
            "F": F,
            "beta1": beta1,
            "gs": gs,
            "lambda0": lambda0,
            "delta1": delta1,
            "delta2": delta2,
            "thickness_term": thickness_term,
        }


    def dielectric_loss_W_per_m(self):
        # I use IEC 60287-1-1 Eq. (15) for capacitance, with d_c taken as the
        # diameter over conductor screen and D_i taken as the diameter over
        # insulation; semiconducting layers are excluded from the dielectric-only
        # geometry selection for this benchmark.
        Di = self.d_ins
        d = self.d_inner_semicon
        if Di <= d:
            raise Exception("Need Di > d")
        C = 2.0 * math.pi * CASE.constants.eps0_f_per_m * self.eps_r / math.log(Di / d)

        # I then apply IEC 60287-1-1 Eq. (14), using U_0 = U_LL / sqrt(3), because
        # dielectric losses in a single-core cable depend on phase-to-earth voltage.
        U0 = self.U_LL_V / math.sqrt(3.0)
        omega = 2.0 * math.pi * self.f
        return omega * C * U0**2 * self.tan_delta

    def calculate_losses(self, T_core_C, T_screen_C, area_mode="analytical"):
        """Single IEC 60287 loss chain used by the analytical and FEM workflows."""
        # I begin with the conductor resistance chain R_dc -> R_ac because conductor
        # Joule loss W_c = I^2 * R_ac is the dominant source term in the benchmark.
        Rac = self.Rac_at_temp(T_core_C)
        Rdc = self.Rdc_at_temp(T_core_C)
        Rs = self.Rs_at_temp(T_screen_C, area_mode=area_mode)

        # I evaluate dielectric loss independently from current because IEC 60287-1-1
        # Eq. (14) depends on voltage, capacitance, and tan(delta), not on load current.
        Wd = self.dielectric_loss_W_per_m()
        Wc = self.I**2 * Rac

        # I compute sheath losses through lambda1_prime and lambda1_doubleprime
        # because IEC 60287-1-1 defines W_s = lambda1 * W_c for the metallic
        # sheath, where lambda1 is assembled from circulating and eddy terms.
        X = self.sheath_reactance_X()
        lambda1_prime = self._lambda1_prime(Rs, Rac, X)

        # I derive the apply_F decision from three central flags on CASE.assumptions:
        #
        # 1) include_sheath_circulating_losses must be True because F represents
        #    the physical suppression of eddy currents by circulating currents.
        #    If circulating losses are off, there is nothing to suppress eddies.
        #
        # 2) include_sheath_eddy_losses must be True because F modifies lambda1''
        #    and lambda1'' is only relevant if eddy losses are included.
        #
        # 3) include_F_factor_for_eddy_reduction must be True because the user
        #    may want to deliberately omit F to study its effect on the rating.
        #    This corresponds to the IEC simplified approach for non-Milliken
        #    conductors, where eddy losses are either neglected entirely or
        #    computed without the F reduction.
        #
        # Reference: IEC 60287-1-1:2023 Section 2.3.5 for the F factor formula.
        # Reference: TB 880 Guidance Point 31 for extending F to all conductor types.
        apply_F = (
            CASE.assumptions.include_sheath_circulating_losses
            and CASE.assumptions.include_sheath_eddy_losses
            and CASE.assumptions.include_F_factor_for_eddy_reduction
        )
        eddy_terms = self._lambda1_doubleprime(Rs, Rac, apply_F=apply_F)
        lambda1_doubleprime = eddy_terms["lambda1_doubleprime"]

        # I assemble lambda1 only from the centralized CASE flags so that
        # sheath circulating and eddy losses can be toggled independently
        # for parametric studies of simplification effects.
        lambda1 = 0.0
        if CASE.assumptions.include_sheath_circulating_losses:
            lambda1 = lambda1 + lambda1_prime
        if CASE.assumptions.include_sheath_eddy_losses:
            lambda1 = lambda1 + lambda1_doubleprime

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
            "beta1": eddy_terms["beta1"],
            "gs": eddy_terms["gs"],
            "lambda0": eddy_terms["lambda0"],
            "delta1": eddy_terms["delta1"],
            "delta2": eddy_terms["delta2"],
            "thickness_term": eddy_terms["thickness_term"],
            "F": eddy_terms["F"],
            "Rac": Rac,
            "Rdc": Rdc,
            "Rs": Rs,
        }

    def q_core(self, T_core_C, T_screen_C, area_mode="fem"):
        # I convert conductor loss from W/m to W/m^3 using the selected area mode,
        # defaulting to the ANSYS-loaded FEM body area for volumetric heating.
        losses = self.calculate_losses(T_core_C, T_screen_C, area_mode=area_mode)
        return losses["core"] / self._get_area("conductor", area_mode)

    def q_screen(self, T_core_C, T_screen_C, area_mode="fem"):
        # I convert sheath loss from W/m to W/m^3 using the selected screen area so
        # ANSYS applies heating to the exact body for the active workflow.
        losses = self.calculate_losses(T_core_C, T_screen_C, area_mode=area_mode)
        return losses["screen"] / self._get_area("screen", area_mode)

    def q_innerins(self, T_core_C, T_screen_C, area_mode="fem"):
        # I convert dielectric loss from W/m to W/m^3 using the selected InnerIns
        # area, which is only available in FEM mode for this benchmark workflow.
        losses = self.calculate_losses(T_core_C, T_screen_C, area_mode=area_mode)
        return losses["dielectric"] / self._get_area("innerins", area_mode)


def create_case0_cables(I_rms_A=None):
    """Return three symmetric trefoil Cable objects built from the centralized CASE data."""
    if I_rms_A is None:
        I_rms_A = CASE.benchmark.i_final_a

    d_core = CASE.geometry.d_cond_m
    d_inner_semicon = CASE.geometry.d_inner_semicon_m
    d_ins = CASE.geometry.d_ins_m
    d_outer_semicon = CASE.geometry.d_outer_semicon_m
    d_screen_mean = CASE.geometry.d_screen_mean_m
    d_screen_out = CASE.geometry.d_screen_out_m
    d_oversheath = CASE.geometry.d_oversheath_m

    A_cond_elec = CASE.geometry.a_cond_elec_m2
    A_screen_elec = CASE.geometry.a_screen_elec_m2
    A_cond_FE_geom = CASE.geometry.a_cond_geom_m2
    A_innerins_FE_geom = CASE.geometry.a_innerins_geom_m2
    A_screen_FE_geom = CASE.geometry.a_screen_geom_m2

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
            A_cond_FE_geom=A_cond_FE_geom,
            A_innerins_FE_geom=A_innerins_FE_geom,
            A_screen_FE_geom=A_screen_FE_geom,
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
        )
    return cables
