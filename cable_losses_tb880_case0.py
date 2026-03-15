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


class Cable:
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

        # Electrical properties
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
        return 1.0 + alpha * (T_C - 20.0)

    def Rdc_at_temp(self, T_cond_C):
        """Conductor DC resistance [ohm/m] at temperature from stored 20°C benchmark input."""
        return self.r_cond_dc_20_ohm_per_m * self._resistance_ratio(self.alpha_cond, T_cond_C)

    def Rac_at_temp(self, T_cond_C):
        Rdc = self.Rdc_at_temp(T_cond_C)

        xs2 = (8.0 * math.pi * self.f * 1e-7 * self.ks) / Rdc
        xs = math.sqrt(xs2)
        if xs <= 2.8:
            ys = xs**4 / (192.0 + 0.8 * xs**4)
        elif xs <= 3.8:
            ys = -0.136 - 0.0177 * xs + 0.0563 * xs**2
        else:
            ys = 0.354 * xs - 0.733

        xp2 = (8.0 * math.pi * self.f * 1e-7 * self.kp) / Rdc
        xp = math.sqrt(xp2)
        g = xp**4 / (192.0 + 0.8 * xp**4)
        s = self.d_oversheath
        r = (self.d_cond / s) ** 2
        yp = g * r * (0.312 * r + 1.18 / (g + 0.27))

        return Rdc * (1.0 + ys + yp)

    def Rs_at_temp(self, T_screen_C):
        rho_screen = self.rho_screen_20 * self._resistance_ratio(self.alpha_screen, T_screen_C)
        return rho_screen / self.A_screen_elec

    def sheath_resistance_at_temp(self, T_screen_C):
        return self.Rs_at_temp(T_screen_C)

    def sheath_reactance_X(self):
        s = self.d_oversheath
        d = self.d_screen_mean
        if (2.0 * s) <= d:
            raise ValueError("Need 2*s > d for ln(2s/d)")
        omega = 2.0 * math.pi * self.f
        return 2.0 * omega * 1e-7 * math.log((2.0 * s) / d)

    def _lambda1_prime(self, Rs, Rac, X):
        return (Rs / Rac) * (1.0 / (1.0 + (Rs / X) ** 2))

    def _lambda1_doubleprime(self, Rs, Rac):
        omega = 2.0 * math.pi * self.f
        m = omega * 1e-7 / Rs

        d = self.d_screen_mean
        s = self.d_oversheath

        lambda0 = 3.0 * (m**2 / (1.0 + m**2)) * (d / (2.0 * s)) ** 2
        delta1 = (1.14 * m**2.45 + 0.33) * (d / (2.0 * s)) ** (0.92 * m + 1.66)
        delta2 = 0.0
        gs = 1.0
        thickness_term = 0.0

        lambda1pp = (Rs / Rac) * (gs * lambda0 * (1.0 + delta1 + delta2) + thickness_term)

        if self.bonding == "solid":
            X = self.sheath_reactance_X()
            M = Rs / X
            F = (4.0 * M**4 + (2.0 * M)**2) / (4.0 * (M**2 + 1.0)**2)
            lambda1pp *= F

        return lambda1pp

    def _include_eddy_for_current_setup(self):
        mode = self.sheath_eddy_policy
        if mode == "exclude":
            return False
        if mode == "include_for_solid":
            return self.bonding == "solid"
        if mode == "include_for_all":
            return True

        # auto -> repository default interpretation from CASE assumptions.
        if self.bonding == "solid":
            return bool(CASE.assumptions.solid_bonded_include_eddy_default)
        return True

    def dielectric_loss_W_per_m(self):
        Di = self.d_outer_semicon
        d = self.d_inner_semicon
        if Di <= d:
            raise Exception("Need Di > d")
        C = 2.0 * math.pi * CASE.constants.eps0_f_per_m * self.eps_r / math.log(Di / d)

        U0 = self.U_LL_V / math.sqrt(3.0)
        omega = 2.0 * math.pi * self.f
        return omega * C * U0**2 * self.tan_delta

    def calculate_losses(self, T_core_C, T_screen_C):
        """Single authoritative loss path for analytical and FEM coupling workflows."""
        Rac = self.Rac_at_temp(T_core_C)
        Rdc = self.Rdc_at_temp(T_core_C)
        Rs = self.Rs_at_temp(T_screen_C)

        Wd = self.dielectric_loss_W_per_m()
        Wc = self.I**2 * Rac

        X = self.sheath_reactance_X()
        lambda1_prime = self._lambda1_prime(Rs, Rac, X) if self.bonding == "solid" else 0.0
        lambda1_doubleprime = self._lambda1_doubleprime(Rs, Rac) if self._include_eddy_for_current_setup() else 0.0

        if self.bonding == "solid":
            lambda1 = lambda1_prime + lambda1_doubleprime
        else:
            lambda1 = lambda1_doubleprime

        Ws = lambda1 * Wc

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
        losses = self.calculate_losses(T_core_C, T_screen_C)
        return losses["core"] / self.A_cond_FE

    def q_screen(self, T_core_C, T_screen_C):
        losses = self.calculate_losses(T_core_C, T_screen_C)
        return losses["screen"] / self.A_screen_FE

    def q_innerins(self, T_core_C, T_screen_C):
        losses = self.calculate_losses(T_core_C, T_screen_C)
        return losses["dielectric"] / self.A_innerins_FE


def create_case0_cables(I_rms_A=None, bonding=None, sheath_eddy_policy="auto"):
    """Return three Cable objects (C01/C02/C03) with centralized TB880 Case #0 data."""
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


cables = create_case0_cables()
