# ============================================================
# cable_losses_tb880_case0.py (REVISED)
# ------------------------------------------------------------
# IEC 60287-1-1 style losses for TB 880 - Case 0 (TREFOIL)
#
# Improvements:
#  Object oriented Cable class per phase.
#  Full sheath losses: Lambda_1 strich (circulating) + Lambda_1 strich strich (eddy), with F factor for solid bonding.
#  Separate electrical area (for loss calc) and FE area (for volumetric source).
#  Stranded conductor: effective radial thermal conductivity.
#  Explicit geometry: inner/outer semicon layers (thermal resistance only, no heat).
#  All calculations follow guidance of TB880 and TB963.
# ============================================================

import math
from tb880_case0_data import CASE


def normalize_bonding_mode(bonding):
    """Normalize bonding mode to: solid, single, cross (with one legacy alias mapped to solid)."""
    mode = str(bonding).strip().lower()
    legacy = {"both": "solid"}
    normalized = legacy.get(mode, mode)
    if normalized not in ("solid", "single", "cross"):
        raise ValueError("Unsupported bonding mode '{0}'. Use: solid, single, cross".format(bonding))
    return normalized

class Cable:
    """
    Represents one phase of a cable. Contains all geometry, material,
    and electrical parameters needed to compute losses according to
    IEC 60287-1-1.
    """

    def __init__(self, cid,
                 # Geometry (all in meters)
                 d_cond,                # conductor outer diameter (physical)
                 d_inner_semicon,        # inner semicon outer diameter
                 d_ins,                  # insulation outer diameter
                 d_outer_semicon,         # outer semicon outer diameter
                 d_screen_mean,           # mean diameter of metallic screen
                 d_screen_out,            # screen outer diameter (for oversheath)
                 d_oversheath,            # oversheath outer diameter

                 # Cross-sectional areas (m^2)
                 A_cond_elec,             # electrical area of conductor (e.g., 630e-6)
                 A_screen_elec,            # electrical area of screen

                 # FE areas (m^2)  for converting W/m -> W/m^3
                 A_cond_FE,
                 A_innerins_FE,
                 A_screen_FE,

                 # Material properties
                 rho_cond_20,             # conducter resistivity at 20 C (ohm * m)
                 alpha_cond,               # temperature coefficient [1/K]
                 ks, kp,                   # skin/proximity constants
                 rho_screen_20,            # screen resistivity at 20 C [Ohm*m]
                 alpha_screen,              # screen temperature coeff. [1/K]
                 eps_r,                    # relative permittivity of insulation
                 tan_delta,                 # dielectric loss factor

                 # Installation / operating
                 f_hz,                      # frequency [Hz]
                 U_LL_V,                    # line line voltage [V]
                 I_rms_A,                    # current [A] (RMS)

                 # Bonding scheme
                 bonding="solid",            # "solid", "single", "cross"
                 # For cross bonding, you may need minor section lengths; simplified here.
                 ):
        self.cid = cid
        self.I = float(I_rms_A)
        self.f = float(f_hz)
        self.U_LL_V = float(U_LL_V)

        # ----- Geometry -------------------------------------------------
        # All diameters in meters
        self.d_cond = d_cond
        self.d_inner_semicon = d_inner_semicon
        self.d_ins = d_ins
        self.d_outer_semicon = d_outer_semicon
        self.d_screen_mean = d_screen_mean
        self.d_screen_out = d_screen_out
        self.d_oversheath = d_oversheath

        # Areas
        self.A_cond_elec = A_cond_elec
        self.A_screen_elec = A_screen_elec
        self.A_cond_FE = A_cond_FE
        self.A_innerins_FE = A_innerins_FE
        self.A_screen_FE = A_screen_FE

        # ----- Material properties ---------------------------------------
        self.rho_cond_20 = rho_cond_20
        self.alpha_cond = alpha_cond
        self.ks = ks
        self.kp = kp
        self.rho_screen_20 = rho_screen_20
        self.alpha_screen = alpha_screen
        self.eps_r = eps_r
        self.tan_delta = tan_delta

        # Derived (will be computed later)
        self.Rac = None          # AC resistance at current temperature
        self.Rs  = None          # screen resistance at current temperature
        self.lambda1 = None      # total sheath loss factor

        # Bonding scheme
        self.bonding = normalize_bonding_mode(bonding)

    # ------------------------------------------------------------------
    # Helper: resistance at temperature
    # ------------------------------------------------------------------
    def _resistance(self, rho20, alpha, T_C):
        return rho20 * (1.0 + alpha * (T_C - 20.0))

    # ------------------------------------------------------------------
    # Conductor AC resistance (IEC 60287-1-1 Chapter 2.1)
    # ------------------------------------------------------------------
    def Rac_at_temp(self, T_cond_C):
        # DC resistance at temperature
        rho_cond = self._resistance(self.rho_cond_20, self.alpha_cond, T_cond_C)
        # Convert resistivity to resistance per meter using electrical area
        Rdc = rho_cond / self.A_cond_elec

        # Skin effect factor ys
        xs2 = (8.0 * math.pi * self.f * 1e-7 * self.ks) / Rdc
        xs = math.sqrt(xs2)
        if xs <= 2.8:
            ys = xs**4 / (192.0 + 0.8 * xs**4)
        elif xs <= 3.8:
            ys = -0.136 - 0.0177 * xs + 0.0563 * xs**2
        else:
            ys = 0.354 * xs - 0.733

        # Proximity effect factor yp
        xp2 = (8.0 * math.pi * self.f * 1e-7 * self.kp) / Rdc
        xp = math.sqrt(xp2)
        g = xp**4 / (192.0 + 0.8 * xp**4)
        # For trefoil, distance s is approx oversheath diameter (touching)
        s = self.d_oversheath   # center to center distance (touching)
        r = (self.d_cond / s) ** 2
        yp = g * r * (0.312 * r + 1.18 / (g + 0.27))

        Rac = Rdc * (1.0 + ys + yp)
        return Rac

    # ------------------------------------------------------------------
    # Sheath resistance at temperature
    # ------------------------------------------------------------------
    def Rs_at_temp(self, T_screen_C):
        rho_screen = self._resistance(self.rho_screen_20, self.alpha_screen, T_screen_C)
        return rho_screen / self.A_screen_elec

    # Backward-compatible alias expected by legacy analytical solver code.
    def sheath_resistance_at_temp(self, T_screen_C):
        return self.Rs_at_temp(T_screen_C)

    # ------------------------------------------------------------------
    # Sheath reactance X (for circulating currents)
    # ------------------------------------------------------------------
    def sheath_reactance_X(self):
        # For trefoil, use mean screen diameter d and spacing s = outer diameter
        s = self.d_oversheath
        d = self.d_screen_mean
        if (2.0 * s) <= d:
            raise ValueError("Need 2*s > d for ln(2s/d)")
        omega = 2.0 * math.pi * self.f
        return 2.0 * omega * 1e-7 * math.log((2.0 * s) / d)

    # ------------------------------------------------------------------
    # Sheath loss factors: Lambda_1 strich (circulating) and Lambda_1 strich strich (eddy)
    # ------------------------------------------------------------------
    def _lambda1_prime(self, Rs, Rac, X):
        # Eq. from IEC 60287-1-1 Chapter 2.3.1
        return (Rs / Rac) * (1.0 / (1.0 + (Rs / X) ** 2))

    def _lambda1_doubleprime(self, Rs, Rac, T_screen_C):
        # Eddy current loss factor for isolated sheaths (IEC Chapter 2.3.6.1)
        # For trefoil formation.
        omega = 2.0 * math.pi * self.f
        m = omega * 1e-7 / Rs   # (using m = omega/Rs * 1e-7), omega  = 2 * pi * f

        d = self.d_screen_mean
        s = self.d_oversheath

        # Base term Lambda_0
        lambda0 = 3.0 * (m**2 / (1.0 + m**2)) * (d / (2.0 * s)) ** 2

        # Corrections Delta_1, Delta_2 (Delta_2 = 0 for trefoil)
        # Coefficients from IEC
        delta1 = (1.14 * m**2.45 + 0.33) * (d / (2.0 * s)) ** (0.92 * m + 1.66)
        delta2 = 0.0

        # Correction factor gs (accounts for sheath thickness)
        # We use simplified gs = 1 (for thin sheaths). For thick Al sheaths,
        # proper formula from IEC could be added, but omitted here for brevity.
        gs = 1.0

        # Term due to sheath thickness (Beta_1 t)^4 / (12*10^12)  usually negligible for thin sheaths
        # For Al sheaths one may need it; here we set zero.
        thickness_term = 0.0

        # Lambda_1 strich strich = (Rs/Rac) * [gs Lambda_0 (1+Delta_1+Delta_2) + thickness_term]
        lambda1pp = (Rs / Rac) * (gs * lambda0 * (1.0 + delta1 + delta2) + thickness_term)

        # For solid bonding, eddy currents are reduced by circulating currents.
        # Apply factor F (IEC Chapter 2.3.5) only if bonding is solid.
        if self.bonding == "solid":
            # Factor F for trefoil: M = N = Rs/X
            X = self.sheath_reactance_X()
            M = Rs / X
            F = (4.0 * M**4 + (2.0 * M)**2) / (4.0 * (M**2 + 1.0)**2)   # simplified when M=N
            lambda1pp *= F

        return lambda1pp

    # ------------------------------------------------------------------
    # Dielectric losses (IEC Chapter 2.2)
    # ------------------------------------------------------------------
    def dielectric_loss_W_per_m(self):
        # Capacitance per meter
        # Inner electrode: conductor screen (d_inner_semicon), outer: insulation screen (d_outer_semicon)
        Di = self.d_outer_semicon
        d  = self.d_inner_semicon
        if Di <= d:
            raise Exception("Need Di > d")
        C = 2.0 * math.pi * CASE.constants.eps0_f_per_m * self.eps_r / math.log(Di / d)

        U0 = self.U_LL_V / math.sqrt(3.0)
        omega = 2.0 * math.pi * self.f
        Wd = omega * C * U0**2 * self.tan_delta
        return Wd

    # ------------------------------------------------------------------
    # Main loss calculation: return dict of losses in W/m
    # ------------------------------------------------------------------
    def calculate_losses(self, T_core_C, T_screen_C):
        """
        Compute all heat sources for this cable at given core and screen temperatures.
        Returns: dict with keys 'core', 'screen', 'dielectric' in W/m.
        """
        # Conductor AC resistance at T_core
        Rac = self.Rac_at_temp(T_core_C)
        self.Rac = Rac

        # Screen resistance at T_screen
        Rs = self.Rs_at_temp(T_screen_C)
        self.Rs = Rs

        # Dielectric loss (constant, temperature independent for XLPE)
        Wd = self.dielectric_loss_W_per_m()

        # Conductor loss
        Wc = self.I**2 * Rac

        # Sheath loss
        X = self.sheath_reactance_X()
        lambda1_prime = self._lambda1_prime(Rs, Rac, X)

        if self.bonding == "solid":   # solid-bonded at both ends
            lambda1_doubleprime = self._lambda1_doubleprime(Rs, Rac, T_screen_C)
            lambda1 = lambda1_prime + lambda1_doubleprime
        else:                                   # single-point or cross-bonded
            lambda1 = self._lambda1_doubleprime(Rs, Rac, T_screen_C)
            lambda1_prime = 0.0   # not used, but keep for clarity

        self.lambda1 = lambda1
        Ws = lambda1 * Wc

        # Store values (optional)
        self.Wc = Wc
        self.Ws = Ws
        self.Wd = Wd

        return {
            "core": Wc,
            "screen": Ws,
            "dielectric": Wd
        }

    # ------------------------------------------------------------------
    # Convert W/m to volumetric heat source (W/m^3) for FE
    # ------------------------------------------------------------------
    def q_core(self, T_core_C, T_screen_C):
        losses = self.calculate_losses(T_core_C, T_screen_C)
        return losses["core"] / self.A_cond_FE

    def q_screen(self, T_core_C, T_screen_C):
        losses = self.calculate_losses(T_core_C, T_screen_C)
        return losses["screen"] / self.A_screen_FE

    def q_innerins(self, T_core_C, T_screen_C):
        losses = self.calculate_losses(T_core_C, T_screen_C)
        return losses["dielectric"] / self.A_innerins_FE


# ----------------------------------------------------------------------
# Helper function to instantiate three cables with common TB880 Case 0 data
# ----------------------------------------------------------------------
def create_case0_cables(I_rms_A=None, bonding=None):
    """
    Returns a dict of three Cable objects for C01, C02, C03 with TB880 Case 0 parameters.
    """
    if I_rms_A is None:
        I_rms_A = CASE.benchmark.i_final_a
    if bonding is None:
        bonding = CASE.installation.bonding
    bonding = normalize_bonding_mode(bonding)

    # Geometry and areas from centralized TB 880 Case 0 data.
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
    # FE area convention is centralized in CASE.assumptions.innerins_area_convention.
    A_innerins_FE = CASE.geometry.a_innerins_geom_m2
    A_screen_FE = CASE.geometry.a_screen_geom_m2

    # Material and operating data from centralized TB 880 Case 0 data.
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
            rho_cond_20=rho_cond_20,
            alpha_cond=alpha_cond,
            ks=ks, kp=kp,
            rho_screen_20=rho_screen_20,
            alpha_screen=alpha_screen,
            eps_r=eps_r,
            tan_delta=tan_delta,
            f_hz=f_hz,
            U_LL_V=U_LL_V,
            I_rms_A=I_rms_A,
            bonding=bonding
        )
    return cables


# Backward-compatible module-level dictionary expected by legacy callers.
cables = create_case0_cables()
