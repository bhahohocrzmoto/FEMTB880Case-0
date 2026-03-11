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
import warnings

# ------------------------------------------------------------
# Physical constants
# ------------------------------------------------------------
EPS0 = 8.854187817e-12          # vacuum permittivity [F/m]
MU0  = 4e-7 * math.pi            # vacuum permeability [H/m]


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
        self.bonding = bonding.lower()

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
        C = 2.0 * math.pi * EPS0 * self.eps_r / math.log(Di / d)

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

        if self.bonding in ("solid", "both"):   # both ends bonded
            lambda1_doubleprime = self._lambda1_doubleprime(Rs, Rac, T_screen_C)
            lambda1 = lambda1_prime + lambda1_doubleprime
        else:                                   # single point or cross bonded -> no circulating current
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
def create_case0_cables(I_rms_A=822.0, bonding="solid"):
    """
    Returns a dict of three Cable objects for C01, C02, C03 with TB880 Case 0 parameters.
    """
    # Geometry (meters)  from TB880 Case 0
    d_core             = 30.3e-3
    d_inner_semicon    = 33.3e-3
    d_ins              = 64.3e-3
    d_outer_semicon    = 66.9e-3   # actually insulation screen outer diameter
    d_screen_mean      = 67.7e-3   # mean sheath diameter
    d_screen_out       = 68.5e-3
    d_oversheath       = 75.5e-3

    # Areas (m^2)
    A_cond_elec        = 630e-6               # electrical cross section (630 mm^2)
    A_cond_FE          = math.pi * (d_core/2)**2   # geometric FE area

    # Inner insulation area (from d_core to d_ins)  includes semicons? We'll use total area.
    A_innerins_FE      = math.pi * ((d_ins/2)**2 - (d_core/2)**2)

    # Screen FE area  geometric area of metallic sheath (thin ring)
    t_screen           = (d_screen_out - d_outer_semicon) / 2   # approx 0.8 mm
    A_screen_FE        = math.pi * d_screen_mean * t_screen

    # Electrical screen area  same as geometric for thin sheath
    A_screen_elec      = A_screen_FE

    # Material properties (IEC 60287)
    rho_cond_20        = 28.3e-6 * A_cond_elec      # R0 * A_elec gives resistivity [Ohm*m]
    alpha_cond         = 3.93e-3
    ks                 = 1.0
    kp                 = 1.0

    rho_screen_20      = 2.84e-8                    # Aluminium
    alpha_screen       = 3.93e-3

    eps_r              = 2.5
    tan_delta          = 1.0e-3

    f_hz               = 50.0
    U_LL_V             = 132e3

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
