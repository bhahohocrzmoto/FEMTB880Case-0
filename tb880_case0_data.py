# tb880_case0_data.py
# Centralized, object-based storage for CIGRE TB 880 Case #0
# -----------------------------------------------------------
# This module contains ONLY case-specific data and benchmark values.
# No solver logic should live here.
#
# References available from uploaded sources:
# - TB 880 (2022), Introductory Case study 0, input data and benchmark results
# - TB 963 (2025), Case study 1 uses TB 880 Case 0 cable parameters
# - IEC 60287-1-1:2023, Table 1 and Table 2
# - IEC 60287-2-1:2023, T1/T3/T4 formula framework
#
# Important design rule:
#   Keep physical input data separate from benchmark results.
#   Solver scripts should import this object and never hard-code constants.

from dataclasses import dataclass, field
import math


# ----------------------------------------------------------------------
# Reference metadata
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class Ref:
    doc: str
    section: str
    page: str
    note: str = ""


# ----------------------------------------------------------------------
# Fundamental constants
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class PhysicalConstants:
    eps0_f_per_m: float = 8.854187817e-12
    mu0_h_per_m: float = 4e-7 * math.pi


# ----------------------------------------------------------------------
# Cable geometry for TB 880 Case 0
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class Geometry:
    # Diameters [m]
    d_cond_m: float
    d_inner_semicon_m: float
    d_ins_m: float
    d_outer_semicon_m: float
    d_screen_mean_m: float
    d_screen_out_m: float
    d_oversheath_m: float

    # Electrical areas [m²]
    a_cond_elec_m2: float
    a_screen_elec_m2: float

    @property
    def screen_thickness_m(self) -> float:
        return 0.5 * (self.d_screen_out_m - self.d_outer_semicon_m)

    @property
    def a_cond_geom_m2(self) -> float:
        return math.pi * (self.d_cond_m / 2.0) ** 2

    @property
    def a_innerins_geom_m2(self) -> float:
        # Convention for this benchmark workflow: annulus from conductor OD to insulation OD.
        # This geometric FE area may include semicon sublayers when they are part of one InnerIns FE body.
        return math.pi * ((self.d_ins_m / 2.0) ** 2 - (self.d_cond_m / 2.0) ** 2)

    @property
    def a_screen_geom_m2(self) -> float:
        # Thin-ring approximation
        return math.pi * self.d_screen_mean_m * self.screen_thickness_m

    @property
    def spacing_touching_trefoil_m(self) -> float:
        # In touching trefoil, centre-to-centre spacing is cable outer diameter.
        return self.d_oversheath_m


# ----------------------------------------------------------------------
# Material and electrical parameters
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class MaterialElectricalData:
    # Benchmark-aligned conductor DC resistance at 20°C [ohm/m]
    # This is a primary Case #0 input for IEC/TB880 benchmark traceability.
    r_cond_dc_20_ohm_per_m: float

    # True material resistivities at 20°C [ohm·m]
    rho_cond_20_ohm_m: float
    rho_screen_20_ohm_m: float

    # Temperature coefficients at 20°C [1/K]
    alpha_cond_20_per_k: float
    alpha_screen_20_per_k: float

    # IEC skin/proximity coefficients
    ks: float
    kp: float

    # Dielectric data
    eps_r: float
    tan_delta: float

    # Thermal resistivities [K·m/W]
    rho_semicon_k_m_per_w: float
    rho_ins_k_m_per_w: float
    rho_oversheath_k_m_per_w: float

# ----------------------------------------------------------------------
# Installation / operating data
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class Installation:
    frequency_hz: float
    u_ll_v: float
    ambient_temp_c: float
    burial_depth_to_trefoil_center_m: float
    soil_thermal_resistivity_k_m_per_w: float
    bonding: str

    @property
    def u0_v(self) -> float:
        return self.u_ll_v / math.sqrt(3.0)


# ----------------------------------------------------------------------
# Benchmark values from TB 880 Case 0-1
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class Benchmark:
    # Thermal resistances [K·m/W]
    t1_k_m_per_w: float
    t2_k_m_per_w: float
    t3_k_m_per_w: float
    t4_k_m_per_w: float

    # Benchmark losses and rating
    wd_w_per_m: float
    i_final_a: float
    wc_final_w_per_m: float
    ws_final_w_per_m: float

    # Final temperatures
    theta_core_final_c: float
    theta_screen_final_c: float
    theta_oversheath_final_c: float

    # Useful first-iteration benchmark quantities
    theta_screen_init_c: float
    r_ac_90c_ohm_per_m: float
    r_screen_80c_ohm_per_m: float
    lambda1_circulating_first_iter: float


# ----------------------------------------------------------------------
# Source map
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class SourceMap:
    geometry: Ref
    material_resistivity: Ref
    skin_proximity: Ref
    iteration_guidance: Ref
    benchmark_t1_t3: Ref
    benchmark_t4_wd: Ref
    benchmark_final: Ref
    fem_alignment: Ref


# ----------------------------------------------------------------------
# Explicit model assumptions / conventions used by this benchmark setup
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class ModelAssumptions:
    # The InnerIns FE body is treated as one annulus from conductor OD to insulation OD.
    # This area is used to convert dielectric W/m to W/m^3 for NS_C0x_InnerIns.
    innerins_area_convention: str

    # Informational traceability flag: in this case setup, screen FE and electrical
    # areas are intentionally treated as equivalent for load conversion and resistance.
    screen_fe_area_equals_electrical_area_assumed: bool

    # Touching trefoil (centre spacing = cable outer diameter) is the default geometry.
    touching_trefoil_spacing_assumed: bool

    # Repository default interpretation for sheath losses in TB 880 Case #0-1 regression.
    sheath_loss_interpretation: str

    # In the default benchmark-faithful path, solid-bonded eddy sheath losses are neglected.
    solid_bonded_include_eddy_default: bool


# ----------------------------------------------------------------------
# Master case object
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class TB880Case0:
    case_id: str
    title: str
    constants: PhysicalConstants
    geometry: Geometry
    material: MaterialElectricalData
    installation: Installation
    benchmark: Benchmark
    sources: SourceMap
    assumptions: ModelAssumptions


TB880_CASE_0 = TB880Case0(
    case_id="TB880_CASE_0",
    title="CIGRE TB 880 Introductory Case Study #0 - 132 kV single-core XLPE cable, touching trefoil, directly buried",
    constants=PhysicalConstants(),
    geometry=Geometry(
        # TB 880 Case 0 cable geometry
        # Figure 25 / Figure 26, Case #0 input data
        d_cond_m=30.3e-3,
        d_inner_semicon_m=33.3e-3,
        d_ins_m=64.3e-3,
        d_outer_semicon_m=66.9e-3,
        d_screen_mean_m=67.7e-3,
        d_screen_out_m=68.5e-3,
        d_oversheath_m=75.5e-3,

        # Electrical cross-sections
        a_cond_elec_m2=630e-6,
        # TB 963 states the TB 880 Case 0 aluminium laminated sheath cross-section as 170 mm².
        a_screen_elec_m2=170e-6,
    ),
    material=MaterialElectricalData(
        # IEC 60228 / TB 880-aligned benchmark primary input for 630 mm² copper conductor.
        r_cond_dc_20_ohm_per_m=2.83e-5,

        # IEC 60287-1-1 Table 1
        rho_cond_20_ohm_m=1.7241e-8,   # Copper
        rho_screen_20_ohm_m=2.84e-8,   # Aluminium

        alpha_cond_20_per_k=3.93e-3,
        alpha_screen_20_per_k=4.03e-3,

        # IEC 60287-1-1 Table 2 for round stranded conductor with extruded insulation
        ks=1.0,
        kp=1.0,

        # TB 880 Case 0 / IEC assumptions
        eps_r=2.5,
        tan_delta=1.0e-3,

        # Thermal resistivities
        rho_semicon_k_m_per_w=2.5,
        rho_ins_k_m_per_w=3.5,
        rho_oversheath_k_m_per_w=3.5,
    ),
    installation=Installation(
        frequency_hz=50.0,
        u_ll_v=132e3,
        ambient_temp_c=20.0,
        burial_depth_to_trefoil_center_m=1.0,
        soil_thermal_resistivity_k_m_per_w=1.0,
        bonding="solid",   # both ends bonded, base Case #0-1
    ),
    benchmark=Benchmark(
        # TB 880 Case #0-1
        t1_k_m_per_w=0.4198714890,
        t2_k_m_per_w=0.0,
        t3_k_m_per_w=0.08671937479,
        t4_k_m_per_w=1.594692892,

        wd_w_per_m=0.3851382172,
        i_final_a=821.7763334392,
        wc_final_w_per_m=26.6895,
        ws_final_w_per_m=7.8442,

        theta_core_final_c=90.0,
        theta_screen_final_c=78.7130,
        theta_oversheath_final_c=75.68,

        theta_screen_init_c=80.0,
        r_ac_90c_ohm_per_m=3.952152638e-5,
        r_screen_80c_ohm_per_m=2.072723957e-4,
        lambda1_circulating_first_iter=0.2928142510,
    ),
    sources=SourceMap(
        geometry=Ref(
            doc="CIGRE TB 880 (2022)",
            section="Case #0 input data, Figure 25 and Figure 26",
            page="p. 63",
            note="Layer diameters for Case 0 cable"
        ),
        material_resistivity=Ref(
            doc="IEC 60287-1-1:2023",
            section="Table 1",
            page="p. 36",
            note="Copper and aluminium resistivity and temperature coefficients"
        ),
        skin_proximity=Ref(
            doc="IEC 60287-1-1:2023",
            section="Table 2",
            page="p. 36",
            note="ks=1 and kp=1 for round stranded conductor with extruded insulation"
        ),
        iteration_guidance=Ref(
            doc="CIGRE TB 880 (2022)",
            section="Case #0 guidelines for iterative calculation",
            page="pp. 66-68",
            note="Initial screen temperature and convergence guidance"
        ),
        benchmark_t1_t3=Ref(
            doc="CIGRE TB 880 (2022)",
            section="Case #0-1 results independent of temperature",
            page="pp. 70-71",
            note="T1 layer split and T3 including 1.6 trefoil factor"
        ),
        benchmark_t4_wd=Ref(
            doc="CIGRE TB 880 (2022)",
            section="Case #0-1 external thermal resistance and dielectric losses",
            page="pp. 71-72",
            note="T4 and Wd"
        ),
        benchmark_final=Ref(
            doc="CIGRE TB 880 (2022)",
            section="Case #0-1 next iterations / final values",
            page="pp. 74-76",
            note="Final current rating and ladder network values"
        ),
        fem_alignment=Ref(
            doc="CIGRE TB 963 (2025)",
            section="Case study 1 input data",
            page="pp. 103-104",
            note="TB 963 explicitly reuses TB 880 Case 0 cable parameters for FEM benchmarking"
        ),
    ),
    assumptions=ModelAssumptions(
        innerins_area_convention="conductor_od_to_insulation_od",
        screen_fe_area_equals_electrical_area_assumed=True,
        touching_trefoil_spacing_assumed=True,
        sheath_loss_interpretation="tb880_case0_iec_base",
        solid_bonded_include_eddy_default=False,
    ),
)

# Convenience alias for scripts expecting a short case object name.
CASE = TB880_CASE_0
