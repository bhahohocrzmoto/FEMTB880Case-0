# -*- coding: utf-8 -*-
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

import math


# ----------------------------------------------------------------------
# Reference metadata
# ----------------------------------------------------------------------
class Ref(object):
    def __init__(self, doc, section, page, note=""):
        self.doc = doc
        self.section = section
        self.page = page
        self.note = note


# ----------------------------------------------------------------------
# Fundamental constants
# ----------------------------------------------------------------------
class PhysicalConstants(object):
    def __init__(self, eps0_f_per_m=8.854187817e-12, mu0_h_per_m=4e-7 * math.pi):
        self.eps0_f_per_m = eps0_f_per_m
        self.mu0_h_per_m = mu0_h_per_m


# ----------------------------------------------------------------------
# Cable geometry for TB 880 Case 0
# ----------------------------------------------------------------------
class Geometry(object):
    def __init__(
        self,
        d_cond_m,
        d_inner_semicon_m,
        d_ins_m,
        d_outer_semicon_m,
        d_screen_mean_m,
        d_screen_out_m,
        d_oversheath_m,
        a_cond_elec_m2,
        a_screen_elec_m2,
    ):
        # Diameters [m]
        self.d_cond_m = d_cond_m
        self.d_inner_semicon_m = d_inner_semicon_m
        self.d_ins_m = d_ins_m
        self.d_outer_semicon_m = d_outer_semicon_m
        self.d_screen_mean_m = d_screen_mean_m
        self.d_screen_out_m = d_screen_out_m
        self.d_oversheath_m = d_oversheath_m

        # Electrical areas [m2]
        self.a_cond_elec_m2 = a_cond_elec_m2
        self.a_screen_elec_m2 = a_screen_elec_m2

    @property
    def screen_thickness_m(self):
        return 0.5 * (self.d_screen_out_m - self.d_outer_semicon_m)

    @property
    def a_cond_geom_m2(self):
        return math.pi * (self.d_cond_m / 2.0) ** 2

    @property
    def a_innerins_geom_m2(self):
        # Convention for this benchmark workflow: annulus from conductor OD to insulation OD.
        # This geometric FE area may include semicon sublayers when they are part of one InnerIns FE body.
        return math.pi * ((self.d_ins_m / 2.0) ** 2 - (self.d_cond_m / 2.0) ** 2)

    @property
    def a_screen_geom_m2(self):
        # Thin-ring approximation
        return math.pi * self.d_screen_mean_m * self.screen_thickness_m

    @property
    def spacing_touching_trefoil_m(self):
        # In touching trefoil, centre-to-centre spacing is cable outer diameter.
        return self.d_oversheath_m


# ----------------------------------------------------------------------
# Material and electrical parameters
# ----------------------------------------------------------------------
class MaterialElectricalData(object):
    def __init__(
        self,
        r_cond_dc_20_ohm_per_m,
        rho_cond_20_ohm_m,
        rho_screen_20_ohm_m,
        alpha_cond_20_per_k,
        alpha_screen_20_per_k,
        ks,
        kp,
        eps_r,
        tan_delta,
        rho_semicon_k_m_per_w,
        rho_ins_k_m_per_w,
        rho_oversheath_k_m_per_w,
    ):
        # Benchmark-aligned conductor DC resistance at 20C [ohm/m]
        # This is a primary Case #0 input for IEC/TB880 benchmark traceability.
        self.r_cond_dc_20_ohm_per_m = r_cond_dc_20_ohm_per_m

        # True material resistivities at 20C [ohm*m]
        self.rho_cond_20_ohm_m = rho_cond_20_ohm_m
        self.rho_screen_20_ohm_m = rho_screen_20_ohm_m

        # Temperature coefficients at 20C [1/K]
        self.alpha_cond_20_per_k = alpha_cond_20_per_k
        self.alpha_screen_20_per_k = alpha_screen_20_per_k

        # IEC skin/proximity coefficients
        self.ks = ks
        self.kp = kp

        # Dielectric data
        self.eps_r = eps_r
        self.tan_delta = tan_delta

        # Thermal resistivities [K*m/W]
        self.rho_semicon_k_m_per_w = rho_semicon_k_m_per_w
        self.rho_ins_k_m_per_w = rho_ins_k_m_per_w
        self.rho_oversheath_k_m_per_w = rho_oversheath_k_m_per_w


# ----------------------------------------------------------------------
# Installation / operating data
# ----------------------------------------------------------------------
class Installation(object):
    def __init__(
        self,
        frequency_hz,
        u_ll_v,
        ambient_temp_c,
        burial_depth_to_trefoil_center_m,
        soil_thermal_resistivity_k_m_per_w,
        bonding,
    ):
        self.frequency_hz = frequency_hz
        self.u_ll_v = u_ll_v
        self.ambient_temp_c = ambient_temp_c
        self.burial_depth_to_trefoil_center_m = burial_depth_to_trefoil_center_m
        self.soil_thermal_resistivity_k_m_per_w = soil_thermal_resistivity_k_m_per_w
        self.bonding = bonding

    @property
    def u0_v(self):
        return self.u_ll_v / math.sqrt(3.0)


# ----------------------------------------------------------------------
# Benchmark values from TB 880 Case 0-1
# ----------------------------------------------------------------------
class Benchmark(object):
    def __init__(
        self,
        t1_k_m_per_w,
        t2_k_m_per_w,
        t3_k_m_per_w,
        t4_k_m_per_w,
        wd_w_per_m,
        i_final_a,
        wc_final_w_per_m,
        ws_final_w_per_m,
        theta_core_final_c,
        theta_screen_final_c,
        theta_oversheath_final_c,
        theta_screen_init_c,
        r_ac_90c_ohm_per_m,
        r_screen_80c_ohm_per_m,
        lambda1_circulating_first_iter,
    ):
        # Thermal resistances [K*m/W]
        self.t1_k_m_per_w = t1_k_m_per_w
        self.t2_k_m_per_w = t2_k_m_per_w
        self.t3_k_m_per_w = t3_k_m_per_w
        self.t4_k_m_per_w = t4_k_m_per_w

        # Benchmark losses and rating
        self.wd_w_per_m = wd_w_per_m
        self.i_final_a = i_final_a
        self.wc_final_w_per_m = wc_final_w_per_m
        self.ws_final_w_per_m = ws_final_w_per_m

        # Final temperatures
        self.theta_core_final_c = theta_core_final_c
        self.theta_screen_final_c = theta_screen_final_c
        self.theta_oversheath_final_c = theta_oversheath_final_c

        # Useful first-iteration benchmark quantities
        self.theta_screen_init_c = theta_screen_init_c
        self.r_ac_90c_ohm_per_m = r_ac_90c_ohm_per_m
        self.r_screen_80c_ohm_per_m = r_screen_80c_ohm_per_m
        self.lambda1_circulating_first_iter = lambda1_circulating_first_iter


# ----------------------------------------------------------------------
# Source map
# ----------------------------------------------------------------------
class SourceMap(object):
    def __init__(
        self,
        geometry,
        material_resistivity,
        skin_proximity,
        iteration_guidance,
        benchmark_t1_t3,
        benchmark_t4_wd,
        benchmark_final,
        fem_alignment,
    ):
        self.geometry = geometry
        self.material_resistivity = material_resistivity
        self.skin_proximity = skin_proximity
        self.iteration_guidance = iteration_guidance
        self.benchmark_t1_t3 = benchmark_t1_t3
        self.benchmark_t4_wd = benchmark_t4_wd
        self.benchmark_final = benchmark_final
        self.fem_alignment = fem_alignment


# ----------------------------------------------------------------------
# Explicit model assumptions / conventions used by this benchmark setup
# ----------------------------------------------------------------------
class ModelAssumptions(object):
    def __init__(
        self,
        innerins_area_convention,
        screen_fe_area_equals_electrical_area_assumed,
        touching_trefoil_spacing_assumed,
        sheath_loss_interpretation,
        solid_bonded_include_eddy_default,
    ):
        # The InnerIns FE body is treated as one annulus from conductor OD to insulation OD.
        # This area is used to convert dielectric W/m to W/m^3 for NS_C0x_InnerIns.
        self.innerins_area_convention = innerins_area_convention

        # Informational traceability flag: in this case setup, screen FE and electrical
        # areas are intentionally treated as equivalent for load conversion and resistance.
        self.screen_fe_area_equals_electrical_area_assumed = screen_fe_area_equals_electrical_area_assumed

        # Touching trefoil (centre spacing = cable outer diameter) is the default geometry.
        self.touching_trefoil_spacing_assumed = touching_trefoil_spacing_assumed

        # Repository default interpretation for sheath losses in TB 880 Case #0-1 regression.
        self.sheath_loss_interpretation = sheath_loss_interpretation

        # In the default benchmark-faithful path, solid-bonded eddy sheath losses are neglected.
        self.solid_bonded_include_eddy_default = solid_bonded_include_eddy_default


# ----------------------------------------------------------------------
# Master case object
# ----------------------------------------------------------------------
class TB880Case0(object):
    def __init__(
        self,
        case_id,
        title,
        constants,
        geometry,
        material,
        installation,
        benchmark,
        sources,
        assumptions,
    ):
        self.case_id = case_id
        self.title = title
        self.constants = constants
        self.geometry = geometry
        self.material = material
        self.installation = installation
        self.benchmark = benchmark
        self.sources = sources
        self.assumptions = assumptions


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
        # TB 963 states the TB 880 Case 0 aluminium laminated sheath cross-section as 170 mm2.
        a_screen_elec_m2=170e-6,
    ),
    material=MaterialElectricalData(
        # IEC 60228 / TB 880-aligned benchmark primary input for 630 mm2 copper conductor.
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
