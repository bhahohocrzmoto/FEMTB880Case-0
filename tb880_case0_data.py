# -*- coding: utf-8 -*-
"""Central TB 880 Case #0 data model for the IEC 60287 benchmark workflow.

I use this module as the single source of truth for the CIGRE TB 880 Case #0
inputs and published reference results that drive both the analytical IEC 60287
solver and the FEM coupling workflow. I store the cable geometry from TB 880
Figure 25 and Figure 26, and I store the material constants that trace back to
IEC 60287-1-1:2023 Table 1 and Table 2 so that every downstream calculation can
be audited against the benchmark literature.
"""
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
        # I keep vacuum permittivity and permeability here so that IEC 60287-1-1
        # Eq. (15) and the sheath reactance expression can be evaluated without
        # duplicating fundamental constants across solver modules.
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
        # I store the conductor diameter from TB 880 Figure 25 because it defines
        # the physical copper region used for both the IEC proximity factor and the
        # geometric FE heat source area.
        self.d_cond_m = d_cond_m
        # I store the diameter over the conductor screen from TB 880 Figure 25
        # because IEC 60287-1-1 Eq. (15) uses it as d_c in the cable capacitance.
        self.d_inner_semicon_m = d_inner_semicon_m
        # I store the diameter over the main insulation from TB 880 Figure 25
        # because this is the D_i used for dielectric capacitance in the TB 880
        # Case #0 benchmark, while also keeping the insulation FE annulus aligned.
        self.d_ins_m = d_ins_m
        # I store the diameter over the insulation screen from TB 880 Figure 25
        # because it is relevant for screen/sheath geometry and thickness recovery,
        # not the D_i used in the dielectric capacitance formula for this benchmark.
        self.d_outer_semicon_m = d_outer_semicon_m
        # I store the mean sheath diameter from TB 880 Figure 25 because IEC
        # 60287-1-1 Section 2.3.1 uses the mean sheath diameter in reactance X.
        self.d_screen_mean_m = d_screen_mean_m
        # I store the outer diameter of the metallic screen from TB 880 Figure 25
        # so that the sheath thickness can be recovered for FE area estimation.
        self.d_screen_out_m = d_screen_out_m
        # I store the overall cable diameter from TB 880 Figure 25 because touching
        # trefoil means the cable centre spacing s equals this outer diameter.
        self.d_oversheath_m = d_oversheath_m

        # I keep the nominal IEC electrical conductor area as 630 mm2 copper,
        # which is distinct from the geometric circular FE area used for heat load
        # conversion in the coupled FEM model.
        self.a_cond_elec_m2 = a_cond_elec_m2
        # I keep the metallic sheath electrical area because TB 963 p. 104 states
        # the reused TB 880 sheath corresponds to a 170 mm2 aluminium laminate.
        self.a_screen_elec_m2 = a_screen_elec_m2

    @property
    def screen_thickness_m(self):
        # I recover the sheath thickness from the TB 880 layer diameters because
        # the thin-ring FE area estimate depends on the screen radial thickness.
        return 0.5 * (self.d_screen_out_m - self.d_outer_semicon_m)

    @property
    def a_cond_geom_m2(self):
        # I use the true geometric conductor circle here because TB 963 requires
        # volumetric heat generation in W/m^3 for the copper FE body.
        return math.pi * (self.d_cond_m / 2.0) ** 2

    @property
    def a_innerins_geom_m2(self):
        # I define this annular area from conductor OD to insulation OD because the
        # repository groups the inner semicon and main insulation into one FE body,
        # and TB 963 converts dielectric W/m into volumetric W/m^3 using that body area.
        return math.pi * ((self.d_ins_m / 2.0) ** 2 - (self.d_cond_m / 2.0) ** 2)

    @property
    def a_screen_geom_m2(self):
        # I use the thin-ring approximation because the sheath is slender relative
        # to its diameter and this gives the FE body area needed for W/m^3 loading.
        return math.pi * self.d_screen_mean_m * self.screen_thickness_m

    @property
    def spacing_touching_trefoil_m(self):
        # I set the touching trefoil spacing to the cable outer diameter because
        # TB 880 assumes the three single-core cables are in touching trefoil.
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
        # I store the IEC 60228 maximum DC resistance at 20 degC for the 630 mm2
        # Class 2 stranded copper conductor because TB 880 uses this benchmark input
        # directly instead of deriving R_dc(20 degC) from resistivity and area.
        self.r_cond_dc_20_ohm_per_m = r_cond_dc_20_ohm_per_m

        # I keep the copper resistivity from IEC 60287-1-1:2023 Table 1 so that
        # the dataset documents the canonical material property even though the
        # benchmark uses the published IEC 60228 resistance as its main input.
        self.rho_cond_20_ohm_m = rho_cond_20_ohm_m
        # I keep the aluminium sheath resistivity from IEC 60287-1-1:2023 Table 1
        # because the sheath resistance is evaluated from resistivity and area.
        self.rho_screen_20_ohm_m = rho_screen_20_ohm_m

        # I use the conductor temperature coefficient from IEC 60287-1-1:2023
        # Table 1 in the Eq. (1) resistance-temperature correction.
        self.alpha_cond_20_per_k = alpha_cond_20_per_k
        # I use the screen temperature coefficient from IEC 60287-1-1:2023 Table 1
        # so that sheath resistance follows the same IEC temperature convention.
        self.alpha_screen_20_per_k = alpha_screen_20_per_k

        # I store k_s from IEC 60287-1-1:2023 Table 2 for a round stranded
        # conductor with extruded insulation because Eq. (5) depends on it.
        self.ks = ks
        # I store k_p from IEC 60287-1-1:2023 Table 2 for the same conductor type
        # because Eq. (9) and Eq. (10) use it in the proximity correction.
        self.kp = kp

        # I treat epsilon_r as the standard IEC XLPE relative permittivity used in
        # IEC 60287-1-1 Eq. (15) for the dielectric loss capacitance.
        self.eps_r = eps_r
        # I treat tan(delta) as the standard IEC XLPE dielectric loss factor used
        # in IEC 60287-1-1 Eq. (14) for benchmark-faithful dielectric heating.
        self.tan_delta = tan_delta

        # I store semicon thermal resistivity so that the published TB 880 thermal
        # resistances can be traced back to the layer-by-layer IEC 60287-2-1 model.
        self.rho_semicon_k_m_per_w = rho_semicon_k_m_per_w
        # I store insulation thermal resistivity for the same reason, namely to
        # preserve the T1 benchmark lineage from TB 880 pp. 70-71.
        self.rho_ins_k_m_per_w = rho_ins_k_m_per_w
        # I store oversheath thermal resistivity because IEC 60287-2-1 Section 2.3
        # and TB 880 pp. 70-71 use it in T3 with the trefoil correction.
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
        # I store the bonding mode explicitly because solid means both sheath ends
        # are bonded, which activates circulating current losses lambda1' per IEC
        # 60287-1-1:2023 Section 2.3.1 for single-core trefoil cables.
        self.bonding = bonding

    @property
    def u0_v(self):
        # I convert line-to-line voltage to phase-to-earth voltage because IEC
        # 60287-1-1 Eq. (14) uses U_0 in the dielectric loss expression.
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
        # I store T1 from TB 880 pp. 70-71 because it is the conductor-to-sheath
        # insulation thermal resistance used in the IEC 60287-2-1 network.
        self.t1_k_m_per_w = t1_k_m_per_w
        # I store T2 for completeness even though it is zero in this benchmark,
        # which reflects the absence of a separate bedding layer in the model.
        self.t2_k_m_per_w = t2_k_m_per_w
        # I store T3 from TB 880 pp. 70-71 because it includes the IEC 60287-2-1
        # Section 2.3 oversheath resistance with the 1.6 trefoil correction.
        self.t3_k_m_per_w = t3_k_m_per_w
        # I store T4 from TB 880 pp. 71-72 because it represents the external soil
        # thermal resistance used in the final buried-trefoil benchmark equation.
        self.t4_k_m_per_w = t4_k_m_per_w

        # I store W_d from TB 880 pp. 71-72 as the published dielectric loss value.
        self.wd_w_per_m = wd_w_per_m
        # I store the final ampacity from TB 880 p. 76 so that analytical and FEM
        # workflows can be checked against the published benchmark target.
        self.i_final_a = i_final_a
        # I store the final conductor loss from TB 880 pp. 74-76 for regression.
        self.wc_final_w_per_m = wc_final_w_per_m
        # I store the final sheath loss from TB 880 pp. 74-76 for regression.
        self.ws_final_w_per_m = ws_final_w_per_m

        # I store the final conductor temperature from TB 880 p. 76 because the
        # benchmark is defined at theta_c = 90 degC for the rated current.
        self.theta_core_final_c = theta_core_final_c
        # I store the final screen temperature from TB 880 pp. 74-76 so that the
        # coupled iterations can be validated against the published solution.
        self.theta_screen_final_c = theta_screen_final_c
        # I store the oversheath temperature from TB 880 pp. 74-76 for completeness
        # even though the present analytical model solves only core and sheath nodes.
        self.theta_oversheath_final_c = theta_oversheath_final_c

        # According to TB 880 pp. 66-68, I use the initial screen temperature guess
        # of 80 degC to start the iterative benchmark workflow.
        self.theta_screen_init_c = theta_screen_init_c
        # I store the final AC conductor resistance from TB 880 p. 74 because it is
        # a sensitive checkpoint for Eq. (2) plus the skin and proximity factors.
        self.r_ac_90c_ohm_per_m = r_ac_90c_ohm_per_m
        # I store the screen resistance at about 80 degC from TB 880 p. 74 because
        # it is the key first-iteration input for the sheath loss calculation.
        self.r_screen_80c_ohm_per_m = r_screen_80c_ohm_per_m
        # I store the first-iteration circulating loss factor from TB 880 p. 74 so
        # that the Section 2.3.1 implementation can be regression checked directly.
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
        # I record that the InnerIns FE body is treated as one annulus from the
        # conductor OD to the insulation OD because that area is used to convert
        # IEC dielectric W/m into volumetric W/m^3 for the FEM body in TB 963.
        self.innerins_area_convention = innerins_area_convention

        # I flag that the screen FE area is intentionally treated as equivalent to
        # the electrical area in this benchmark so that the FEM load conversion and
        # the IEC resistance model remain aligned with the simplified dataset.
        self.screen_fe_area_equals_electrical_area_assumed = screen_fe_area_equals_electrical_area_assumed

        # I capture the touching trefoil assumption because IEC 60287 Eq. (10) and
        # the sheath reactance expression both depend on the spacing convention.
        self.touching_trefoil_spacing_assumed = touching_trefoil_spacing_assumed

        # I state which sheath loss interpretation is considered benchmark-faithful
        # so that alternate comparison modes can be distinguished from the baseline.
        self.sheath_loss_interpretation = sheath_loss_interpretation

        # I document the design decision that solid-bonded eddy losses are excluded
        # by default because the published TB 880 Case #0-1 benchmark results match
        # the circulating-current-only interpretation rather than a combined value.
        self.solid_bonded_include_eddy_default = solid_bonded_include_eddy_default


# ----------------------------------------------------------------------
# Runtime FEM area data read from the active ANSYS model
# ----------------------------------------------------------------------
class RuntimeFEAreas(object):
    def __init__(self):
        # I keep the actual FE body areas read from Mechanical separate from the
        # diameter-derived reference areas so the FEM source conversion can be
        # audited without changing the benchmark geometry inputs.
        self.a_model_read_m2_by_region = {}


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
        self.runtime_fe_areas = RuntimeFEAreas()


TB880_CASE_0 = TB880Case0(
    case_id="TB880_CASE_0",
    title="CIGRE TB 880 Introductory Case Study #0 - 132 kV single-core XLPE cable, touching trefoil, directly buried",
    constants=PhysicalConstants(),
    geometry=Geometry(
        # According to TB 880 Figure 25 and Figure 26, I use these diameters as the
        # authoritative Case #0 layer geometry for the entire benchmark workflow.
        d_cond_m=30.3e-3,
        d_inner_semicon_m=33.3e-3,
        d_ins_m=64.3e-3,
        d_outer_semicon_m=66.9e-3,
        d_screen_mean_m=67.7e-3,
        d_screen_out_m=68.5e-3,
        d_oversheath_m=75.5e-3,

        # I use the IEC nominal electrical cross-sections from the published case
        # because the benchmark electrical model is defined by nominal areas.
        a_cond_elec_m2=630e-6,
        # TB 963 p. 104 states the reused TB 880 aluminium laminated sheath has an
        # electrical cross-section of 170 mm2, so I store that directly here.
        a_screen_elec_m2=170e-6,
    ),
    material=MaterialElectricalData(
        # I keep the IEC 60228 maximum conductor resistance at 20 degC as a primary
        # TB 880 input rather than recomputing it from rho and nominal area.
        r_cond_dc_20_ohm_per_m=2.83e-5,

        # IEC 60287-1-1:2023 Table 1 material resistivities at 20 degC.
        rho_cond_20_ohm_m=1.7241e-8,
        rho_screen_20_ohm_m=2.84e-8,

        # IEC 60287-1-1:2023 Table 1 temperature coefficients at 20 degC.
        alpha_cond_20_per_k=3.93e-3,
        alpha_screen_20_per_k=4.03e-3,

        # IEC 60287-1-1:2023 Table 2 gives k_s = 1.0 and k_p = 1.0 for a round
        # stranded conductor with extruded insulation, which matches Case #0.
        ks=1.0,
        kp=1.0,

        # I use the standard IEC XLPE dielectric constants for Eq. (14) and Eq. (15).
        eps_r=2.5,
        tan_delta=1.0e-3,

        # I retain the thermal resistivities used to produce the published TB 880
        # thermal resistances and temperature ladder for Case #0-1.
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
        # I use solid bonding here, meaning both sheath ends are bonded and the IEC
        # 60287-1-1 Section 2.3.1 circulating-current term lambda1' is active.
        bonding="solid",
    ),
    benchmark=Benchmark(
        # TB 880 pp. 70-71 publish the temperature-independent thermal resistances.
        t1_k_m_per_w=0.4198714890,
        t2_k_m_per_w=0.0,
        t3_k_m_per_w=0.08671937479,
        t4_k_m_per_w=1.594692892,

        # TB 880 pp. 71-72 and pp. 74-76 publish the benchmark losses and rating.
        wd_w_per_m=0.3851382172,
        i_final_a=821.7763334392,
        wc_final_w_per_m=26.6895,
        ws_final_w_per_m=7.8442,

        # TB 880 pp. 74-76 publish the final conductor, screen, and oversheath temperatures.
        theta_core_final_c=90.0,
        theta_screen_final_c=78.7130,
        theta_oversheath_final_c=75.68,

        # TB 880 pp. 66-68 starts the iterative process with an 80 degC screen guess.
        theta_screen_init_c=80.0,
        # TB 880 p. 74 publishes these intermediate electrical checkpoints.
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
