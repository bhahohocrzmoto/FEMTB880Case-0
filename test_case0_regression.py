# -*- coding: utf-8 -*-
"""Traceability self-check against the published TB 880 Case #0-1 benchmark.

I use this script to verify that the analytical IEC 60287 implementation still
reproduces the CIGRE TB 880 Case #0-1 benchmark values. The checks preserve a
clear audit trail between the Python code, the centralized case dataset, and the
published reference values that TB 963 later reuses for FEM benchmarking.

Run:
    python test_case0_regression.py
"""

import math

from analytical_solver_case0 import solve_case0
from cable_losses_tb880_case0 import create_case0_cables
from tb880_case0_data import CASE


def assert_close(name, value, ref, rel_tol=0.01, abs_tol=1e-9):
    if not math.isclose(value, ref, rel_tol=rel_tol, abs_tol=abs_tol):
        raise AssertionError(
            "{0} out of tolerance: value={1}, ref={2}, rel_tol={3}".format(name, value, ref, rel_tol)
        )


def main():
    default_circulating = CASE.assumptions.include_sheath_circulating_losses
    default_eddy = CASE.assumptions.include_sheath_eddy_losses
    default_F = CASE.assumptions.include_F_factor_for_eddy_reduction
    default_skin = CASE.assumptions.include_skin_effect
    default_prox = CASE.assumptions.include_proximity_effect
    default_diel = CASE.assumptions.include_dielectric_losses

    try:
        default_result = solve_case0(verbose=False)

        required = [
            "Tc_final_c", "Ts_final_c", "Wc_w_per_m", "Ws_w_per_m", "Wd_w_per_m", "I_max_a",
            "Rac_ohm_per_m", "Rs_ohm_per_m", "lambda1_prime", "lambda1_doubleprime",
        ]
        for k in required:
            if k not in default_result:
                raise AssertionError("Missing output key: {0}".format(k))
            if not math.isfinite(default_result[k]):
                raise AssertionError("Non-finite output for key: {0}".format(k))

        # I check the final AC conductor resistance against TB 880 p. 74 because this
        # validates the IEC 60287-1-1 Eq. (2) implementation at 90 degC.
        assert_close("Rac_90C", default_result["Rac_ohm_per_m"], CASE.benchmark.r_ac_90c_ohm_per_m, rel_tol=0.005)
        # I check the screen resistance against TB 880 p. 74 because the first-iteration
        # sheath loss factor depends directly on this benchmark quantity.
        assert_close("Rs_screen", default_result["Rs_ohm_per_m"], CASE.benchmark.r_screen_80c_ohm_per_m, rel_tol=0.03)
        # I check the final conductor temperature against TB 880 p. 76, where the rated
        # solution is defined by theta_c = 90 degC.
        assert_close("Tc_final_c", default_result["Tc_final_c"], CASE.benchmark.theta_core_final_c, rel_tol=0.001)
        # I check the final screen temperature against TB 880 pp. 74-76 because it is a
        # key convergence target for the benchmark iteration procedure.
        assert_close("Ts_final_c", default_result["Ts_final_c"], CASE.benchmark.theta_screen_final_c, rel_tol=0.01)
        # I check the final conductor loss against TB 880 pp. 74-76 to preserve the
        # published Joule heating reference value.
        assert_close("Wc_w_per_m", default_result["Wc_w_per_m"], CASE.benchmark.wc_final_w_per_m, rel_tol=0.01)
        # I check the final sheath loss against TB 880 pp. 74-76 because the
        # circulating-current result is the most sensitive benchmark output.
        assert_close("Ws_w_per_m", default_result["Ws_w_per_m"], CASE.benchmark.ws_final_w_per_m, rel_tol=0.03)
        # I check the dielectric loss against TB 880 pp. 71-72, where W_d is published
        # independently of the temperature iteration. After fixing the dielectric
        # geometry to use d_c at the conductor screen and D_i at the insulation OD,
        # the benchmark mismatch drops to about 0.14%, so 0.5% keeps this traceability
        # check tight enough to catch geometry regressions without being brittle.
        assert_close("Wd_w_per_m", default_result["Wd_w_per_m"], CASE.benchmark.wd_w_per_m, rel_tol=0.005)
        # I check the final current rating against TB 880 p. 76, where the published
        # Case #0-1 ampacity is 821.78 A.
        assert_close("I_max_a", default_result["I_max_a"], CASE.benchmark.i_final_a, rel_tol=0.01)

        # I enforce the design intent that the default benchmark path must use the
        # centralized sheath flags: circulating included and eddy excluded.
        if not CASE.assumptions.include_sheath_circulating_losses:
            raise AssertionError("Default benchmark path must include circulating sheath losses")
        if CASE.assumptions.include_sheath_eddy_losses:
            raise AssertionError("Default benchmark path must exclude eddy sheath losses")
        if default_result["lambda1_prime"] <= 0.0:
            raise AssertionError("Default benchmark path expected lambda1_prime > 0")
        if default_result["lambda1"] != default_result["lambda1_prime"]:
            raise AssertionError("Default benchmark path expected lambda1 = lambda1_prime")
        if default_result["lambda1_doubleprime"] <= 0.0:
            raise AssertionError("Default benchmark path should still compute lambda1_doubleprime")

        # I verify the IEC 60287-1-1:2023 Sec 5.3.7.1 eddy intermediates against the
        # TB 880 p. 79 benchmark values at the published 90/80 degC operating point.
        # I test the eddy-only scenario by setting the circulating flag to
        # False and the eddy flag to True. With the circulating flag off,
        # the F factor is automatically not applied regardless of the
        # dedicated F flag, because F only acts when both loss mechanisms
        # coexist. This configuration represents a sheath loss assembly
        # where only eddy currents contribute to lambda1.
        CASE.assumptions.include_sheath_circulating_losses = False
        CASE.assumptions.include_sheath_eddy_losses = True
        cable_eddy = create_case0_cables(I_rms_A=CASE.benchmark.i_final_a)["C02"]
        eddy_terms = cable_eddy.calculate_losses(
            CASE.benchmark.theta_core_final_c,
            80.0,
        )
        assert_close("beta1", eddy_terms["beta1"], 105.8022420, rel_tol=0.002)
        assert_close("gs", eddy_terms["gs"], 1.002449783, rel_tol=0.002)
        assert_close("lambda0", eddy_terms["lambda0"], 0.01354, rel_tol=0.01)
        assert_close("delta1", eddy_terms["delta1"], 0.08056, rel_tol=0.01)
        assert_close("delta2", eddy_terms["delta2"], 0.0, rel_tol=0.0, abs_tol=0.0)
        assert_close("lambda1_doubleprime_single", eddy_terms["lambda1_doubleprime"], 0.07696, rel_tol=0.02)
        if eddy_terms["thickness_term"] <= 0.0:
            raise AssertionError("IEC thickness term should now be computed and positive")
        # I verify that F = 1.0, meaning the correction is not applied, because
        # there are no circulating currents in this eddy-only scenario.
        if eddy_terms["F"] != 1.0:
            raise AssertionError(
                "Eddy-only mode must not apply F factor, but F={0}".format(
                    eddy_terms["F"]
                )
            )

        # I also check the optional comparison mode, where the central flags include
        # both circulating and eddy losses and therefore reduce the current rating.
        CASE.assumptions.include_sheath_circulating_losses = True
        CASE.assumptions.include_sheath_eddy_losses = True
        eddy_result = solve_case0(verbose=False)
        if eddy_result["lambda1_doubleprime"] <= 0.0:
            raise AssertionError("Eddy comparison mode expected lambda1_doubleprime > 0")
        if not math.isclose(
            eddy_result["lambda1"],
            eddy_result["lambda1_prime"] + eddy_result["lambda1_doubleprime"],
            rel_tol=1e-12,
            abs_tol=1e-12,
        ):
            raise AssertionError("Combined sheath mode expected lambda1 = lambda1_prime + lambda1_doubleprime")
        if not (eddy_result["I_max_a"] < default_result["I_max_a"]):
            raise AssertionError("Combined sheath mode should produce lower current rating")
        # I verify that F < 1.0 when both circulating and eddy losses are active,
        # confirming that the IEC 60287-1-1 Section 2.3.5 correction is applied
        # to reduce eddy losses in the presence of circulating currents.
        cable_combined_check = create_case0_cables(
            I_rms_A=CASE.benchmark.i_final_a
        )["C02"]
        combined_check_losses = cable_combined_check.calculate_losses(
            CASE.benchmark.theta_core_final_c,
            80.0,
        )
        if combined_check_losses["F"] >= 1.0:
            raise AssertionError(
                "Combined circ+eddy mode must apply F < 1.0, but F={0}".format(
                    combined_check_losses["F"]
                )
            )

        # I verify that setting include_F_factor_for_eddy_reduction = False
        # causes F to not be applied even when both circulating and eddy
        # losses are active. This lets the user study the impact of omitting
        # the F correction from the simplified non-Milliken approach.
        CASE.assumptions.include_F_factor_for_eddy_reduction = False
        cable_no_F = create_case0_cables(
            I_rms_A=CASE.benchmark.i_final_a
        )["C02"]
        losses_no_F = cable_no_F.calculate_losses(
            CASE.benchmark.theta_core_final_c,
            80.0,
        )
        if losses_no_F["F"] != 1.0:
            raise AssertionError(
                "F flag False must produce F=1.0, but F={0}".format(
                    losses_no_F["F"]
                )
            )
        # I also check that lambda1'' without F is larger than with F,
        # confirming the F factor reduces eddy losses as expected.
        if not (
            losses_no_F["lambda1_doubleprime"]
            > combined_check_losses["lambda1_doubleprime"]
        ):
            raise AssertionError(
                "lambda1'' without F must be larger than with F"
            )
        CASE.assumptions.include_F_factor_for_eddy_reduction = default_F

        # I verify the direct loss API too, so the reusable loss library stays
        # aligned with the central flag model when both sheath terms are disabled.
        CASE.assumptions.include_sheath_circulating_losses = False
        CASE.assumptions.include_sheath_eddy_losses = False
        cable_no_sheath = create_case0_cables(I_rms_A=CASE.benchmark.i_final_a)["C02"]
        losses_no_sheath = cable_no_sheath.calculate_losses(
            CASE.benchmark.theta_core_final_c,
            CASE.benchmark.theta_screen_final_c,
        )
        if losses_no_sheath["lambda1"] != 0.0:
            raise AssertionError("Both sheath flags false must produce lambda1 = 0")
        if losses_no_sheath["screen"] != 0.0:
            raise AssertionError("Both sheath flags false must produce Ws = 0")

        # I verify that area_mode correctly routes area-dependent quantities
        # through the full loss chain. I set A_screen_FE_model to a deliberately
        # different value from A_screen_elec so that any function using the
        # screen area will produce a different result depending on area_mode.
        # This confirms that Rs (from Rs_at_temp) and lambda1_doubleprime (from
        # the full eddy-current chain) both respect the area_mode flag. I also
        # directly test _beta1 with a fixed Rs value so I can confirm that the
        # helper itself honors area_mode when recovering sheath resistivity.
        cable_fem_test = create_case0_cables(I_rms_A=CASE.benchmark.i_final_a)["C02"]
        cable_fem_test.A_screen_FE_model = 200e-6
        losses_analytical = cable_fem_test.calculate_losses(90.0, 80.0, area_mode="analytical")
        losses_fem = cable_fem_test.calculate_losses(90.0, 80.0, area_mode="fem")
        if losses_analytical["Rs"] == losses_fem["Rs"]:
            raise AssertionError(
                "area_mode='fem' must produce different Rs when "
                "A_screen_FE_model != A_screen_elec"
            )
        if losses_analytical["lambda1_doubleprime"] == losses_fem["lambda1_doubleprime"]:
            raise AssertionError(
                "area_mode='fem' must produce different lambda1_doubleprime when "
                "A_screen_FE_model != A_screen_elec"
            )
        beta1_from_fixed_rs_analytical = cable_fem_test._beta1(
            losses_analytical["Rs"], area_mode="analytical"
        )
        beta1_from_fixed_rs_fem = cable_fem_test._beta1(
            losses_analytical["Rs"], area_mode="fem"
        )
        if beta1_from_fixed_rs_analytical == beta1_from_fixed_rs_fem:
            raise AssertionError(
                "_beta1 must honor area_mode when recovering sheath resistivity"
            )

        # -- Skin-effect flag test --
        # I verify that disabling the skin-effect flag produces a lower
        # Rac than the default, because ys >= 0 always increases Rac.
        CASE.assumptions.include_skin_effect = False
        cable_no_skin = create_case0_cables(I_rms_A=CASE.benchmark.i_final_a)["C02"]
        losses_no_skin = cable_no_skin.calculate_losses(90.0, 80.0)
        CASE.assumptions.include_skin_effect = True
        cable_with_skin = create_case0_cables(I_rms_A=CASE.benchmark.i_final_a)["C02"]
        losses_with_skin = cable_with_skin.calculate_losses(90.0, 80.0)
        if not (losses_no_skin["Rac"] < losses_with_skin["Rac"]):
            raise AssertionError(
                "Disabling skin effect must reduce Rac"
            )

        # -- Proximity-effect flag test --
        # I verify that disabling the proximity-effect flag produces a lower
        # Rac than the default, because yp >= 0 always increases Rac.
        CASE.assumptions.include_proximity_effect = False
        cable_no_prox = create_case0_cables(I_rms_A=CASE.benchmark.i_final_a)["C02"]
        losses_no_prox = cable_no_prox.calculate_losses(90.0, 80.0)
        CASE.assumptions.include_proximity_effect = True
        if not (losses_no_prox["Rac"] < losses_with_skin["Rac"]):
            raise AssertionError(
                "Disabling proximity effect must reduce Rac"
            )

        # -- Dielectric-loss flag test --
        # I verify that disabling the dielectric-loss flag produces Wd = 0.
        CASE.assumptions.include_dielectric_losses = False
        cable_no_diel = create_case0_cables(I_rms_A=CASE.benchmark.i_final_a)["C02"]
        losses_no_diel = cable_no_diel.calculate_losses(90.0, 80.0)
        CASE.assumptions.include_dielectric_losses = True
        if losses_no_diel["dielectric"] != 0.0:
            raise AssertionError(
                "Disabling dielectric losses must produce Wd = 0"
            )

        print("TB880 Case #0 regression checks passed.")
    finally:
        CASE.assumptions.include_sheath_circulating_losses = default_circulating
        CASE.assumptions.include_sheath_eddy_losses = default_eddy
        CASE.assumptions.include_F_factor_for_eddy_reduction = default_F
        CASE.assumptions.include_skin_effect = default_skin
        CASE.assumptions.include_proximity_effect = default_prox
        CASE.assumptions.include_dielectric_losses = default_diel



if __name__ == "__main__":
    main()
