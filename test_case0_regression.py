"""Regression/self-check for TB880 Case #0 benchmark-faithful workflow.

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

    assert_close("Rac_90C", default_result["Rac_ohm_per_m"], CASE.benchmark.r_ac_90c_ohm_per_m, rel_tol=0.005)
    assert_close("Rs_screen", default_result["Rs_ohm_per_m"], CASE.benchmark.r_screen_80c_ohm_per_m, rel_tol=0.03)
    assert_close("Tc_final_c", default_result["Tc_final_c"], CASE.benchmark.theta_core_final_c, rel_tol=0.001)
    assert_close("Ts_final_c", default_result["Ts_final_c"], CASE.benchmark.theta_screen_final_c, rel_tol=0.01)
    assert_close("Wc_w_per_m", default_result["Wc_w_per_m"], CASE.benchmark.wc_final_w_per_m, rel_tol=0.01)
    assert_close("Ws_w_per_m", default_result["Ws_w_per_m"], CASE.benchmark.ws_final_w_per_m, rel_tol=0.03)
    assert_close("Wd_w_per_m", default_result["Wd_w_per_m"], CASE.benchmark.wd_w_per_m, rel_tol=0.08)
    assert_close("I_max_a", default_result["I_max_a"], CASE.benchmark.i_final_a, rel_tol=0.01)

    # Default benchmark-faithful path for solid bonding must neglect eddy losses.
    if default_result["lambda1_doubleprime"] != 0.0:
        raise AssertionError("Default solid-bonded benchmark path must not include eddy losses")

    # Optional comparison mode: include eddy in solid bonding -> lower ampacity.
    eddy_result = solve_case0(verbose=False, sheath_eddy_policy="include_for_solid")
    if eddy_result["lambda1_doubleprime"] <= 0.0:
        raise AssertionError("Eddy comparison mode expected lambda1_doubleprime > 0")
    if not (eddy_result["I_max_a"] < default_result["I_max_a"]):
        raise AssertionError("Eddy comparison mode should produce lower current rating")

    # Also verify direct loss API default does not add eddy.
    cable_default = create_case0_cables(I_rms_A=CASE.benchmark.i_final_a)["C02"]
    losses_default = cable_default.calculate_losses(CASE.benchmark.theta_core_final_c, CASE.benchmark.theta_screen_final_c)
    if losses_default["lambda1_doubleprime"] != 0.0:
        raise AssertionError("Loss API default must keep lambda1_doubleprime = 0 for solid bonding")

    print("TB880 Case #0 regression checks passed.")


if __name__ == "__main__":
    main()
