"""Lightweight regression/self-check for TB880 Case 0 analytical workflow.

Run:
    python test_case0_regression.py
"""

import math

from analytical_solver_case0 import solve_case0
from tb880_case0_data import CASE


def assert_close(name, value, ref, rel_tol=0.03, abs_tol=1e-6):
    if not math.isclose(value, ref, rel_tol=rel_tol, abs_tol=abs_tol):
        raise AssertionError(
            "{0} out of tolerance: value={1}, ref={2}, rel_tol={3}".format(name, value, ref, rel_tol)
        )


def main():
    result = solve_case0(verbose=False)

    required = ["Tc_final_c", "Ts_final_c", "Wc_w_per_m", "Ws_w_per_m", "Wd_w_per_m", "I_max_a"]
    for k in required:
        if k not in result:
            raise AssertionError("Missing output key: {0}".format(k))
        if not math.isfinite(result[k]):
            raise AssertionError("Non-finite output for key: {0}".format(k))

    if result["Wc_w_per_m"] <= 0 or result["Ws_w_per_m"] <= 0 or result["Wd_w_per_m"] <= 0:
        raise AssertionError("Expected positive core/screen/dielectric losses")

    # Benchmark alignment checks (tolerances keep this robust for minor numerical drift)
    assert_close("Tc_final_c", result["Tc_final_c"], CASE.benchmark.theta_core_final_c, rel_tol=0.03)
    assert_close("Ts_final_c", result["Ts_final_c"], CASE.benchmark.theta_screen_final_c, rel_tol=0.04)
    assert_close("Wc_w_per_m", result["Wc_w_per_m"], CASE.benchmark.wc_final_w_per_m, rel_tol=0.05)
    assert_close("Ws_w_per_m", result["Ws_w_per_m"], CASE.benchmark.ws_final_w_per_m, rel_tol=0.25)
    assert_close("Wd_w_per_m", result["Wd_w_per_m"], CASE.benchmark.wd_w_per_m, rel_tol=0.08)
    assert_close("I_max_a", result["I_max_a"], CASE.benchmark.i_final_a, rel_tol=0.03)

    print("TB880 Case 0 regression self-check passed.")


if __name__ == "__main__":
    main()
