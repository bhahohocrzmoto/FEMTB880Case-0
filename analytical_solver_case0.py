# -*- coding: utf-8 -*-
# ============================================================
# analytical_solver_case0.py
# ------------------------------------------------------------
# Standalone analytical iterative solver for TB 880 Case #0.
# Uses the centralized loss API (Cable.calculate_losses) as
# the single source of truth for electrical losses.
# ============================================================

import math
import cable_losses_tb880_case0 as loss_mod
from tb880_case0_data import CASE


def solve_case0(verbose=True, sheath_eddy_policy="auto"):
    T1 = CASE.benchmark.t1_k_m_per_w
    T2 = CASE.benchmark.t2_k_m_per_w
    T3 = CASE.benchmark.t3_k_m_per_w
    T4 = CASE.benchmark.t4_k_m_per_w

    T_amb = CASE.installation.ambient_temp_c
    I_load = CASE.benchmark.i_final_a

    cables = loss_mod.create_case0_cables(
        I_rms_A=I_load,
        bonding=CASE.installation.bonding,
        sheath_eddy_policy=sheath_eddy_policy,
    )
    cable = cables["C02"]

    if verbose:
        print("--- IEC 60287 Iterative Temperature Solver ---")
        print("Load Current: {0} A | Ambient Temp: {1} C".format(I_load, T_amb))

    Tc = T_amb
    Ts = CASE.benchmark.theta_screen_init_c

    max_iter = 100
    tolerance = 1e-4

    converged_iterations = None
    losses = None
    for iteration in range(max_iter):
        losses = cable.calculate_losses(Tc, Ts)
        Wc = losses["core"]
        Ws = losses["screen"]
        Wd = losses["dielectric"]

        delta_T_ins = (Wc + 0.5 * Wd) * T1
        delta_T_sheath = (Wc + Ws + Wd) * (T2 + T3)
        delta_T_soil = (Wc + Ws + Wd) * T4

        Tc_new = T_amb + delta_T_ins + delta_T_sheath + delta_T_soil
        Ts_new = Tc_new - delta_T_ins

        error_Tc = abs(Tc_new - Tc)
        Tc, Ts = Tc_new, Ts_new

        if error_Tc < tolerance:
            converged_iterations = iteration + 1
            if verbose:
                print("Converged in {0} iterations.".format(converged_iterations))
                print("Final Conductor Temp (Tc) = {0:.3f} C".format(Tc))
                print("Final Sheath Temp (Ts)    = {0:.3f} C".format(Ts))
                print(
                    "Losses: Wc = {0:.2f} W/m, Wd = {1:.2f} W/m, Ws = {2:.2f} W/m\n".format(
                        losses["core"], losses["dielectric"], losses["screen"]
                    )
                )
            break

    if verbose:
        print("--- IEC 60287 Analytical Ampacity (For Tc = 90 C) ---")

    Tc_max = CASE.benchmark.theta_core_final_c
    Ts_guess = CASE.benchmark.theta_screen_init_c
    lambda2 = 0.0

    ampacity_losses = None
    for _ in range(5):
        ampacity_losses = cable.calculate_losses(Tc_max, Ts_guess)
        Rac_90 = ampacity_losses["Rac"]
        lambda1_90 = ampacity_losses["lambda1"]
        Wd = ampacity_losses["dielectric"]

        delta_theta = Tc_max - T_amb
        numerator = delta_theta - Wd * (0.5 * T1 + T2 + T3 + T4)

        denominator = (
            Rac_90 * T1
            + Rac_90 * (1 + lambda1_90) * T2
            + Rac_90 * (1 + lambda1_90 + lambda2) * (T3 + T4)
        )

        I_max = math.sqrt(numerator / denominator)

        cable.I = I_max
        ampacity_losses = cable.calculate_losses(Tc_max, Ts_guess)
        Wc_max = ampacity_losses["core"]
        Ts_guess = Tc_max - (Wc_max + 0.5 * Wd) * T1

    if verbose:
        print("Maximum Current Rating (I) = {0:.2f} A".format(I_max))
        print("Sheath Temp at Max Load    = {0:.3f} C".format(Ts_guess))

    return {
        "Tc_final_c": Tc,
        "Ts_final_c": Ts,
        "Wc_w_per_m": losses["core"],
        "Ws_w_per_m": losses["screen"],
        "Wd_w_per_m": losses["dielectric"],
        "lambda1": losses["lambda1"],
        "lambda1_prime": losses["lambda1_prime"],
        "lambda1_doubleprime": losses["lambda1_doubleprime"],
        "Rac_ohm_per_m": losses["Rac"],
        "Rdc_ohm_per_m": losses["Rdc"],
        "Rs_ohm_per_m": losses["Rs"],
        "I_max_a": I_max,
        "Ts_at_Imax_c": Ts_guess,
        "iterations": converged_iterations,
        "sheath_eddy_policy": sheath_eddy_policy,
    }


def run_analytical_solver():
    return solve_case0(verbose=True)


if __name__ == "__main__":
    run_analytical_solver()
