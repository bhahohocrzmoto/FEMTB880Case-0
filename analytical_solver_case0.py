# -*- coding: utf-8 -*-
"""Iterative IEC 60287 analytical solver for TB 880 Case #0.

I implement the benchmark procedure described in CIGRE TB 880 pp. 66-68. This
module performs two related calculations: first, it iterates the IEC 60287-2-1
thermal network to recover the steady-state conductor and sheath temperatures at
the benchmark load current; second, it solves the IEC general rating equation to
recover the ampacity that corresponds to theta_c = 90 degC.
"""
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


def solve_case0(verbose=True):
    T1 = CASE.benchmark.t1_k_m_per_w
    T2 = CASE.benchmark.t2_k_m_per_w
    T3 = CASE.benchmark.t3_k_m_per_w
    T4 = CASE.benchmark.t4_k_m_per_w

    T_amb = CASE.installation.ambient_temp_c
    I_load = CASE.benchmark.i_final_a

    cables = loss_mod.create_case0_cables(
        I_rms_A=I_load,
        bonding=CASE.installation.bonding,
    )
    cable = cables["C02"]

    if verbose:
        print("--- IEC 60287 Iterative Temperature Solver ---")
        print("Load Current: {0} A | Ambient Temp: {1} C".format(I_load, T_amb))

    Tc = T_amb
    # According to TB 880 pp. 66-68, I start with the published 80 degC screen
    # temperature guess because the losses depend on sheath temperature.
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

        # I apply the IEC 60287-2-1 thermal ladder exactly as TB 880 pp. 66-68
        # describes it: half of W_d heats inward through T1, while the full sum of
        # W_c + W_s + W_d crosses the outer thermal resistances.
        delta_T_ins = (Wc + 0.5 * Wd) * T1
        delta_T_sheath = (Wc + Ws + Wd) * (T2 + T3)
        delta_T_soil = (Wc + Ws + Wd) * T4

        # I recover the new conductor temperature as ambient plus the complete
        # temperature rise decomposition from the IEC 60287-2-1 general equation.
        Tc_new = T_amb + delta_T_ins + delta_T_sheath + delta_T_soil
        # I update the screen temperature as theta_s = theta_c - delta_T_ins because
        # the screen sits outside the insulation drop in the thermal network.
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
    # I set lambda2 = 0 because this cable has no armour, so the IEC general rating
    # equation includes only conductor and sheath losses in the denominator.
    lambda2 = 0.0

    ampacity_max_iter = 100
    ampacity_tol_I = 1e-6
    ampacity_tol_Ts = 1e-6
    ampacity_losses = None
    I_prev = None
    error_I = float("inf")
    error_Ts = float("inf")
    converged_ampacity_iterations = None
    for iteration in range(ampacity_max_iter):
        Ts_prev = Ts_guess
        ampacity_losses = cable.calculate_losses(Tc_max, Ts_guess)
        Rac_90 = ampacity_losses["Rac"]
        lambda1_90 = ampacity_losses["lambda1"]
        Wd = ampacity_losses["dielectric"]

        delta_theta = Tc_max - T_amb
        # I use the IEC 60287-2-1 general rating equation numerator as the available
        # temperature rise after subtracting the dielectric heating burden.
        numerator = delta_theta - Wd * (0.5 * T1 + T2 + T3 + T4)

        # I use the denominator form from the IEC rating equation, where conductor
        # and sheath loss multipliers are weighted by T1, T2, T3, and T4.
        denominator = (
            Rac_90 * T1
            + Rac_90 * (1 + lambda1_90) * T2
            + Rac_90 * (1 + lambda1_90 + lambda2) * (T3 + T4)
        )

        I_max = math.sqrt(numerator / denominator)

        cable.I = I_max
        ampacity_losses = cable.calculate_losses(Tc_max, Ts_guess)
        Wc_max = ampacity_losses["core"]
        # I iterate the rated operating point until both the current rating and the
        # sheath temperature are converged in this fixed-point refinement.
        Ts_new = Tc_max - (Wc_max + 0.5 * Wd) * T1

        if I_prev is not None:
            error_I = abs(I_max - I_prev)
            error_Ts = abs(Ts_new - Ts_prev)
            if error_I < ampacity_tol_I and error_Ts < ampacity_tol_Ts:
                Ts_guess = Ts_new
                converged_ampacity_iterations = iteration + 1
                break

        I_prev = I_max
        Ts_guess = Ts_new
    else:
        raise RuntimeError(
            "Ampacity iteration did not converge within "
            f"{ampacity_max_iter} iterations. "
            f"Last error_I={error_I:.3e}, error_Ts={error_Ts:.3e}"
        )

    if verbose:
        print("Maximum Current Rating (I) = {0:.2f} A".format(I_max))
        print("Sheath Temp at Max Load    = {0:.3f} C".format(Ts_guess))
        print("Ampacity iterations used   = {0}".format(converged_ampacity_iterations))
        print("Final ampacity error_I     = {0:.3e} A".format(error_I))
        print("Final ampacity error_Ts    = {0:.3e} C".format(error_Ts))

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
        "ampacity_iterations": converged_ampacity_iterations,
    }


def run_analytical_solver():
    return solve_case0(verbose=True)


if __name__ == "__main__":
    run_analytical_solver()
