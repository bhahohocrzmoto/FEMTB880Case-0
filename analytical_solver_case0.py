# ============================================================
# analytical_solver_case0.py
# ------------------------------------------------------------
# Standalone analytical iterative solver for TB 880 Case 0.
# Uses the Cable class to calculate temperature-dependent losses
# and pushes them through the IEC thermal resistance network.
# ============================================================

import math
import cable_losses_tb880_case0 as loss_mod
from tb880_case0_data import CASE


def solve_case0(verbose=True):
    # --------------------------------------------------------
    # 1. INPUT PARAMETERS (centralized TB 880 Case 0 data)
    # --------------------------------------------------------
    T1 = CASE.benchmark.t1_k_m_per_w
    T2 = CASE.benchmark.t2_k_m_per_w
    T3 = CASE.benchmark.t3_k_m_per_w
    T4 = CASE.benchmark.t4_k_m_per_w

    T_amb = CASE.installation.ambient_temp_c
    I_load = CASE.benchmark.i_final_a

    # Initialize the Cable (using the central cable C02 as the hottest)
    cables = loss_mod.create_case0_cables(I_rms_A=I_load, bonding=CASE.installation.bonding)
    cable = cables["C02"]

    # --------------------------------------------------------
    # 2. ITERATIVE TEMPERATURE SOLVER
    # --------------------------------------------------------
    if verbose:
        print("--- IEC 60287 Iterative Temperature Solver ---")
        print("Load Current: {0} A | Ambient Temp: {1} C".format(I_load, T_amb))

    # Initial guesses
    Tc = T_amb
    Ts = CASE.benchmark.theta_screen_init_c

    max_iter = 100
    tolerance = 1e-4

    converged_iterations = None
    for iteration in range(max_iter):
        # a. Calculate temperature-dependent electrical resistances
        Rac = cable.Rac_at_temp(Tc)
        Rs = cable.sheath_resistance_at_temp(Ts)

        # b. Calculate Losses [W/m]
        Wc = (I_load ** 2) * Rac
        Wd = cable.dielectric_loss_W_per_m()

        # Sheath losses (Circulating + Eddy with F factor)
        lambda1_prime = cable._lambda1_prime(Rs, Rac, cable.sheath_reactance_X())
        lambda1_dprime = cable._lambda1_doubleprime(Rs, Rac, Ts)
        lambda1 = lambda1_prime + lambda1_dprime
        Ws = Wc * lambda1

        # Armour losses (lambda2) are 0 for this case
        lambda2 = 0.0

        # c. Calculate Temperature Drops across the thermal network
        delta_T_ins = (Wc + 0.5 * Wd) * T1
        delta_T_sheath = (Wc + Ws + Wd) * (T2 + T3)
        delta_T_soil = (Wc + Ws + Wd) * T4

        # d. Update Temperatures
        Tc_new = T_amb + delta_T_ins + delta_T_sheath + delta_T_soil
        Ts_new = Tc_new - delta_T_ins

        # e. Check Convergence
        error_Tc = abs(Tc_new - Tc)

        Tc = Tc_new
        Ts = Ts_new

        if error_Tc < tolerance:
            converged_iterations = iteration + 1
            if verbose:
                print("Converged in {0} iterations.".format(converged_iterations))
                print("Final Conductor Temp (Tc) = {0:.3f} C".format(Tc))
                print("Final Sheath Temp (Ts)    = {0:.3f} C".format(Ts))
                print("Losses: Wc = {0:.2f} W/m, Wd = {1:.2f} W/m, Ws = {2:.2f} W/m\n".format(Wc, Wd, Ws))
            break

    # --------------------------------------------------------
    # 3. ANALYTICAL AMPACITY (CURRENT RATING)
    # --------------------------------------------------------
    if verbose:
        print("--- IEC 60287 Analytical Ampacity (For Tc = 90 C) ---")

    Tc_max = CASE.benchmark.theta_core_final_c
    Ts_guess = CASE.benchmark.theta_screen_init_c

    for _ in range(5):
        Rac_90 = cable.Rac_at_temp(Tc_max)
        Rs_guess = cable.sheath_resistance_at_temp(Ts_guess)

        lambda1_prime_90 = cable._lambda1_prime(Rs_guess, Rac_90, cable.sheath_reactance_X())
        lambda1_dprime_90 = cable._lambda1_doubleprime(Rs_guess, Rac_90, Ts_guess)
        lambda1_90 = lambda1_prime_90 + lambda1_dprime_90

        delta_theta = Tc_max - T_amb
        numerator = delta_theta - Wd * (0.5 * T1 + T2 + T3 + T4)

        denominator = Rac_90 * T1 + Rac_90 * (1 + lambda1_90) * T2 + Rac_90 * (1 + lambda1_90 + lambda2) * (T3 + T4)

        I_max = math.sqrt(numerator / denominator)

        Wc_max = (I_max ** 2) * Rac_90
        Ts_guess = Tc_max - (Wc_max + 0.5 * Wd) * T1

    if verbose:
        print("Maximum Current Rating (I) = {0:.2f} A".format(I_max))
        print("Sheath Temp at Max Load    = {0:.3f} C".format(Ts_guess))

    return {
        "Tc_final_c": Tc,
        "Ts_final_c": Ts,
        "Wc_w_per_m": Wc,
        "Ws_w_per_m": Ws,
        "Wd_w_per_m": Wd,
        "I_max_a": I_max,
        "Ts_at_Imax_c": Ts_guess,
        "iterations": converged_iterations,
    }


def run_analytical_solver():
    return solve_case0(verbose=True)


if __name__ == "__main__":
    run_analytical_solver()
