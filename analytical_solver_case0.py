# ============================================================\n
# analytical_solver_case0.py
# ------------------------------------------------------------
# Standalone analytical iterative solver for TB 880 Case 0.
# Uses the Cable class to calculate temperature-dependent losses
# and pushes them through the IEC thermal resistance network.
# ============================================================\n

import math
# Import your revised loss module
import cable_losses_tb880_case0 as loss_mod

def run_analytical_solver():
    # --------------------------------------------------------
    # 1. INPUT PARAMETERS (From Thesis Section 2.2.1)
    # --------------------------------------------------------
    # TODO: Replace these with the exact values calculated in your thesis!
    T1 = 0.12   # Thermal resistance of insulation [K.m/W]
    T2 = 0.00   # Bedding (0 for Case 0 as sheath is directly over screen)
    T3 = 0.04   # Oversheath [K.m/W]
    T4 = 0.80   # External/Soil thermal resistance for the hottest cable [K.m/W]
    
    T_amb = 20.0       # Ambient soil temperature [deg C]
    I_load = 900.0     # Benchmark load current [A] (Adjust to your Case 0 load)

    # Initialize the Cable (using the central cable C02 as the hottest)
    cable = loss_mod.cables["C02"]

    # --------------------------------------------------------
    # 2. ITERATIVE TEMPERATURE SOLVER
    # --------------------------------------------------------
    print(f"--- IEC 60287 Iterative Temperature Solver ---")
    print(f"Load Current: {I_load} A | Ambient Temp: {T_amb} C")
    
    # Initial guesses
    Tc = T_amb
    Ts = T_amb
    
    max_iter = 100
    tolerance = 1e-4
    
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
        # Drop across insulation
        delta_T_ins = (Wc + 0.5 * Wd) * T1
        
        # Drop across oversheath (T2 is 0 here)
        delta_T_sheath = (Wc + Ws + Wd) * (T2 + T3)
        
        # Drop across soil
        delta_T_soil = (Wc + Ws + Wd) * T4
        
        # d. Update Temperatures
        Tc_new = T_amb + delta_T_ins + delta_T_sheath + delta_T_soil
        Ts_new = Tc_new - delta_T_ins  # Sheath is cooler than core by the drop across T1
        
        # e. Check Convergence
        error_Tc = abs(Tc_new - Tc)
        
        Tc = Tc_new
        Ts = Ts_new
        
        if error_Tc < tolerance:
            print(f"Converged in {iteration + 1} iterations.")
            print(f"Final Conductor Temp (Tc) = {Tc:.3f} C")
            print(f"Final Sheath Temp (Ts)    = {Ts:.3f} C")
            print(f"Losses: Wc = {Wc:.2f} W/m, Wd = {Wd:.2f} W/m, Ws = {Ws:.2f} W/m\n")
            break

    # --------------------------------------------------------
    # 3. ANALYTICAL AMPACITY (CURRENT RATING)
    # --------------------------------------------------------
    print(f"--- IEC 60287 Analytical Ampacity (For Tc = 90 C) ---")
    
    Tc_max = 90.0
    # We must guess Ts to find Rs, or iterate. A simple 2-step iteration is enough.
    Ts_guess = 80.0 
    
    for _ in range(5):
        Rac_90 = cable.Rac_at_temp(Tc_max)
        Rs_guess = cable.sheath_resistance_at_temp(Ts_guess)
        
        lambda1_prime_90 = cable._lambda1_prime(Rs_guess, Rac_90, cable.sheath_reactance_X())
        lambda1_dprime_90 = cable._lambda1_doubleprime(Rs_guess, Rac_90, Ts_guess)
        lambda1_90 = lambda1_prime_90 + lambda1_dprime_90
        
        # IEC Eq 1 Ampacity Numerator
        delta_theta = Tc_max - T_amb
        numerator = delta_theta - Wd * (0.5 * T1 + T2 + T3 + T4)
        
        # IEC Eq 1 Ampacity Denominator
        denominator = Rac_90 * T1 + Rac_90 * (1 + lambda1_90) * T2 + Rac_90 * (1 + lambda1_90 + lambda2) * (T3 + T4)
        
        I_max = math.sqrt(numerator / denominator)
        
        # Update Ts_guess based on this new current for the next loop
        Wc_max = (I_max ** 2) * Rac_90
        Ts_guess = Tc_max - (Wc_max + 0.5 * Wd) * T1

    print(f"Maximum Current Rating (I) = {I_max:.2f} A")
    print(f"Sheath Temp at Max Load    = {Ts_guess:.3f} C")

if __name__ == "__main__":
    run_analytical_solver()