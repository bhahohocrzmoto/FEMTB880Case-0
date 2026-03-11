# ============================================================
# Thermal_Trefoil_Try.py (REVISED)
# ------------------------------------------------------------
# ANSYS Mechanical IronPython main script
#
# Uses object oriented Cable class from cable_losses_tb880_case0.
# Handles:
#   - finding named selections / loads / temperature results
#   - running outer iteration: Solve -> Read T -> Update loads
#   - convergence monitoring (unchanged)
# ============================================================

import os, sys, math

# ------------------------------------------------------------
# 0) IMPORT LOSS MODULE (must be in same folder as this script)
# ------------------------------------------------------------
SCRIPT_DIRS = []
try:
    SCRIPT_DIRS.append(os.path.dirname(__file__))
except:
    pass
try:
    SCRIPT_DIRS.append(ExtAPI.DataModel.Project.ProjectDirectory)
except:
    pass
try:
    SCRIPT_DIRS.append(os.getcwd())
except:
    pass

for d in SCRIPT_DIRS:
    try:
        if d and (d not in sys.path) and os.path.isdir(d):
            sys.path.append(d)
    except:
        pass

# Import the revised loss module
import cable_losses_tb880_case0 as loss_mod
from tb880_case0_data import CASE
try:
    import imp
    imp.reload(loss_mod)
except:
    pass

# Optional convergence monitor
try:
    import temp_convergence_monitor as tmon
    try:
        imp.reload(tmon)
    except:
        pass
except:
    tmon = None

# ============================================================
# 1) USER SETTINGS (can be adjusted here)
# ============================================================
LENGTH_UNIT = "m"        # "m" or "mm" affects Quantity unit

CABLE_IDS = ["C01", "C02", "C03"]
PARTS_WITH_HEAT = ["Core", "InnerIns", "Screen"]
TEMP_REFERENCE_PARTS = ["Core", "Screen"]   # parts whose temperature we read

# Electrothermal coupling option: which scalar from FE temperature results is fed back to losses.
TEMP_READ_MODE = "average"   # supported: "average", "maximum", "minimum"

# Outer iteration parameters
T_guess_C = CASE.installation.ambient_temp_c
max_iter  = 20
tol_dT    = 1e-3          # convergence tolerance on core temperature change

# Current (same for all phases)
I_RMS_A = CASE.benchmark.i_final_a
bonding = CASE.installation.bonding          # "solid" (both ends), "single", or "cross"

# Create the three Cable objects using the helper function
cables = loss_mod.create_case0_cables(I_rms_A=I_RMS_A, bonding=bonding)

# ------------------------------------------------------------
# 2) SET UP CONVERGENCE MONITOR (optional, unchanged)
# ------------------------------------------------------------
mon = None
if tmon is not None:
    base_dir = None
    for d in SCRIPT_DIRS:
        if d and os.path.isdir(d):
            base_dir = d
            break
    if base_dir is None:
        base_dir = os.getcwd()
    mon = tmon.TempConvergenceMonitor(
        base_dir,
        CABLE_IDS,
        csv_name="temp_history.csv",
        clear_csv=False,
        live_plot=True
    )

# ============================================================
# 3) HELPER FUNCTIONS (largely unchanged, but adapted to use cables dict)
# ============================================================
def ns_name(cable_id, part):
    return "NS_{0}_{1}".format(cable_id, part)

def hg_name(cable_id, part):
    return "HG_{0}_{1}".format(cable_id, part)

def t_name(cable_id, part):
    return "T_{0}_{1}".format(cable_id, part)

def find_by_name_and_type(name, type_contains):
    objs = ExtAPI.DataModel.GetObjectsByName(name)
    for o in objs:
        try:
            t = o.GetType().Name
        except:
            t = ""
        if type_contains in t:
            return o
    return None

def get_named_selection(name):
    ns = find_by_name_and_type(name, "NamedSelection")
    if ns is None:
        raise Exception("Missing Named Selection: '{0}'".format(name))
    return ns

def q_to_model_units(q_m3):
    if LENGTH_UNIT.lower() == "mm":
        return (q_m3 / 1e9, "W mm^-1 mm^-1 mm^-1")
    else:
        return (q_m3, "W m^-1 m^-1 m^-1")

def set_internal_heat_generation(hg_obj, q_m3):
    q_val, q_unit = q_to_model_units(q_m3)
    hg_obj.Magnitude.Output.SetDiscreteValue(0, Quantity(q_val, q_unit))

def ensure_internal_heatgen(analysis, name, location_ns):
    hg = find_by_name_and_type(name, "InternalHeatGeneration")
    if hg is None:
        hg = analysis.AddInternalHeatGeneration()
        hg.Name = name
    hg.Location = location_ns
    return hg

def find_result_in_solution(solution, name, type_contains):
    stack = [solution]
    while stack:
        node = stack.pop()
        try:
            children = node.Children
        except:
            children = None
        if children is None:
            continue
        for ch in children:
            try:
                if ch.Name == name and (type_contains in ch.GetType().Name):
                    return ch
            except:
                pass
            stack.append(ch)
    return None

def ensure_temperature_result(solution, name, location_ns):
    tr = find_result_in_solution(solution, name, "Temperature")
    if tr is None:
        tr = solution.AddTemperature()
        tr.Name = name
    try:
        tr.ClearGeneratedData()
    except Exception as e:
        print("Warning: could not clear generated data for result '{0}': {1}".format(name, e))
    try:
        tr.Location = location_ns
    except:
        pass
    return tr

def read_temperature_value(temp_result, mode=TEMP_READ_MODE):
    read_map = {
        "average": lambda r: float(r.Average.Value),
        "maximum": lambda r: float(r.Maximum.Value),
        "minimum": lambda r: float(r.Minimum.Value),
    }
    mode_key = str(mode).strip().lower()
    if mode_key not in read_map:
        raise Exception("Unsupported TEMP_READ_MODE='{0}'. Use average, maximum, or minimum.".format(mode))

    # First honor the configured read mode; only fall back if that field is unavailable on this result.
    try:
        return read_map[mode_key](temp_result)
    except:
        pass

    for fb in [m for m in ["average", "maximum", "minimum"] if m != mode_key]:
        try:
            return read_map[fb](temp_result)
        except:
            pass
    raise Exception("Could not read temperature from '{0}'".format(temp_result.Name))

# ============================================================
# 4) MAIN ITERATION
# ============================================================
analysis = Model.Analyses[0]
solution = analysis.Solution

# ---- 4.1 Collect named selections ----
ns_map = {}
missing = []
for cid in CABLE_IDS:
    all_needed_parts = list(set(PARTS_WITH_HEAT + TEMP_REFERENCE_PARTS))
    for part in all_needed_parts:
        key = (cid, part)
        name = ns_name(cid, part)
        try:
            ns_map[key] = get_named_selection(name)
        except:
            missing.append(name)

uniq = []
for m in missing:
    if m not in uniq:
        uniq.append(m)
if uniq:
    raise Exception("Create these Named Selections first:\n  " + "\n  ".join(uniq))

# ---- 4.2 Ensure loads and temperature results exist ----
hg_map = {}   # (cid, part) -> internal heat generation object
t_map  = {}   # (cid, part) -> temperature result object

for cid in CABLE_IDS:
    # Temperature results for core and screen
    for part in TEMP_REFERENCE_PARTS:
        t_map[(cid, part)] = ensure_temperature_result(
            solution,
            t_name(cid, part),
            ns_map[(cid, part)]
        )
    # Heat generation loads
    for part in PARTS_WITH_HEAT:
        hg_map[(cid, part)] = ensure_internal_heatgen(
            analysis,
            hg_name(cid, part),
            ns_map[(cid, part)]
        )

# ---- 4.3 Initialise loads with guess temperatures ----
for cid in CABLE_IDS:
    cable = cables[cid]
    # Use guess temperature for both core and screen (first iteration)
    losses = cable.calculate_losses(T_guess_C, T_guess_C)
    set_internal_heat_generation(hg_map[(cid, "Core")],      losses["core"] / cable.A_cond_FE)
    # Dielectric W/m -> W/m^3 uses CASE.assumptions.innerins_area_convention via cable.A_innerins_FE.
    set_internal_heat_generation(hg_map[(cid, "InnerIns")],  losses["dielectric"] / cable.A_innerins_FE)
    set_internal_heat_generation(hg_map[(cid, "Screen")],    losses["screen"] / cable.A_screen_FE)

print("Initialised loads with T_guess = {0} C (applied to Core & Screen)".format(T_guess_C))

# ---- 4.4 Iteration loop ----
T_prev_core = {}
converged = False

def solve_and_eval():
    solution.Solve(True)
    try:
        solution.EvaluateAllResults()
    except:
        pass

for it in range(1, max_iter + 1):
    solve_and_eval()

    # Read temperatures
    T_curr_core = {}
    T_curr_screen = {}
    for cid in CABLE_IDS:
        T_curr_core[cid]   = read_temperature_value(t_map[(cid, "Core")])
        T_curr_screen[cid] = read_temperature_value(t_map[(cid, "Screen")])

    # Compute max core temperature change
    max_dT = 0.0
    if T_prev_core:
        for cid in CABLE_IDS:
            dT = abs(T_curr_core[cid] - T_prev_core.get(cid, T_curr_core[cid]))
            if dT > max_dT:
                max_dT = dT

    # Update loads using the Cable objects
    for cid in CABLE_IDS:
        cable = cables[cid]
        losses = cable.calculate_losses(T_curr_core[cid], T_curr_screen[cid])
        set_internal_heat_generation(hg_map[(cid, "Core")],      losses["core"] / cable.A_cond_FE)
        # Dielectric W/m -> W/m^3 uses CASE.assumptions.innerins_area_convention via cable.A_innerins_FE.
        set_internal_heat_generation(hg_map[(cid, "InnerIns")],  losses["dielectric"] / cable.A_innerins_FE)
        set_internal_heat_generation(hg_map[(cid, "Screen")],    losses["screen"] / cable.A_screen_FE)

    # Logging
    msg = "Iter {0:02d}: ".format(it)
    for cid in CABLE_IDS:
        msg += "{0} Tc={1:.2f} Ts={2:.2f} | ".format(cid, T_curr_core[cid], T_curr_screen[cid])
    msg += "max_dT={0:.4f}".format(max_dT)
    print(msg)

    if mon is not None:
        mon.log(it, "iter", T_curr_core, max_dT)

    # Check convergence
    if T_prev_core and max_dT < tol_dT:
        print("Converged -> running FINAL solve with updated losses...")
        solve_and_eval()
        # Read final temperatures
        T_final_core = {}
        T_final_screen = {}
        for cid in CABLE_IDS:
            T_final_core[cid]   = read_temperature_value(t_map[(cid, "Core")])
            T_final_screen[cid] = read_temperature_value(t_map[(cid, "Screen")])
        # Print final
        msg = "Final: "
        for cid in CABLE_IDS:
            msg += "{0} Tc={1:.4f} Ts={2:.4f} | ".format(cid, T_final_core[cid], T_final_screen[cid])
        print(msg)
        if mon is not None:
            mon.log(it, "final", T_final_core, 0.0)
            mon.save_png("temp_history.png")
        converged = True
        break

    T_prev_core = T_curr_core

if not converged:
    print("Reached max_iter without convergence -> running FINAL solve with last updated losses...")
    solve_and_eval()
    T_final_core = {}
    T_final_screen = {}
    for cid in CABLE_IDS:
        T_final_core[cid]   = read_temperature_value(t_map[(cid, "Core")])
        T_final_screen[cid] = read_temperature_value(t_map[(cid, "Screen")])
    msg = "Final (Non-Converged): "
    for cid in CABLE_IDS:
        msg += "{0} Tc={1:.4f} Ts={2:.4f} | ".format(cid, T_final_core[cid], T_final_screen[cid])
    print(msg)
    if mon is not None:
        mon.log(it, "final", T_final_core, 0.0)
        mon.save_png("temp_history.png")

print("Done.")
