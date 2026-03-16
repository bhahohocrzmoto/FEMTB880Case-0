# ============================================================
# Thermal_Trefoil_Try.py
# ------------------------------------------------------------
# ANSYS Mechanical IronPython main script.
#
# Outer electrothermal Picard iteration:
#   1) Solve thermal FEM with current volumetric heat loads,
#   2) Read core/screen temperatures (degC),
#   3) Recompute IEC losses from readback temperatures,
#   4) Update volumetric loads and repeat to convergence.
# ============================================================

import os
import sys
import importlib

SCRIPT_DIRS = []
try:
    SCRIPT_DIRS.append(os.path.dirname(__file__))
except Exception:
    pass
try:
    SCRIPT_DIRS.append(ExtAPI.DataModel.Project.ProjectDirectory)
except Exception:
    pass
try:
    SCRIPT_DIRS.append(os.getcwd())
except Exception:
    pass

for d in SCRIPT_DIRS:
    try:
        if d and (d not in sys.path) and os.path.isdir(d):
            sys.path.append(d)
    except Exception:
        pass

import cable_losses_tb880_case0 as loss_mod
from tb880_case0_data import CASE
loss_mod = importlib.reload(loss_mod)

try:
    import temp_convergence_monitor as tmon
    tmon = importlib.reload(tmon)
except Exception as exc:
    print("Info: temp_convergence_monitor unavailable: {0}".format(exc))
    tmon = None

# ============================================================
# 1) USER SETTINGS
# ============================================================
# Set explicitly to match Mechanical model length units.
# No implicit fallback: use only "m" or "mm".
LENGTH_UNIT = "m"

CABLE_IDS = ["C01", "C02", "C03"]
PARTS_WITH_HEAT = ["Core", "InnerIns", "Screen"]
# Core and Screen temperatures are fed back into IEC loss update in next outer iteration.
TEMP_REFERENCE_PARTS = ["Core", "Screen"]

TEMP_READ_MODE = "average"   # "average", "maximum", "minimum"

T_guess_C = CASE.installation.ambient_temp_c
max_iter = 20
tol_dT = 1e-3

I_RMS_A = CASE.benchmark.i_final_a
bonding = CASE.installation.bonding

cables = loss_mod.create_case0_cables(I_rms_A=I_RMS_A, bonding=bonding)

if LENGTH_UNIT not in ("m", "mm"):
    raise ValueError("LENGTH_UNIT must be explicitly set to 'm' or 'mm'. Got: {0}".format(LENGTH_UNIT))

# ============================================================
# 2) OPTIONAL CONVERGENCE MONITOR
# ============================================================
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
        live_plot=True,
    )


# ============================================================
# 3) HELPERS
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
        except Exception:
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
        children = getattr(node, "Children", None)
        if children is None:
            continue
        for ch in children:
            try:
                if ch.Name == name and (type_contains in ch.GetType().Name):
                    return ch
            except Exception:
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
    tr.Location = location_ns
    return tr


def _extract_unit_string(value_obj):
    for attr in ("Unit", "Units"):
        if hasattr(value_obj, attr):
            unit_obj = getattr(value_obj, attr)
            try:
                text = str(unit_obj)
                if text:
                    return text
            except Exception:
                pass
    return ""


def _to_celsius(value_obj, result_name):
    raw_value = float(value_obj.Value)
    unit_text = _extract_unit_string(value_obj).strip().lower().replace(chr(176), "")

    if unit_text in ("", "c", "degc", "celsius", "degree c"):
        return raw_value
    if unit_text in ("k", "kelvin"):
        return raw_value - 273.15

    raise RuntimeError(
        "Unsupported temperature unit for '{0}': '{1}'. "
        "Expected Celsius or Kelvin. Aborting to avoid unit-misuse.".format(result_name, unit_text)
    )


def read_temperature_celsius(temp_result, mode=TEMP_READ_MODE):
    read_map = {
        "average": lambda r: r.Average,
        "maximum": lambda r: r.Maximum,
        "minimum": lambda r: r.Minimum,
    }
    mode_key = str(mode).strip().lower()
    if mode_key not in read_map:
        raise Exception("Unsupported TEMP_READ_MODE='{0}'. Use average, maximum, or minimum.".format(mode))

    value_obj = read_map[mode_key](temp_result)
    return _to_celsius(value_obj, temp_result.Name)


# ============================================================
# 4) MAIN ITERATION
# ============================================================
analysis = Model.Analyses[0]
solution = analysis.Solution

ns_map = {}
missing = []
for cid in CABLE_IDS:
    for part in list(set(PARTS_WITH_HEAT + TEMP_REFERENCE_PARTS)):
        key = (cid, part)
        name = ns_name(cid, part)
        try:
            ns_map[key] = get_named_selection(name)
        except Exception:
            missing.append(name)

if missing:
    uniq = []
    for m in missing:
        if m not in uniq:
            uniq.append(m)
    raise Exception("Create these Named Selections first:\n  " + "\n  ".join(uniq))

hg_map = {}
t_map = {}
for cid in CABLE_IDS:
    for part in TEMP_REFERENCE_PARTS:
        t_map[(cid, part)] = ensure_temperature_result(solution, t_name(cid, part), ns_map[(cid, part)])
    for part in PARTS_WITH_HEAT:
        hg_map[(cid, part)] = ensure_internal_heatgen(analysis, hg_name(cid, part), ns_map[(cid, part)])

# Initial IEC losses from ambient guess, then W/m -> W/m^3 for FEM loads.
for cid in CABLE_IDS:
    cable = cables[cid]
    losses = cable.calculate_losses(T_guess_C, T_guess_C)
    set_internal_heat_generation(hg_map[(cid, "Core")], losses["core"] / cable.A_cond_FE)
    set_internal_heat_generation(hg_map[(cid, "InnerIns")], losses["dielectric"] / cable.A_innerins_FE)
    set_internal_heat_generation(hg_map[(cid, "Screen")], losses["screen"] / cable.A_screen_FE)

print("Initialised loads with T_guess = {0} C (Core+Screen readback seeds).".format(T_guess_C))


def solve_and_eval():
    solution.Solve(True)
    solution.EvaluateAllResults()


T_prev_core = {}
converged = False

for it in range(1, max_iter + 1):
    solve_and_eval()

    T_curr_core = {}
    T_curr_screen = {}
    for cid in CABLE_IDS:
        T_curr_core[cid] = read_temperature_celsius(t_map[(cid, "Core")])
        T_curr_screen[cid] = read_temperature_celsius(t_map[(cid, "Screen")])

    max_dT = 0.0
    if T_prev_core:
        for cid in CABLE_IDS:
            dT = abs(T_curr_core[cid] - T_prev_core.get(cid, T_curr_core[cid]))
            if dT > max_dT:
                max_dT = dT

    for cid in CABLE_IDS:
        cable = cables[cid]
        losses = cable.calculate_losses(T_curr_core[cid], T_curr_screen[cid])
        try:
            set_internal_heat_generation(hg_map[(cid, "Core")], losses["core"] / cable.A_cond_FE)
            set_internal_heat_generation(hg_map[(cid, "InnerIns")], losses["dielectric"] / cable.A_innerins_FE)
            set_internal_heat_generation(hg_map[(cid, "Screen")], losses["screen"] / cable.A_screen_FE)
        except Exception as exc:
            raise RuntimeError("Failed applying heat generation for {0}: {1}".format(cid, exc))

    msg = "Iter {0:02d}: ".format(it)
    for cid in CABLE_IDS:
        msg += "{0} Tc={1:.2f} Ts={2:.2f} | ".format(cid, T_curr_core[cid], T_curr_screen[cid])
    msg += "max_dT={0:.4f}".format(max_dT)
    print(msg)

    if mon is not None:
        mon.log(it, "iter", T_curr_core, max_dT)

    if T_prev_core and max_dT < tol_dT:
        print("Converged -> running FINAL solve with updated losses...")
        solve_and_eval()
        T_final_core = {}
        T_final_screen = {}
        for cid in CABLE_IDS:
            T_final_core[cid] = read_temperature_celsius(t_map[(cid, "Core")])
            T_final_screen[cid] = read_temperature_celsius(t_map[(cid, "Screen")])
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
        T_final_core[cid] = read_temperature_celsius(t_map[(cid, "Core")])
        T_final_screen[cid] = read_temperature_celsius(t_map[(cid, "Screen")])
    msg = "Final (Non-Converged): "
    for cid in CABLE_IDS:
        msg += "{0} Tc={1:.4f} Ts={2:.4f} | ".format(cid, T_final_core[cid], T_final_screen[cid])
    print(msg)
    if mon is not None:
        mon.log(it, "final", T_final_core, 0.0)
        mon.save_png("temp_history.png")

print("Done.")
