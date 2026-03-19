# -*- coding: utf-8 -*-
"""ANSYS Mechanical IronPython driver for the TB 963 Picard coupling workflow.

I use this script as the main ANSYS Mechanical driver for the fixed-point
Picard electrothermal strategy described in CIGRE TB 963. The algorithm seeds
initial FEM body loads from an ambient-temperature IEC estimate, solves the
thermal FEM model, reads back body temperatures, recomputes IEC 60287 losses,
updates the volumetric heat generation loads, and repeats until the cable core
temperatures converge.
"""
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
        # I append candidate project directories to sys.path because the ANSYS
        # Mechanical scripting host does not reliably start in the project folder,
        # so companion modules must be discoverable explicitly.
        if d and (d not in sys.path) and os.path.isdir(d):
            sys.path.append(d)
    except Exception:
        pass

import cable_losses_tb880_case0 as loss_mod
from tb880_case0_data import CASE


def _reload_module(mod):
    # I reload modules because the ANSYS Mechanical Python interpreter persists
    # between script runs, and stale cached module state would otherwise mask edits.
    try:
        return reload(mod)
    except NameError:
        import importlib
        return importlib.reload(mod)


loss_mod = _reload_module(loss_mod)

try:
    import temp_convergence_monitor as tmon
    tmon = _reload_module(tmon)
except Exception as exc:
    print("Info: temp_convergence_monitor unavailable: {0}".format(exc))
    tmon = None

# ============================================================
# 1) USER SETTINGS
# ============================================================
# I force this setting because ANSYS Mechanical may use m or mm internally, and
# the volumetric heat generation units must match the active model length unit.
# When the model uses mm, W/m^3 must be converted to W/mm^3 by dividing by 1e9.
LENGTH_UNIT = "m"

CABLE_IDS = ["C01", "C02", "C03"]
PARTS_WITH_HEAT = ["Core", "InnerIns", "Screen"]
# I feed back only the core and screen temperatures because those are the state
# variables required by the IEC loss model in the next Picard iteration.
TEMP_REFERENCE_PARTS = ["Core", "Screen"]

# I use the volume-averaged body temperature because that is the most appropriate
# FEM quantity to couple back into the IEC lumped-parameter loss formulation.
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
    # I enforce the NS_ prefix convention so the script can locate Named Selections
    # that identify each FE body in the Mechanical model tree.
    return "NS_{0}_{1}".format(cable_id, part)


def hg_name(cable_id, part):
    # I enforce the HG_ prefix convention so each Internal Heat Generation object
    # can be found or created deterministically for every cable part.
    return "HG_{0}_{1}".format(cable_id, part)


def t_name(cable_id, part):
    # I enforce the T_ prefix convention for Temperature result objects that read
    # back the coupled body temperatures from the FEM solution.
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
    # I convert the IEC-derived volumetric load into the active Mechanical unit
    # system because ANSYS applies body heat generation in model-length units.
    # The upstream conversion from W/m to W/m^3 comes from dividing by FE area,
    # and a mm-based model requires the additional 1e9 scaling to W/mm^3.
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


def _model_length_scale_to_m():
    # I convert the active Mechanical model length unit into a multiplier to meters
    # so that a native 2D body area can always be converted to m^2 before it is
    # stored in CASE or used in the W/m to W/m^3 source conversion.
    unit = LENGTH_UNIT.strip().lower()
    scale_map = {
        "m": 1.0,
        "mm": 1e-3,
    }
    # I fail immediately for any unsupported unit string because a wrong area-unit
    # assumption would silently corrupt every volumetric heat generation value.
    if unit not in scale_map:
        raise RuntimeError("Unsupported LENGTH_UNIT for area conversion: '{0}'".format(LENGTH_UNIT))
    return scale_map[unit]


def _selection_entity_ids(ns_obj):
    # I collect entity ids from the Named Selection object using the most common
    # Mechanical access paths so the helper remains compatible with different
    # ACT / IronPython object layouts seen across ANSYS versions.
    candidates = []
    location = getattr(ns_obj, "Location", None)
    if location is not None and hasattr(location, "Ids"):
        candidates.append(location.Ids)
    if hasattr(ns_obj, "Ids"):
        candidates.append(ns_obj.Ids)

    # I return the first non-empty id list because the heat-load Named Selection
    # is expected to resolve to a concrete scoped set of model entities.
    for ids in candidates:
        if ids is None:
            continue
        entity_ids = [int(entity_id) for entity_id in ids]
        if entity_ids:
            return entity_ids

    # I stop with a clear error instead of guessing because the FEM area split must
    # be traceable to actual model bodies, not inferred from incomplete selection data.
    raise RuntimeError(
        "Named Selection '{0}' has no accessible entity ids for FE body area readout.".format(ns_obj.Name)
    )


def _geo_entity_by_id(entity_id):
    # I resolve each entity id back to a geometry object so the script can query the
    # actual body area represented by the Named Selection rather than using stored
    # reference diameters. ANSYS APIs vary, so I check both common access names.
    geo_data = ExtAPI.DataModel.GeoData
    if hasattr(geo_data, "GeoEntityById"):
        return geo_data.GeoEntityById(entity_id)
    if hasattr(geo_data, "GetGeoEntityById"):
        return geo_data.GetGeoEntityById(entity_id)
    raise RuntimeError("Unable to resolve geometry entities by id in this Mechanical session.")


def _as_float_scalar(value_obj, context_text):
    # ANSYS objects sometimes return plain numeric values and sometimes return value
    # wrappers with a .Value field. I normalize both cases to a Python float so the
    # downstream area summation code can stay explicit and version-tolerant.
    if hasattr(value_obj, "Value"):
        return float(value_obj.Value)
    return float(value_obj)


def _body_area_native_units(geo_entity, ns_name_text):
    # I verify that each resolved entity is a body because this helper is intended
    # only for 2D body areas that receive Internal Heat Generation, not for edges,
    # faces, or convection boundaries.
    type_name = ""
    try:
        type_name = geo_entity.GetType().Name
    except Exception:
        type_name = str(type(geo_entity))
    if "Body" not in type_name:
        raise RuntimeError(
            "Named Selection '{0}' includes non-body entity type '{1}'. Expected body-scoped heat regions.".format(
                ns_name_text, type_name
            )
        )

    # I then read the native body area using the common ANSYS property / method names.
    # The value is still in native model units here and is converted to m^2 later.
    for attr in ("Area", "area"):
        if hasattr(geo_entity, attr):
            area_value = getattr(geo_entity, attr)
            return _as_float_scalar(area_value, "area attribute")
    for method_name in ("GetArea", "getArea"):
        if hasattr(geo_entity, method_name):
            return _as_float_scalar(getattr(geo_entity, method_name)(), "area method")

    # I fail loudly if the area accessor is unavailable because using an uncertain
    # fallback here would break the thesis traceability requirement for heat loading.
    raise RuntimeError(
        "Unable to read body area for Named Selection '{0}' from geometry entity '{1}'.".format(
            ns_name_text, type_name
        )
    )


def get_named_selection_body_area_m2(ns_obj):
    # I sum the body areas represented by the Named Selection because a single heat
    # region may legitimately contain more than one FE body in the Mechanical model.
    area_native_total = 0.0
    entity_ids = _selection_entity_ids(ns_obj)
    for entity_id in entity_ids:
        geo_entity = _geo_entity_by_id(entity_id)
        area_native_total += _body_area_native_units(geo_entity, ns_obj.Name)

    # Zero or negative total area is treated as a hard error because it means the
    # selected heat-load region cannot safely support W/m to W/m^3 conversion.
    if area_native_total <= 0.0:
        raise RuntimeError(
            "Named Selection '{0}' resolved to zero body area. Check the scoped FE bodies.".format(ns_obj.Name)
        )

    # I convert the summed native area to m^2 so the runtime store and the cable FE
    # model fields are always kept in the same SI unit system as the reference areas.
    length_scale_to_m = _model_length_scale_to_m()
    return area_native_total * (length_scale_to_m ** 2)


def collect_model_read_fe_areas_m2(ns_map):
    # I collect only the regions that actually receive internal heat generation so
    # the runtime area map stays focused on the FEM source-conversion workflow.
    area_map = {}
    for key, ns_obj in ns_map.items():
        cid, part = key
        if part not in PARTS_WITH_HEAT:
            continue
        area_map[(cid, part)] = get_named_selection_body_area_m2(ns_obj)
    return area_map


def _assign_model_read_fe_areas_to_cables(cables, area_map):
    # I keep an explicit mapping from Mechanical region names to the cable FE-model
    # attributes so the geometry-reference fields remain untouched and readable.
    attr_by_part = {
        "Core": "A_cond_FE_model",
        "InnerIns": "A_innerins_FE_model",
        "Screen": "A_screen_FE_model",
    }
    geom_attr_by_part = {
        "Core": "A_cond_FE_geom",
        "InnerIns": "A_innerins_FE_geom",
        "Screen": "A_screen_FE_geom",
    }

    # I store the model-read areas centrally on CASE for traceability and possible
    # reuse elsewhere in the Mechanical workflow.
    CASE.runtime_fe_areas.a_model_read_m2_by_region = dict(area_map)

    # I overwrite each cable object's FEM model areas with the actual ANSYS body area
    # that will receive the Internal Heat Generation load, then print a clear startup
    # comparison against the geometry-derived reference area for thesis documentation.
    for cid in CABLE_IDS:
        cable = cables[cid]
        for part in PARTS_WITH_HEAT:
            key = (cid, part)
            if key not in area_map:
                raise RuntimeError("Missing model-read FE area for {0} {1}.".format(cid, part))
            model_area_m2 = area_map[key]
            setattr(cable, attr_by_part[part], model_area_m2)
            geom_area_m2 = getattr(cable, geom_attr_by_part[part])
            rel_diff_pct = 100.0 * (model_area_m2 - geom_area_m2) / geom_area_m2
            print(
                "FE area {0} {1}: geom={2:.9e} m^2, model={3:.9e} m^2, rel_diff={4:+.3f}%".format(
                    cid, part, geom_area_m2, model_area_m2, rel_diff_pct
                )
            )


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

model_read_fe_areas_m2 = collect_model_read_fe_areas_m2(ns_map)
_assign_model_read_fe_areas_to_cables(cables, model_read_fe_areas_m2)

hg_map = {}
t_map = {}
for cid in CABLE_IDS:
    for part in TEMP_REFERENCE_PARTS:
        t_map[(cid, part)] = ensure_temperature_result(solution, t_name(cid, part), ns_map[(cid, part)])
    for part in PARTS_WITH_HEAT:
        hg_map[(cid, part)] = ensure_internal_heatgen(analysis, hg_name(cid, part), ns_map[(cid, part)])

# I seed the first FEM solve with ambient-temperature IEC losses for both the core
# and the screen, then let the Picard loop progressively correct those source terms.
for cid in CABLE_IDS:
    cable = cables[cid]
    losses = cable.calculate_losses(T_guess_C, T_guess_C)
    set_internal_heat_generation(hg_map[(cid, "Core")], losses["core"] / cable.A_cond_FE_model)
    set_internal_heat_generation(hg_map[(cid, "InnerIns")], losses["dielectric"] / cable.A_innerins_FE_model)
    set_internal_heat_generation(hg_map[(cid, "Screen")], losses["screen"] / cable.A_screen_FE_model)

print("Initialised loads with T_guess = {0} C (Core+Screen readback seeds).".format(T_guess_C))


def solve_and_eval():
    solution.Solve(True)
    solution.EvaluateAllResults()


T_prev_core = {}
converged = False

for it in range(1, max_iter + 1):
    # I follow the TB 963 fixed-point sequence: solve the thermal FEM field, read
    # back temperatures, recompute IEC losses, update heat loads, and test convergence.
    solve_and_eval()

    T_curr_core = {}
    T_curr_screen = {}
    for cid in CABLE_IDS:
        T_curr_core[cid] = read_temperature_celsius(t_map[(cid, "Core")])
        T_curr_screen[cid] = read_temperature_celsius(t_map[(cid, "Screen")])

    max_dT = 0.0
    if T_prev_core:
        for cid in CABLE_IDS:
            # I check convergence on the maximum absolute change in cable core
            # temperature because TB 963 recommends monitoring the coupled state
            # variable that drives the next IEC resistance update.
            dT = abs(T_curr_core[cid] - T_prev_core.get(cid, T_curr_core[cid]))
            if dT > max_dT:
                max_dT = dT

    for cid in CABLE_IDS:
        cable = cables[cid]
        losses = cable.calculate_losses(T_curr_core[cid], T_curr_screen[cid])
        try:
            set_internal_heat_generation(hg_map[(cid, "Core")], losses["core"] / cable.A_cond_FE_model)
            set_internal_heat_generation(hg_map[(cid, "InnerIns")], losses["dielectric"] / cable.A_innerins_FE_model)
            set_internal_heat_generation(hg_map[(cid, "Screen")], losses["screen"] / cable.A_screen_FE_model)
        except Exception as exc:
            raise RuntimeError("Failed applying heat generation for {0}: {1}".format(cid, exc))

    msg = "Iter {0:02d}: ".format(it)
    for cid in CABLE_IDS:
        msg += "{0} Tc={1:.2f} Ts={2:.2f} | ".format(cid, T_curr_core[cid], T_curr_screen[cid])
    msg += "max_dT={0:.4f}".format(max_dT)
    print(msg)

    if mon is not None:
        mon.log(it, "iter", T_curr_core, max_dT)

    # I stop when the maximum core-temperature update falls below the Picard
    # tolerance, which is the practical convergence check recommended in TB 963.
    if T_prev_core and max_dT < tol_dT:
        # I run one additional FEM solve with the converged heat loads so that the
        # reported field is fully consistent with the final IEC loss update.
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
    # If convergence is not reached, I still report the last updated solution but I
    # flag it clearly so the user can inspect mesh quality or the iteration settings.
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
