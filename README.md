# FEMTB880Case-0

Python scripts for **CIGRE TB 880 Case #0** (132 kV single-core XLPE, touching trefoil, directly buried), with centralized benchmark data and an ANSYS Mechanical outer electrothermal iteration workflow.

## Default benchmark interpretation (repository default)

This repository defaults to the **TB 880 Case #0-1 benchmark-faithful IEC base case**.
Physics simplification flags are controlled centrally in `tb880_case0_data.py`
via boolean flags on `CASE.assumptions`:

**Conductor loss flags:**
- `include_skin_effect = True` (ys factor, IEC 60287-1-1 Eq. 5-8)
- `include_proximity_effect = True` (yp factor, IEC 60287-1-1 Eq. 9-10)

**Dielectric loss flag:**
- `include_dielectric_losses = True` (Wd, IEC 60287-1-1 Eq. 14-15 / TB 880 GP 7)

**Sheath loss flags:**
- `include_sheath_circulating_losses = True` (lambda1_prime, IEC 60287-1-1 Sec 2.3.1)
- `include_sheath_eddy_losses = False` (lambda1_doubleprime, IEC 60287-1-1 Sec 5.3.7.1)
- `include_F_factor_for_eddy_reduction = True` (F factor, IEC 60287-1-1 Sec 2.3.5 / TB 880 GP 31)

All flags default to True except `include_sheath_eddy_losses` which defaults
to False to match the IEC simplified benchmark (TB 880 Case 0-1).

Each flag independently controls one physics simplification in the IEC 60287
loss chain. There is no separate sheath-mode parameter. The user selects the
desired loss assembly directly through these flags.

Stored benchmark targets remain:

- `I_final = 821.7763334392 A`
- `Wc_final = 26.6895 W/m`
- `Ws_final = 7.8442 W/m`
- `theta_core_final = 90.0 C`
- `theta_screen_final = 78.7130 C`

To model different sheath loss scenarios, set the flags as follows:

| Scenario                                    | circ  | eddy  | F flag |
|---------------------------------------------|-------|-------|--------|
| IEC simplified (circulating only, no eddy)  | True  | False | True   |
| TB 880 recommended (both, with F reduction) | True  | True  | True   |
| Eddy-only (no circulating currents)         | False | True  | True   |
| Study: omit F correction                    | True  | True  | False  |
| No sheath losses                            | False | False | True   |

The flags fully control which sheath loss terms are assembled into lambda1.
There is no separate sheath-mode parameter. The user selects the desired
physics directly through these flags.

## Current repository layout (top-level files)

- `tb880_case0_data.py`  
  Centralized Case #0 input data, benchmark values, source traceability, and explicit model assumptions.

- `cable_losses_tb880_case0.py`  
  IEC 60287-style electrical loss model in a `Cable` class. `calculate_losses()` is the authoritative path for conductor/screen/dielectric losses and lambda factors.

- `analytical_solver_case0.py`  
  Standalone analytical iterative solver for Case #0 that uses the common loss API in both thermal iteration and ampacity calculation.

- `Thermal_Trefoil_Try.py`  
  ANSYS Mechanical script for outer electrothermal iteration:
  1) solve thermal FEM,
  2) read temperatures in **Celsius**,
  3) update IEC losses,
  4) convert W/m -> W/m^3,
  5) re-apply FEM heat loads,
  6) iterate to convergence.

- `temp_convergence_monitor.py`  
  Optional helper for text-mode CSV logging (Python 3-safe) and optional live convergence plot.

- `test_case0_regression.py`  
  Regression checks for the default benchmark flags, the combined circulating-plus-eddy case, and the zero-sheath-loss case.

## Technical notes

- Conductor DC resistance at 20 degC is treated as a **primary benchmark input** (`r_cond_dc_20_ohm_per_m`) for IEC/TB880 traceability.
- Material resistivity is still stored for documentation and related calculations, but benchmark DC conductor resistance is not reconstructed from `rho/A` in the default path.
- FE area convention for `InnerIns` remains explicitly documented and unchanged.
- The `area_mode` parameter on `calculate_losses()` controls which cross-sectional
  areas are used for all area-dependent quantities in the IEC loss chain (Rs, beta1,
  and the W/m to W/m3 conversion). The analytical solver uses the default
  `area_mode="analytical"`, which routes to the IEC nominal/electrical areas stored
  in the central data file. The FEM driver passes `area_mode="fem"`, which routes
  to the FEM body areas read from ANSYS Mechanical at runtime. This separation allows
  the user to cross-check that the drawn FEM geometry matches the intended analytical
  cable dimensions. Any mismatch between analytical and FEM areas will produce
  different loss values, making geometry errors visible and traceable.
- ANSYS script now guards against model length-unit misuse by requiring explicit `LENGTH_UNIT` validation.

## Quick run commands

- Analytical solver:
  ```bash
  python analytical_solver_case0.py
  ```

- Regression checks:
  ```bash
  python test_case0_regression.py
  ```

- ANSYS Mechanical script:
  Run `Thermal_Trefoil_Try.py` inside ANSYS Mechanical (IronPython environment) with required named selections already defined.
