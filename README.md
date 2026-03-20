# FEMTB880Case-0

Python scripts for **CIGRE TB 880 Case #0** (132 kV single-core XLPE, touching trefoil, directly buried), with centralized benchmark data and an ANSYS Mechanical outer electrothermal iteration workflow.

## Default benchmark interpretation (repository default)

This repository defaults to the **TB 880 Case #0-1 benchmark-faithful IEC base case**.
Sheath-loss inclusion is now controlled centrally in `tb880_case0_data.py` by two
boolean flags on `CASE.assumptions`:

- `include_sheath_circulating_losses = True`
- `include_sheath_eddy_losses = False`

The default benchmark therefore keeps circulating sheath losses included while
excluding eddy-current sheath losses.

Stored benchmark targets remain:

- `I_final = 821.7763334392 A`
- `Wc_final = 26.6895 W/m`
- `Ws_final = 7.8442 W/m`
- `theta_core_final = 90.0 C`
- `theta_screen_final = 78.7130 C`

To study alternate sheath-loss combinations, update those two central flags before running the analytical or FEM workflows.

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

- Conductor DC resistance at 20°C is treated as a **primary benchmark input** (`r_cond_dc_20_ohm_per_m`) for IEC/TB880 traceability.
- Material resistivity is still stored for documentation and related calculations, but benchmark DC conductor resistance is not reconstructed from `rho/A` in the default path.
- FE area convention for `InnerIns` remains explicitly documented and unchanged.
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
