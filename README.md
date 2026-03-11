# FEMTB880Case-0

Python scripts for **CIGRE TB 880 Case 0** (132 kV single-core XLPE, touching trefoil, directly buried), with a centralized data model and an Ansys Mechanical electrothermal iteration workflow.

## Current repository layout (top-level files)

- `tb880_case0_data.py`  
  Centralized Case 0 input data, benchmark values, sources, and explicit model assumptions/conventions used by this repository.

- `cable_losses_tb880_case0.py`  
  IEC 60287-style electrical loss model (conductor, screen/sheath, dielectric) wrapped in a `Cable` class. Uses centralized Case 0 data.

- `analytical_solver_case0.py`  
  Standalone analytical iterative solver for Case 0 using the loss model + IEC thermal resistances. Also exposes `solve_case0()` for programmatic checks.

- `Thermal_Trefoil_Try.py`  
  Ansys Mechanical IronPython script for electrothermal coupling iterations:
  1) solve thermal model, 2) read temperatures, 3) update heat-generation loads from analytical losses, 4) repeat to convergence.

- `temp_convergence_monitor.py`  
  Optional helper for logging/plotting temperature iteration convergence history.

- `test_case0_regression.py`  
  Lightweight regression/self-check entry point to validate that analytical Case 0 outputs remain aligned with benchmark values within tolerances.

## Workflow summary

1. **Centralized data**: `CASE` is defined once in `tb880_case0_data.py`.
2. **Analytical losses**: `cable_losses_tb880_case0.py` computes W/m losses from temperatures and bonding mode.
3. **Mechanical iteration**: `Thermal_Trefoil_Try.py` applies volumetric heat loads (W/m³), solves, reads temperatures, and updates losses.
4. **Convergence monitoring**: optional logging/plotting via `temp_convergence_monitor.py`.
5. **Regression check**: `test_case0_regression.py` quickly guards benchmark behavior after edits.

## Scope and benchmark-specific assumptions

This repository currently targets the **TB880 Case 0 touching-trefoil workflow**.
Some modeling conventions are intentionally benchmark-specific and are kept explicit in centralized data (for example, FE area conventions used to convert dielectric W/m to W/m³ and touching-trefoil spacing assumptions).

## Quick run commands

- Analytical solver:
  ```bash
  python analytical_solver_case0.py
  ```

- Regression/self-check:
  ```bash
  python test_case0_regression.py
  ```

- Ansys Mechanical script:
  Run `Thermal_Trefoil_Try.py` inside Ansys Mechanical (IronPython environment) with required named selections already defined.
