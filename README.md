# Master's Thesis: CIGRE TB 880 Case 0 Ampacity & FEM Benchmark

Welcome to the official repository for my Master's Thesis project. This repository contains the Python-based codebase designed to facilitate Analytical Ampacity Calculations and Finite Element Method (FEM) automation (via Ansys) for underground power cables. 

## 📖 Project Overview

The primary objective of this thesis is to replicate the **CIGRE TB 880 Case 0** benchmark—a three-phase, 132 kV cable circuit in a touching trefoil formation, directly buried in homogeneous soil. 

Using Ansys Mechanical for the FEM thermal analysis and the analytical loss methods described in the **IEC 60287** standard, this project establishes a "Best Practice" baseline model. From this baseline, the code systematically conducts a sensitivity analysis by selectively deviating from key modeling guidance points outlined in **CIGRE TB 640, TB 880, and TB 963**.

The ultimate goal is to provide engineers with an evidence-based framework that quantifies the trade-offs between temperature accuracy and computational time when making specific modeling choices in thermal FEM analysis.

## 🎯 Key Features & Sensitivity Analysis

The code is structured to automatically generate models, calculate losses, and extract results to evaluate the following modeling choices:

* **Domain Size (TB 963 GP 14):** Quantifies the thermal error caused by using an insufficiently large soil box.
* **Mesh Refinement and Element Order (TB 963 GP 50, 54, 56):** Generates mesh convergence data to evaluate the trade-off between accuracy and solver time (comparing linear vs. quadratic elements and coarse vs. fine meshes).
* **Eddy Current Losses (TB 880 GP 6):** Evaluates the temperature deviation when eddy current losses are neglected, highlighting the importance of including all physical heat sources.
* **Geometric Symmetry (TB 963 GP 46):** Validates the use of symmetry by comparing the results and computational cost of a half-model against a full-geometry model.

## 🚫 Out of Scope

To ensure a highly focused and feasible study, the following features are explicitly excluded from this codebase and the corresponding thesis:
* Transient or cyclic load analysis (steady-state 100% load factor only).
* Fluid dynamics or Computational Fluid Dynamics (CFD); convective effects are handled purely via boundary conditions.
* Fully coupled electromagnetic-thermal analysis (all heat loads/losses are calculated analytically via IEC 60287 and applied as inputs).
* Multiple cable configurations (e.g., ducts, troughs, flat formations); the scope is strictly confined to the direct-buried trefoil case.

## 🗂️ Repository Structure

*(Note: Adjust the folder names below to match your actual repository structure)*

* `docs/`: Documentation, reference materials, and the proposed thesis structure.
* `analytical_iec/`: Python scripts for calculating analytical cable losses according to IEC 60287-1-1.
* `ansys_automation/`: Python scripts (e.g., PyAnsys or APDL generation scripts) used to build the geometry, assign material properties, mesh, and apply boundary conditions in Ansys Mechanical.
* `sensitivity_studies/`: Execution scripts for running the parametric sweeps (Domain Size, Mesh Order, etc.).
* `results_processing/`: Scripts for parsing Ansys output files, calculating temperature deviations, and plotting convergence charts.

## ⚙️ Prerequisites & Setup

To run the scripts in this repository, you will need:
* **Python 3.8+**
* **Ansys Mechanical** (A valid license is required to run the FEM solver)
* **Required Python Libraries:** `numpy`, `pandas`, `matplotlib`, `scipy` *(and `ansys-dpf-core` or `pyansys` depending on your implementation)*.

**Installation:**
```bash
git clone [https://github.com/bhahohocrzmoto/FEMTB880Case-0.git](https://github.com/bhahohocrzmoto/FEMTB880Case-0.git)
cd FEMTB880Case-0
pip install -r requirements.txt
