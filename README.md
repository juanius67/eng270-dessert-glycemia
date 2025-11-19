# ENG-270 — Dessert Glycemia Simulator

Nutrition-label–driven simulator built for ENG-270. A single Python Command Line Interface (CLI) maps dessert macronutrient labels to a dual-exponential gut appearance model and calls a C RK4 implementation of the Bergman minimal model via `ctypes`.

The goal is to quantify how different dessert compositions (carbs, sugars, fibre, fat, protein) affect postprandial glucose excursions in a nominal non-diabetic adult, under a reproducible and fully scripted workflow.

---

## Project Description

From a bird's eye view, the project does three distinct things in the following order :

1. **Maps dessert labels to inputs of the Bergman minimal model**  
   Per-portion macronutrients from YAML files (namely: `carbs_g`, `sugars_g`, `fiber_g`, `fat_g`, `protein_g`, `portion_g`) are turned into:
- A two-pool (fast vs slow absorbing) gut appearance term :  
```math
D(t) = A_{\text{fast}} e^{-k_{\text{fast}} t} + A_{\text{slow}} e^{-k_{\text{slow}} t}
```
- An exponential insulin stimulus function ( $u(t)$ ) driven by protein.

2. **Simulates glucose–insulin dynamics in C**  
   This project's C backend implements the classical Bergman minimal model with RK4 integration (described in proposal and final report), parameterised for a nominal adult based on varied scientific literature (see references).

3. **Summarises and visualises outcomes**  
   Lastly, the Python portion performs simulations, produces figures (two types: zoom and full-day horizons), and writes CSV tables of peak glycemia, timing, and integrated exposure (iAUC) across all frozen dessert configs.

### Input files (“data”)

ALL user-supplied inputs live under `configs/`; this section describes are what they mean and how the project uses them :

- `configs/params.yaml`  
  Global physiological and appearance parameters with unit-encoded (i.e. the names of the parameters explicitly include their units of measurement) keys, these were taken from the literature :

  `Gb_mg_dL, Ib_uU_mL, S_G_min1, p2_min1, p3_min1_per_uU_mL, n_min1, Vd_dL, f_hep, f_app0, beta_fiber_per10g, beta_fat_per10g, k_fast0_min1, k_slow0_min1, alpha_prot_uU_mL_per_g, k_prot_min1, dt_min, t_end_min`.

- `configs/frozen/*.yaml`  
  One YAML file **per dessert**, each storing **per-portion** macros only. Below is an example for how a simple donut would be encoded :

  ```yaml
  name: sample_donut
  carbs_g: 45.0
  sugars_g: 25.0
  fiber_g: 2.0
  fat_g: 12.0
  protein_g: 4.0
  portion_g: 85.0
``

These YAML files are treated as **frozen inputs for grading and testing**. The “client” can add new desserts by copying an existing file, adjusting the macros and name, and re-running the pipeline.

**Everything** is derived from labels and `params.yaml`.

### Output files (“results”)

Running the pipeline fills the set of output directories below (all ignored by git so as to not upload junk) :

* `figures/zoom/`
  All figures with a time window of 0–240 min : overlays of glucose, insulin, and gut appearance for each dessert.

* `figures/full/`
  All figures with full-day (0–1440 min) overlays, used for long-time-horizon comparisons.

* `tables/summary.csv`
  CSV tables with one row per dessert with key metrics, for example :

  * `peak_delta_G` (mg/dL above baseline)
  * `t_peak_min` (time of peak glucose)
  * `iAUC_0_240` (mg·dL⁻¹·min, 0–240 min)
  * return-to-baseline indicators

* `tables/sensitivity.csv`
  Written by `--sensitivity` (this command will be described later) : A table with ±20 % perturbations of fat, fibre, and protein per dessert with resulting peaks/timings/iAUCs. This is essentially a *stress test* of the model, made to test if the simulator’s predictions remain meaningful even when inputs vary, highlighting which nutrients most strongly impact glycemia.


* `build/env.json`
  Here lies the environment metadata : Python version, NumPy version, OS, etc. Moreso used for bookkeeping than anything else; having this ensures that anyone reproducing results knows exactly what environment was used.

* `build/build.json`
  This the C backend build metadata : it is a compiler command and file hash of `src/model.c` for traceability.

The artefacts above are the **only** artefacts used in the final report; every figure or table in the report is sourced from these directories.

### Report

The final report is stored as:

* `report/eng270-dessert-glycemia-report.pdf`  

The report uses:

* Zoom and full-horizon plots from `figures/zoom/` and `figures/full/`,
* Summary statistics from `tables/summary.csv`,
* Sensitivity metrics from `tables/sensitivity.csv` (optional section on nutrient perturbations).

To regenerate all artefacts used in the report on a fresh machine:

```bash
pip install -r requirements.txt
python run.py --reproduce                 # regenerate baseline figures/tables
python run.py --reproduce --calibrate     # optional equal-amplitude comparison set
python run.py --reproduce --sensitivity   # refresh sensitivity table
python -m pytest -q                       # verify minimal tests
```

## Running the program

### Dependencies

* **Python:** 3.11+ (So far, tested and works on 3.11 and 3.12; on Linux/VDI use `python3` / `pip3` if `python` still points to a variation of 2.x)

* **Python packages:** pinned in `requirements.txt` (as aforementionned, install with `pip install -r requirements.txt`).
  Main libraries:

  * `numpy` for numerical work
  * `matplotlib` for plotting
  * `pyyaml` for configuration
  * `pytest` for tests

* **C compiler :** `gcc` or `clang` available on `PATH`. So as to minimise commands, the C backend is automatically compiled on first run, meaning no manual `make` step is required to make the program run.

#### Windows notes (MSYS2 + gcc)

On Windows, the recommended route to running the pipeline is MSYS2, to do so I reccommend the same steps I took :

1. Install MSYS2 from the official website.

2. Open **“MSYS2 MSYS”** and run:

   ```bash
   pacman -S --needed base-devel mingw-w64-x86_64-toolchain
   ```

3. Add `C:\msys64\mingw64\bin` to your **User PATH** (System Properties → Environment Variables).

4. Just in case, o a new terminal and verify:

   ```bash
   gcc --version
   ```

After this, `python run.py --reproduce` will be able to compile `src/model.c`.
Although not personally tested, it seems that a Chocolatey-managed MinGW (e.g., `choco install mingw`) can alternitavely be used, as long as the compiler is on `PATH`.

### Build

No explicit manual build step is needed:

* `run.py` loads `src/bindings.py`, which in turn:

  * Detects whether the shared library for `model.c` exists under `build/`,
  * Invokes the system C compiler with the appropriate flags if not,
  * Records the compiler command and file hash in `build/build.json`.

This keeps the build process **transparent and reproducible** while staying within a single CLI entry point.

### Execute

Core commands for the client’s engineering team:

```bash
# Install pinned dependencies
pip install -r requirements.txt

# Full baseline pipeline (all frozen desserts)
python run.py --reproduce

# Optional: equal-amplitude calibration set
python run.py --reproduce --calibrate

# Optional: ±20% perturbations of fat/fiber/protein (writes tables/sensitivity.csv)
python run.py --reproduce --sensitivity

# Quick headless smoke test (no figures, just tables + metadata)
python run.py --reproduce --summary-only --no-plots

# Minimal numerical tests (steady state + dt-halving)
python -m pytest -q
```

Flags:

* `--reproduce`
  Cleans prior outputs, rebuilds the C core if necessary, and replays every YAML in `configs/frozen/`. Produces all figures, summary tables, and build metadata.

* `--calibrate`
  For each dessert, recomputes appearance parameters so that **fast and slow pools have equal amplitudes at (t=0)** while:

  * preserving each pool’s decay rate, and
  * preserving the **total appearance dose**.

  Mathematically, amplitudes are set to
  (A_\text{eq} = \text{dose}*\text{tot} / (1/k*\text{fast} + 1/k_\text{slow})).
  This provides a controlled “what if both pools start equally strong?” comparison.

* `--sensitivity`
  Applies factors 0.8 and 1.2 to `fat_g`, `fiber_g`, and `protein_g` for each dessert, then records the resulting peaks, timings, and iAUCs in `tables/sensitivity.csv`.

* `--window {0-120,0-240,0-1440}`
  Zoom window for the “zoom” plots.

* `--dpi DPI`
  Plot resolution (default 150).

* `--no-plots`, `--summary-only`
  Speed up runs by skipping figures and/or only writing CSV tables.

* `--outdir PATH`
  Optional override for the base output directory (default: repository root).

---

## Model overview (for the engineering team)

### Bergman minimal model

The core is the classical Bergman minimal model of glucose regulation, with state variables:

* (G(t)): plasma glucose (mg/dL)
* (X(t)): “remote” insulin effect
* (I(t)): plasma insulin (µU/mL)

and dynamics:

[
\begin{aligned}
\frac{dG}{dt} &= -(S_G + X),(G - G_b) + D(t),\
\frac{dX}{dt} &= -p_2,X + p_3,(I - I_b),\
\frac{dI}{dt} &= -n,(I - I_b) + u(t),
\end{aligned}
]

where (D(t)) is gut appearance and (u(t)) is insulin input.
The C backend:

* Stores parameters internally after a `set_params` call,
* Integrates with a fixed-step classical RK4,
* Exposes a `step` function bridged to Python via `ctypes`.

### Mapping from labels to appearance and insulin

In `run.py`:

1. **Available carbohydrates**

   (C_\text{avail} = \max(0, \text{carbs}_g - 0.5 \cdot \text{fiber}_g))

   (simple fibre discount on total carbs).

2. **Fast fraction and appearance fraction**

   * Fast fraction (f_\text{fast} = \text{sugars}_g / \max(\text{carbs}_g, 10^{-9})).
   * Base appearance fraction `f_app0` reduced by fibre via `beta_fiber_per10g`.

3. **Rate modifiers from fat and fibre**

   A multiplicative factor:

   [
   k_\text{mod} = \frac{1}{1 + \beta_\text{fat} \cdot \frac{\text{fat}*g}{10} + \beta*\text{fiber} \cdot \frac{\text{fiber}_g}{10}},
   ]

   scales both `k_fast0_min1` and `k_slow0_min1`.

4. **Dose and amplitudes**

   * Total post-hepatic appearance dose (mg/dL·min equivalent):

     [
     \text{dose}*\text{tot}
     = \left(\frac{1000,C*\text{avail}}{V_d}\right),f_\text{app},(1 - f_\text{hep}).
     ]

   * Dose split into fast/slow pools by (f_\text{fast}), then:

     [
     A_\text{fast} = k_\text{fast} ,\text{dose}*\text{fast},\quad
     A*\text{slow} = k_\text{slow} ,\text{dose}_\text{slow}.
     ]

5. **Protein-driven insulin pulse**

   * Amplitude (A_\text{prot} = \alpha_\text{prot} \cdot \text{protein}_g).
   * Decay (k_\text{prot}) from `params.yaml`.
   * Stimulus (u(t) = A_\text{prot} e^{-k_\text{prot} t}).

These choices are motivated by literature on gastric emptying, intestinal transport, incretin response, and protein co-ingestion.

---

## Repository layout

```text
configs/
  ├─ params.yaml      # Physiology + appearance knobs (unit-encoded)
  └─ frozen/          # Immutable dessert YAMLs consumed by --reproduce
run.py                # Single CLI entry point
src/
  ├─ model.c,h        # RK4 implementation of Bergman minimal model
  └─ bindings.py      # ctypes loader + auto-build + param bridge
tests/
  └─ test_sanity.py   # Steady-state + dt-halving checks
requirements.txt      # numpy, matplotlib, pyyaml, pytest (pinned)
README.md             # This file
LICENSE               # Open-source license
CITATION.cff          # How to cite this repository
```

Git ignores `build/`, `figures/`, `tables/`, compiled libraries, and cache directories so that the repository stays lean while all results remain regenerable.

---

## Testing and reproducibility

The minimal test suite is designed to reassure the client’s engineering team that the numerical implementation behaves sensibly:

* `tests/test_sanity.py`:

  * **Fasting steady state:** no-meal simulation remains within 0.5 mg/dL (glucose) and 0.5 µU/mL (insulin) of baseline after 120 min.
  * **Time-step convergence:** peak ΔG under a standard dessert changes by less than 1 % when `dt` is halved.

Run all tests with:

```bash
python -m pytest -q
```

For a full end-to-end regeneration on a fresh machine, the recommended sequence is:

```bash
pip install -r requirements.txt
python run.py --reproduce --summary-only --no-plots
python -m pytest -q
python run.py --reproduce              # full figures/tables
python run.py --reproduce --calibrate  # optional calibration set
python run.py --reproduce --sensitivity
```

---

## Contributors

* **Juan Lucas de Oliveira** — EPFL SIE student, project design, implementation, and analysis.

(If additional collaborators join, list their names and contributions here.)

---

## Acknowledgments

### Data sources

* Dessert macro data are taken from nutrition labels of commercially available products; values are encoded manually into `configs/frozen/*.yaml`.

### Code

* Repository structure and workflow were inspired by the official ENG-270 project template (`stakahama/sie-eng270-project-template`).
* Numerical and physiological modelling draws on the literature listed in the **References** section below.
* The project uses `numpy`, `matplotlib`, and `pyyaml`, which are cited implicitly via the Python ecosystem.

For citation of this repository itself, see `CITATION.cff`.

---

## References

The modelling choices and interpretation are grounded in the following references (from the original project proposal):

1. Marathe, C. S., Rayner, C. K., Jones, K. L., & Horowitz, M. (2013). Relationships between gastric emptying, postprandial glycemia, and incretin hormones. *Diabetes Care*, 36(5), 1396–1405.

2. Koepsell, H. (2020). Glucose transporters in the small intestine in health and disease. *Pflügers Archiv – European Journal of Physiology*, 472(9), 1207–1248.

3. Park, S. Y., Gautier, J.-F., & Chon, S. (2021). Assessment of insulin secretion and insulin resistance in humans. *Diabetes & Metabolism Journal*, 45(5), 641–654.

4. Bergman, R. N. (2021). Origins and history of the minimal model of glucose regulation. *Frontiers in Endocrinology*, 11, 583016.

5. Bergman, R. N., Ider, Y. Z., Bowden, C. R., & Cobelli, C. (1979). Quantitative estimation of insulin sensitivity. *American Journal of Physiology*, 236(6), E667–E677.

6. Bergman, M. (2021). The 1-hour plasma glucose: Common link across the glycemic spectrum. *Frontiers in Endocrinology*, 12, 1–9.

7. Cobelli, C., Dalla Man, C., Sparacino, G., Magni, L., De Nicolao, G., & Kovatchev, B. P. (2009). Diabetes: Models, signals, and control. *IEEE Reviews in Biomedical Engineering*, 2, 54–96.

8. Taylor, P. N., Collins, K. S., Lam, A., et al. (2023). C-peptide and metabolic outcomes in trials of disease-modifying therapy in new-onset type 1 diabetes: An individual participant meta-analysis. *The Lancet Diabetes & Endocrinology*, 11(12), 915–925.

9. Slavin, J. (2013). Fiber and prebiotics: Mechanisms and health benefits. *Nutrients*, 5(4), 1417–1435.

10. Camilleri, M. (2006). Integrated upper gastrointestinal response to food intake. *Gastroenterology*, 131(2), 640–658.

11. Frid, A. H., Nilsson, M., Holst, J. J., & Björck, I. M. E. (2005). Effect of whey on blood glucose and insulin responses to composite meals in type 2 diabetic subjects. *American Journal of Clinical Nutrition*, 82(1), 69–75.

```
