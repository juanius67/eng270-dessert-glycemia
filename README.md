# ENG-270 — Dessert Glycemia Simulator

Nutrition-label–driven simulator built for ENG-270. A single Python Command Line Interface (CLI) maps dessert macronutrient labels to a dual-exponential gut appearance model and calls a C RK4 implementation of the Bergman minimal model via `ctypes`.

The goal is to quantify how different dessert compositions (carbs, sugars, fiber, fat, protein) affect postprandial glucose excursions in a nominal, non-diabetic adult, under a reproducible and fully scripted workflow.

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


These YAML files are treated as **frozen inputs for testing**. The “client” can add new desserts by copying an existing file, adjusting the macros and name, and re-running the pipeline.

**Everything** is derived from labels and `params.yaml`.

### Output files (“results”)

Running the pipeline creates and fills the set of output directories below :

* `figures/zoom/`
  All figures with a time window of 0–240 min : overlays of glucose, insulin, and gut appearance for each dessert.

* `figures/full/`
  All figures with full-day (0–1440 min) overlays, used for long-time-horizon comparisons.
  
* `figures/summary_dashboard.png`
A consolidated bar-chart strip comparing Peak G, iAUC, and Timings across all desserts side-by-side.

* `tables/summary.csv`
  CSV tables with one row per dessert with key metrics, for example :

  * `peak_delta_G` (mg/dL above baseline)
  * `t_peak_min` (time of peak glucose)
  * `iAUC_0_240` (mg·dL⁻¹·min, 0–240 min)
  * return-to-baseline indicators

* `tables/sensitivity.csv`
  Written by `--sensitivity` (this command will be described later) : A table with ±20 % perturbations of fat, fiber, and protein per dessert with resulting peaks/timings/iAUCs. This is essentially a *stress test* of the model, made to test if the simulator’s predictions remain meaningful even when inputs vary, highlighting which nutrients most strongly impact glycemia.


* `build/env.json` 
  Here lies the environment metadata : Python version, NumPy version, OS, etc. Moreso used for bookkeeping than anything else; having this ensures that anyone reproducing results knows exactly what environment was used.

* `build/build.json`
  This the C backend build metadata : it is a compiler command and file hash of `src/model.c` for traceability.

Except for the two json files created through AI-generated commands for bookkeeping and cross-referencing across devices, the artefacts above are the **only** artefacts used in the final report; every figure or table in the report is sourced from these directories.

### Report

The final report is stored as:

* `report/Juan_Lucas_de_Oliveira-ENG-270_Project_Report.pdf`  

The report uses:

* Zoom and full-horizon plots from `figures/zoom/` and `figures/full/`,
* Summary statistics from `tables/summary.csv`

To regenerate all artefacts used in the report on a fresh machine:

```bash
pip install -r requirements.txt

# Single full pipeline: baseline + sensitivity (matches report)
python run.py --reproduce --sensitivity

# Optional: regenerate equal-amplitude calibration figures
# (WARNING: overwrites figures/summary for the chosen outdir)
python run.py --reproduce --calibrate

# Numerical sanity tests (same C RK4 core)
python -m pytest -q

```

## Running the program

### Dependencies

* **Python:** 3.11+ (So far, tested and works on 3.11 and 3.12; on Linux/VDI use `python3` / `pip3` and `python3 -m pip` for all commands if `python` still points to a variation of 2."x")

* **Python packages:** pinned in `requirements.txt` (as aforementionned, install with `pip install -r requirements.txt`).
  Main libraries:

  * `numpy` for numerical work
  * `matplotlib` for plotting
  * `pyyaml` for configuration
  * `pytest` for tests

* **C compiler :** `gcc` or `clang` available on `PATH`. So as to minimise commands, the C backend is automatically compiled on first run, meaning no manual `make` step is required to make the program run.

#### Windows notes (MSYS2 + gcc)

On Windows, the recommended route to running the pipeline is MSYS2, to do so I reccommend the same steps I took on my personal machine (ARM CPU Windows Laptop) :

1. Install MSYS2 from the official website.

2. Open **“MSYS2 MSYS”** and run:

   ```bash
   pacman -S --needed base-devel mingw-w64-x86_64-toolchain
   ```

3. Add `C:\msys64\mingw64\bin` to your **User PATH** (System Properties → Environment Variables).

4. Just in case, open a new terminal and verify:

   ```bash
   gcc --version
   ```

After this, `python run.py --reproduce` should be able to compile `src/model.c`.

### Build

No explicit manual build step is needed:

* `run.py` loads `src/bindings.py`, which in turn:

  * Detects whether the shared (dynamic) library (built by the pipeline) for `model.c` exists under `build/`,
  * Invokes the system C compiler with the appropriate flags if not,
  * Records the compiler command and file hash in `build/build.json` (this is done for documentation and error analysis reasons).

All of this happens behind the run commands, so the client only ever needs to run the CLI

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

# Quick headless smoke test (no figures, just tables + metadata. A quick command to see if everything is up and running)
python run.py --reproduce --summary-only --no-plots

# Minimal numerical tests (steady state + dt-halving)
python -m pytest -q
```

Flags:
Below is a list of flags a user can use; these small modifiers change behaviour without requiring additional arguments:

* `--reproduce`
  Cleans prior outputs, rebuilds the C core (only if necessary), and reruns for every YAML in `configs/frozen/`. Produces all figures, summary tables, and build metadata.

* `--calibrate`
  For each dessert, recomputes appearance parameters so that **fast and slow pools have equal amplitudes at (t=0)** whilst preserving each pool’s decay rate all the while preserving the **total appearance dose**.

  Mathematically, amplitudes are set to:
  ```math
  (A_\text{eq} = \text{dose}_\text{tot} / (1/k_\text{fast} + 1/k_\text{slow})).
  ```
  This provides a controlled “what if both pools start equally strong?” comparison.

* `--sensitivity`
  Applies factors 0.8 and 1.2 (±20%) to `fat_g`, `fiber_g`, and `protein_g` for each dessert, then records the resulting peaks, timings, and iAUCs in `tables/sensitivity.csv`. 

* `--window {0-120,0-240,0-1440}`
  Zoomed-in window for the “zoom” plots.

* `--dpi DPI`
  Plot resolution. The default is 150.

* `--no-plots`, `--summary-only`
  Speed up runs by skipping figures and/or only writing CSV tables.

* `--outdir PATH`
  Optional override for the base output directory (default: repository root).

---

## Model overview 

### Bergman minimal model

The core is the classical Bergman minimal model of glucose regulation, with state variables:

* (G(t)): plasma glucose (mg/dL)
* (X(t)): "remote" insulin effect
* (I(t)): plasma insulin (µU/mL)

and dynamics:
```math

\begin{aligned}
\frac{dG}{dt} &= -(S_G + X),(G - G_b) + D(t),\
\frac{dX}{dt} &= -p_2,X + p_3,(I - I_b),\
\frac{dI}{dt} &= -n,(I - I_b) + u(t),
\end{aligned}

```
where (D(t)) is gut appearance and (u(t)) is insulin input.
The C backend:

* Stores parameters internally after a `set_params` call,
* Integrates with a, fixed-step, classical RK4,
* Makes a `step` (dt) function visible and bridged to Python via `ctypes`.

### Mapping from labels to appearance and insulin

In `run.py`:

1. **Available carbohydrates**

   (C_\text{avail} = \max(0, \text{carbs}_g - 0.5 \cdot \text{fiber}_g))

   (simple fiber discount on total carbs).

2. **Fast fraction and appearance fraction**

   * Fast fraction:
```math
 (f_\text{fast} = \text{sugars}_g / \max(\text{carbs}_g, 10^{-9})).
```
   * Base appearance fraction `f_app0` reduced by fiber via `beta_fiber_per10g`.

3. **Rate modifiers from fat and fiber**

   A multiplicative factor:

```math
   k_\text{mod} = \frac{1}{1 + \beta_\text{fat} \cdot \frac{\text{fat}*g}{10} + \beta_\text{fiber} \cdot \frac{\text{fiber}_g}{10}},
```

   scales both `k_fast0_min1` and `k_slow0_min1`.

4. **Dose and amplitudes**

   * Total post-hepatic (post liver) appearance dose (mg/dL·min equivalent):

```math
     \text{dose}_\text{tot}
     = \left(\frac{1000,C_\text{avail}}{V_d}\right),f_\text{app},(1 - f_\text{hep}).
```
   * Dose split into fast/slow pools by (f_\text{fast}), then:

```math
     A_\text{fast} = k_\text{fast} ,\text{dose}_\text{fast},\quad
     A_\text{slow} = k_\text{slow} ,\text{dose}_\text{slow}.
```

5. **Protein-driven insulin pulse**

   * Amplitude (A_\text{prot} = \alpha_\text{prot} \cdot \text{protein}_g).
   * Decay (k_\text{prot}) from `params.yaml`.
   * Stimulus (u(t) = A_\text{prot} e^{-k_\text{prot} t}).

These choices are motivated by literature on gastric emptying, intestinal transport, incretin response, and protein co-ingestion (see report)

---


Git ignores `build/`, `figures/`, `tables/`, compiled libraries, and cache directories so that the repository stays lean while all results remain regenerable.

---

## Testing and reproducibility

The minimal test suite is designed to reassure the client’s engineering team that the numerical implementation behaves sensibly:

`tests/test_sanity.py`:

* **Fasting steady-state test:**  
  Runs a "no-meal" simulation and checks that glucose and insulin remain close to their baseline values after 120 min (within ±0.5 mg/dL for glucose and ±0.5 µU/mL for insulin); this reflects what should be the case in real life.

* **Time-step convergence test:**  
  Simulates a standard dessert with time-step `dt`, then repeats the simulation with `dt/2`.  
  The test passes if the peak glucose rise (ΔG) differs by less than 1 % between the two runs.


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

* **Juan Lucas de Oliveira**

### AI Agents and thier uses:
* **ChatGPT / ChatGPT Codex 5.0 & 5.1  (OpenAI)** — Used as a coding assistant for refactoring, infrastructure helpers (metadata, JSON, error handling), and documentation suggestions; all AI-generated code was reviewed, adapted over time, and cited.
* **“Jules” (Google AI) (Gemini 3.0)** — Used as an auxiliary assistant for brainstorming implementation options but mostly for intensive documentation; any outputs were treated as suggestions and integrated only after **human** review.


---

## Acknowledgments

### Data sources

* Dessert macro data are taken from nutrition labels of commercially available products on OpenFoodFacts.org; values are encoded manually into `configs/frozen/*.yaml`.

### Code

* Numerical and physiological modelling draws on the literature listed in the **References** section below.
* The project uses `numpy`, `matplotlib`, and `pyyaml`.


---

## References

The modelling choices, constant values and interpretation are grounded in the following references (some from the original project proposal):

[1] Michael Bergman. The 1-Hour Plasma Glucose: Common Link Across the Glycemic Spectrum. Frontiers in Endocrinology, 12, September 2021. ISSN 1664-2392. doi: 10.3389/fendo.2021.752329. URL https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2021.752329/full. Publisher: Frontiers.

[2] Richard N. Bergman. Origins and History of the Minimal Model of Glucose Regulation. Frontiers in Endocrinology, 11:583016, February 2021. ISSN 1664-2392. doi: 10.3389/fendo.2020.583016. URL https://pmc.ncbi.nlm.nih.gov/articles/PMC7917251/.

[3] R. N. Bergman, Y. Z. Ider, C. R. Bowden, and C. Cobelli. Quantitative estimation of insulin sensitivity. The American Journal of Physiology, 236(6):E667–677, June 1979. ISSN 0002-9513. doi: 10.1152/ajpendo.1979.236.6.E667.

[4] Michael Camilleri. Integrated Upper Gastrointestinal Response to Food Intake. Gastroenterology, 131(2):640–658, August 2006. ISSN 0016-5085, 1528-0012. doi: 10.1053/j.gastro.2006.03.023. URL https://www.gastrojournal.org/article/S0016-5085(06)00572-5/fulltext. Publisher: Elsevier.

[5] Claudio Cobelli, Chiara Dalla Man, Giovanni Sparacino, Lalo Magni, Giuseppe De Nicolao, and Boris P. Kovatchev. Diabetes: Models, Signals, and Control. IEEE reviews in biomedical engineering, 2:54–96, January 2009. ISSN 1937-3333. doi: 10.1109/RBME.2009.2036073. URL https://pmc.ncbi.nlm.nih.gov/articles/PMC2951686/.

[6] Anders H. Frid, Mikael Nilsson, Jens Juul Holst, and Inger ME Björck. Effect of whey on blood glucose and insulin responses to composite breakfast and lunch meals in type 2 diabetic subjects. The American Journal of Clinical Nutrition, 82(1):69–75, July 2005. ISSN 0002-9165. doi: 10.1093/ajcn/82.1.69. URL https://www.sciencedirect.com/science/article/pii/S0002916523295118.

[7] Hermann Koepsell. Glucose transporters in the small intestine in health and disease. Pflugers Archiv, 472(9):1207–1248, 2020. ISSN 0031-6768. doi: 10.1007/s00424-020-02439-5. URL https://pmc.ncbi.nlm.nih.gov/articles/PMC7462918/.

[8] Chinmay S. Marathe, Christopher K. Rayner, Karen L. Jones, and Michael Horowitz. Relationships between gastric emptying, postprandial glycemia, and incretin hormones. Diabetes Care, 36(5):1396–1405, May 2013. ISSN 1935-5548. doi: 10.2337/dc12-1609.

[9] So Young Park, Jean-François Gautier, and Suk Chon. Assessment of Insulin Secretion and Insulin Resistance in Human. Diabetes & Metabolism Journal, 45(5):641–654, September 2021. ISSN 2233-6087. doi: 10.4093/dmj.2021.0220.

[10] Joanne Slavin. Fiber and Prebiotics: Mechanisms and Health Benefits. Nutrients, 5(4):1417–1435, April 2013. ISSN 2072-6643. doi: 10.3390/nu5041417. URL https://www.mdpi.com/2072-6643/5/4/1417. Publisher: Multidisciplinary Digital Publishing Institute.

[11] Peter N. Taylor, Kimberly S. Collins, Anna Lam, Stephen R. Karpen, Brianna Greeno, Frank Walker, Alejandro Lozano, Elnaz Atabakhsh, Simi T. Ahmed, Marjana Marinac, Esther Latres, Peter A. Senior, Mark Rigby, Peter A. Gottlieb, Colin M. Dayan, and Trial Outcome Markers Initiative collaboration. C-peptide and metabolic outcomes in trials of disease modifying therapy in new-onset type 1 diabetes: an individual participant meta-analysis. The Lancet. Diabetes & Endocrinology, 11(12):915–925, December 2023. ISSN 2213-8595. doi: 10.1016/S2213-8587(23)00267-X.

[12] Acai berry bowl. Open Food Facts. URL https://world.openfoodfacts.org/product/9120054751695/acai-berry-bowl. [Accessed: 2025-11-29].

[13] Alfajor super dulce de leche – Havanna – 70 g. Open Food Facts. May 2024. URL https://world.openfoodfacts.org/product/7791875101598/alfajor-super-dulce-de-leche-havanna. [Accessed: 2025-11-29].

[14] Brigadeiro. Open Food Facts. July 2024. URL https://world.openfoodfacts.org/product/7896434920662/brigadeiro. [Accessed: 2025-11-29].

[15] Chocotorta – Arcor – 53 g. Open Food Facts. September 2025. URL https://world.openfoodfacts.org/product/7790580142179/chocotorta-arcor. [Accessed: 2025-11-29].

