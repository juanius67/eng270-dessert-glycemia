# ENG-270 — Dessert Glycemia Simulator

Simulates post-prandial glucose/insulin responses for South American desserts using a dual-pool appearance model driven by macronutrient YAML files. The Bergman minimal model is solved with either a C shared library or a pure-Python RK4 fallback. Outputs include PNG figures in `figures/` and CSV summaries in `tables/`.

---

## Repository layout

```text
C/                    # C backend (model.c / model.h)
configs/
  ├─ params.yaml      # physiological parameters + nutrition→kinetics knobs
  ├─ desserts.yaml    # editable dessert catalog (macros per portion)
  └─ frozen/          # frozen YAML inputs used by --reproduce
run.py                # single CLI entry point
requirements.txt      # Python dependencies
LICENSE, CITATION.cff
tests/                # pytest sanity checks
```

Generated artifacts land in `build/`, `figures/`, and `tables/` (all ignored by git).

---

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate  # or .\.venv\Scripts\activate on Windows
pip install -r requirements.txt
python run.py              # uses Python backend if C build unavailable
```

The first run automatically attempts to compile the C solver; if no compiler is found the Python integrator is used instead. Results appear in:

* `figures/zoom/` and `figures/full/` – glucose, insulin, and D(t) plots;
* `tables/summary.csv` – dessert-level metrics (peaks, AUCs, kinetics);
* `build/env.json`, `build/build.json`, `build/sanity.json` – environment and sanity diagnostics.

---

## Grader quick start

```bash
python run.py --reproduce
```

This command wipes prior outputs, rebuilds the C backend, and re-runs every frozen YAML under `configs/frozen/`. It emits the exact figures and tables expected by the ENG-270 grader:

* `figures/zoom/`, `figures/full/`
* `tables/summary.csv`
* `build/env.json`, `build/build.json`, `build/sanity.json`

---

## Configuration

* `configs/params.yaml` stores the minimal-model parameters with units (`S_G_min1`, `p2_min1`, `p3_min1_per_uU_per_mL`, `phi_G_uU_mL_min1_per_mg_dL`, `G_thr_mg_dL`, `n_min1`, `Gb_mg_dL`, `Ib_uU_mL`) plus integration settings and nutrition modifiers. Inline comments document units.
* `configs/desserts.yaml` maps dessert names to macronutrients per portion (`carbs_g`, `sugars_g`, `fiber_g`, `fat_g`, `protein_g`, `portion_g`, optional `barcode`).
* `configs/frozen/*.yaml` are immutable per-dessert snapshots (one file per dessert) used exclusively by `--reproduce` for grader-proof reruns.

The nutrition pipeline converts each dessert YAML into a dual-exponential appearance profile `(A_fast, k_fast, A_slow, k_slow)` with a protein-driven insulin term `(A_prot, k_prot)`. Fat, fiber, and protein adjust appearance fractions and kinetics per the tunables in `params.yaml`.

---

## CLI reference

```
python run.py [--calibrate] [--sensitivity] [--window MIN] [--dpi DPI]
              [--no-plots] [--summary-only] [--reproduce]
```

* `--reproduce` – clean outputs, rebuild C, and run every frozen YAML.
* `--calibrate` – rescale appearance amplitudes so glucose peaks reach ~50 mg/dL.
* `--sensitivity` – perform ±20% fat/fiber/protein sweeps and write `tables/sensitivity.csv`.
* `--window` – override the zoom plot window in minutes (default 240).
* `--dpi` – plot resolution (default 150).
* `--no-plots` – skip PNG generation.
* `--summary-only` – emit CSV tables only (implies `--no-plots`).

All runs load physiology from `configs/params.yaml` and dessert macros from YAML; no hard-coded `(A, k)` paths remain.

---

## Testing

```bash
pytest
```

`tests/test_sanity.py` verifies:

* Steady-state stability when `D(t) = 0`.
* Solver invariance when halving the timestep (`<1%` change in key metrics).
* Parity between C and Python integrators when the shared library is available.

---

## Reproducibility notes

* `build/env.json` captures OS, Python, NumPy, Matplotlib, and git metadata.
* `build/build.json` records compiler and source hashes when the C backend is rebuilt.
* `build/sanity.json` stores dt-halving and end-state checks for the latest run.
* Frozen YAML files in `configs/frozen/` provide the exact nutrition inputs used for grading.
