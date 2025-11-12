# ENG-270 — Dessert Glycemia Simulator

Nutrition labels in YAML drive a dual-exponential gut appearance model that feeds the Bergman minimal model. A C RK4 core (exposed via `ctypes`) updates glucose, insulin, and remote insulin effect for a nominal, non-diabetic adult. The command-line interface in `run.py` is the only entry point and produces figures, tables, and build metadata.

## Grader quick start

```bash
pip install -r requirements.txt
python run.py --reproduce
python run.py --calibrate --sensitivity
python -m pytest -q
```

`--reproduce` wipes prior artifacts, rebuilds the C core if necessary, and replays every frozen YAML under `configs/frozen/`. Outputs are written to:

* `figures/zoom/` (default 0–240 min overlays for glucose/insulin/appearance)
* `figures/full/` (0–1440 min day-long overlays)
* `tables/summary.csv` (peak ΔG, times, iAUCs, baseline return, insulin peak)
* `build/env.json` (Python, NumPy, OS)
* `build/build.json` (compiler command and SHA256 of `src/model.c`)

`--calibrate` equalises the fast and slow appearance amplitudes (keeping decay rates fixed). `--sensitivity` performs ±20 % sweeps of fat, fiber, and protein and writes `tables/sensitivity.csv`.

## Repository layout

```text
configs/
  ├─ params.yaml      # unit-encoded physiology + appearance knobs
  ├─ desserts.yaml    # editable dessert catalog (per-portion macros)
  └─ frozen/          # immutable YAMLs consumed by --reproduce
run.py                # single CLI
src/model.c,h         # RK4 implementation of the Bergman minimal model
src/bindings.py       # ctypes loader with auto-build logic
requirements.txt      # numpy, matplotlib, pyyaml, pytest
LICENSE, CITATION.cff
tests/test_sanity.py  # steady-state + dt-halving checks
```

Git ignores `build/`, `figures/`, `tables/`, compiled libraries, and cache directories.

## Parameters and inputs

`configs/params.yaml` stores the defaults used by the simulator. Keys encode units explicitly and match the specification:

* `Gb_mg_dL: 90`
* `Ib_uU_mL: 7`
* `S_G_min1: 0.025`
* `p2_min1: 0.025`
* `p3_min1_per_uU_mL: 1.3e-3`
* `n_min1: 0.14`
* `Vd_dL: 110`
* `f_hep: 0.25`
* `f_app0: 0.30`
* `beta_fiber_per10g: 0.05`
* `beta_fat_per10g: 0.06`
* `k_fast0_min1: 0.35`
* `k_slow0_min1: 0.07`
* `alpha_prot_uU_mL_per_g: 0.06`
* `k_prot_min1: 0.05`
* `dt_min: 0.5`
* `t_end_min: 1440`

Legacy `p1`–`p4` keys are mapped to the new names with a warning for backward compatibility.

Each dessert YAML provides per-portion macros only—no hard-coded `(A, k)` pairs remain. A minimal example:

```yaml
name: sample_donut
carbs_g: 45.0
sugars_g: 25.0
fiber_g: 2.0
fat_g: 12.0
protein_g: 4.0
portion_g: 85.0
```

During a run, carbohydrates determine the total appearance dose after hepatic extraction (`f_hep`). Sugars route to the fast pool, the remainder to the slow pool. Fat and fiber adjust the appearance fractions and decay rates; protein induces an exponential insulin stimulus.

## CLI usage

```
python run.py [--reproduce] [--calibrate] [--sensitivity]
              [--window {0-120,0-240,0-1440}] [--dpi DPI]
              [--no-plots] [--summary-only] [--outdir PATH]
```

* `--reproduce` – clean outputs and rerun frozen desserts only.
* `--calibrate` – set equal fast/slow amplitudes while keeping decay rates.
* `--sensitivity` – write `tables/sensitivity.csv` after ±20 % fat/fiber/protein sweeps.
* `--window` – adjust the zoom-plot window (default `0-240`).
* `--dpi` – change figure resolution (default `150`).
* `--no-plots` – skip figure creation.
* `--summary-only` – emit tables and build metadata only (implies no plots).
* `--outdir` – redirect outputs (defaults to the repository root).

The CLI always loads constants from `configs/params.yaml`, pushes them into the C layer via `set_params`, and iterates the RK4 integrator (`step(...)`) with the computed appearance and insulin stimuli.

## Testing

```bash
python -m pytest -q
```

`tests/test_sanity.py` checks that the fasting state remains within ±0.5 mg/dL (glucose) and ±0.5 μU/mL (insulin) for 120 min, and that halving the timestep changes the peak glucose excursion by less than 1 %.
