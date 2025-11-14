# ENG-270 — Dessert Glycemia Simulator
### Windows build notes (gcc)
To compile the C backend on Windows, install MSYS2 and add gcc to PATH:

1. Install MSYS2: https://www.msys2.org/
2. Open “MSYS2 MSYS” and run:
```

pacman -S --needed base-devel mingw-w64-x86_64-toolchain

```
3. Add `C:\msys64\mingw64\bin` to your **User PATH** (System Properties → Environment Variables).
4. Open a new terminal and verify:
```

gcc --version

```
Now `python run.py --reproduce` will be able to compile `src/model.c`.
```

(If you prefer Chocolatey/MinGW, a one-liner works too:
`choco install mingw` then add `C:\ProgramData\chocolatey\bin` to PATH.)
--------------------------------
Quick start / Grader walkthrough
--------------------------------
1. `pip install -r requirements.txt`
2. `python run.py --reproduce` – compiles `src/model.c` if needed and regenerates figures, tables, and build metadata from `configs/frozen/*.yaml`.
3. `python run.py --calibrate` – replays frozen desserts with equal fast/slow amplitudes while keeping decay rates fixed.
4. `python run.py --sensitivity` – writes `tables/sensitivity.csv` after ±20% fat/fiber/protein sweeps.
5. `python -m pytest -q` – validates steady state and dt-halving convergence with the shared RK4 integrator.

**Inputs.** The pipeline reads *only* frozen per-dessert YAML files in `configs/frozen/` at grading time. Each file stores per-portion: `carbs_g, sugars_g, fiber_g, fat_g, protein_g, portion_g`. These are mapped to a dual-pool appearance
(D(t) = A_\text{fast} e^{-k_\text{fast} t} + A_\text{slow} e^{-k_\text{slow} t}) with fat/fiber modifiers, and a protein-driven insulin pulse (u(t)).

**Parameters.** Model constants live in `configs/params.yaml` with unit-encoded keys:
`Gb_mg_dL, Ib_uU_mL, S_G_min1, p2_min1, p3_min1_per_uU_mL, n_min1, Vd_dL, f_hep, f_app0, beta_fiber_per10g, beta_fat_per10g, k_fast0_min1, k_slow0_min1, alpha_prot_uU_mL_per_g, k_prot_min1, dt_min, t_end_min`.

Nutrition labels in YAML drive a dual-exponential gut appearance model that feeds the Bergman minimal model. A C RK4 core (exposed via `ctypes`) updates glucose, insulin, and remote insulin effect for a nominal, non-diabetic adult. The command-line interface in `run.py` is the only entry point and produces figures, tables, and build metadata.

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

`configs/params.yaml` stores the defaults used by the simulator. Keys encode units explicitly and match the specification above.

Each dessert YAML provides per-portion macros only. A minimal example (stored in `configs/frozen/`):

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
