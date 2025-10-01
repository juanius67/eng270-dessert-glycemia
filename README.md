
# ENG-270 — Dessert Glycemia Simulator (Python + C)

Simulates post-prandial glucose/insulin responses to four South-American desserts (chocotorta, brigadeiro, alfajor, açaí) using a Bergman minimal-model ODE with a dessert-specific glucose appearance $D(t)=A e^{-kt}$.
Outputs: **PNG figures** in `figures/` and a **summary CSV** in `tables/`.

---

## Repo structure

```
C/               # C RK4 solver (exports simulate)
  ├─ model.c
  ├─ model.h
  └─ Makefile
py/              # Python orchestration
  ├─ run_all.py  # runs all desserts → figures/ + tables/
  ├─ interface.py  # ctypes bridge into C solver
  ├─ desserts.py   # A,k per dessert
  ├─ metrics.py    # peaks, AUCs
  └─ plotting.py   # saves PNGs into ../figures
configs/
  └─ params.yaml   # physiological params + horizon
figures/          # generated (not tracked)
tables/           # generated (not tracked)
requirements.txt  # numpy, matplotlib, pyyaml
README.md
```

---

## Requirements

* Python 3.10+
* `numpy`, `matplotlib`, `pyyaml`, `scipy`, `pytest` (installed below)
* Optional C build: **gcc** (MSYS2 MinGW64 on Windows) or any C11 compiler

---

## Quick start (Codex environment)

**Container image:** `universal`
**Setup script (Manual):**

```bash
pip install -r requirements.txt
make -C C
python -m py.run_all
```

Results: `figures/*.png`, `figures/*_overlay.png`, `tables/summary.csv`.

---

## Quick start (local, VS Code on Windows)

### 1) Python env & deps

```powershell
cd C:\Users\<you>\path\to\eng270-dessert-glycemia
py -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -U pip
pip install -r requirements.txt
```

### 2) Build C (fast path) — pick ONE:

**A. Using MSYS2 MinGW64 gcc (recommended):**

```powershell
# one-time: add to PATH for this session
$env:Path += ";C:\msys64\mingw64\bin"
gcc -shared -O2 -std=c11 -o C\model.dll C\model.c
```

**B. Using Make (if available):**

```powershell
make -C C
```

(Produces `C\model.dll` on Windows; `C/libmodel.so` on Linux/macOS.)

### 3) Run

```powershell
py -m py.run_all
```

Outputs appear in `figures\` and `tables\`.

> No C compiler? Install gcc/clang (see above) or reuse a prebuilt `C/libmodel.so` / `C/model.dll`; the pipeline requires the shared library.

---

## Reproduce everything (single command)

From repo root (venv active):

```bash
python -m py.run_all
```

---

## Parameters & scenarios

### Model/solver settings

`configs/params.yaml`:

```yaml
p1: 0.03
p2: 0.012
p3: 4.92e-6
p4: 0.0039
p5: 80
p6: 0.265
Gb: 80
Ib: 7
dt: 0.5
t_end: 1440
```

### Dessert inputs (edit to tune)

`py/desserts.py`:

```python
desserts = {
    "chocotorta": (250.0, 0.08),
    "brigadeiro": (250.0, 0.50),
    "alfajor":    (200.0, 0.10),
    "acai":       (300.0, 0.15),
}
# (A,k) control dose and speed; reduce A to lower peaks; reduce k to slow absorption.
```

---

## Outputs

* `figures/<dessert>_glucose.png`
* `figures/<dessert>_insulin.png`
* `figures/glucose_overlay.png`
* `figures/insulin_overlay.png`
* `tables/summary.csv` with dessert-level metrics (peaks, peak times, AUCs)
  feeding the Results/Discussion tables.

---

## Testing

Run the fast sanity suite (after building the shared library):

```bash
pytest
```

---

## Commands

* Run tests: `python -m pytest -q`
* Reproduce figures/tables:
  * Dose-driven: `python run.py --all`
  * Equalized peaks: `python run.py --all --calibrate`
  * Include LaTeX export: `python run.py --all --latex`
  * 24 h run (dose-driven, nutrition-aware): `python run.py --all`
  * Equalized peaks (appendix): `python run.py --all --calibrate`
  * Disable nutrition mapping (legacy): `python run.py --all --no-nutrition`
  * Horizon/step are set in `configs/params.yaml` (`t_end=1440`, `dt=0.5`).

---

## How it works (brief)

* C exports `simulate(G,I,nsteps,dt,Params*)` (fixed-step RK4).
* Python loads `C/model.dll` or `C/libmodel.so` via `ctypes` and streams results back to NumPy.
* `run_all.py` iterates desserts → runs sim → saves plots & metrics.

---

## Troubleshooting

**Figures/CSV not appearing**
Ensure `py/plotting.py` and `py/metrics.py` use repo-relative paths (already fixed). Run from repo root:

```bash
python -m py.run_all
```

**`make: command not found` (Windows)**
Use direct gcc build:

```powershell
$env:Path += ";C:\msys64\mingw64\bin"
gcc -shared -O2 -std=c11 -o C\model.dll C\model.c
```

**`pip` blocked (PEP 668 / externally managed)**
Use a venv:

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows: .\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

**`yaml` module missing**
`pip install pyyaml` in the active venv.

**Using PowerShell: execution policy**

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
```

---

## Clean build

```bash
# remove generated outputs
git clean -dfX figures tables
# rebuild C (if using make)
make -C C clean && make -C C
```

---

## Citation & attribution

* Bergman minimal model; dessert $D(t)$ modeled as exponential appearance.
* Any external or AI-assisted code should be cited in the report/README section “Attribution”.

---

## License

Academic/educational use for ENG-270 coursework.

---

## How to reproduce

- Install: `pip install -r requirements.txt`
- Build (optional, Windows x64 example): `cl /nologo /LD C\model.c /Fe:C\model.dll`
- Run (dose-driven): `python run.py --all`
- Run (equalized peaks): `python run.py --all --calibrate`
- Sanity checks: `python run.py --sanity`
- Outputs: figures/*.png and tables/summary.csv (regenerated each run).

Graders only need `run.py` + `configs/params.yaml` + `C/` to reproduce.

## Data sources
- Havanna Alfajor Super dulce de leche (EAN 7791875101598), per-100g label:
  ~393 kcal; fat 17 g (sat 7 g); carbs 52 g (sugars 38 g); protein 7 g; sodium ~107 mg.
  Source: OpenFoodFacts product page and retail listing with nutrition per 100 g.
- Per-unit cross-check (70 g): ~275 kcal; fat ~12 g; carbs ~37 g; protein ~5 g.
  Source: FatSecret (Havanna Alfajor Súper Dulce de Leche, 70 g).

## Usage

- 24h run (nutrition-aware): `python run.py --all`
- Equalized peaks (appendix): `python run.py --all --calibrate`
- Import from OpenFoodFacts by barcode: `python run.py --import-off alfajor 7791875101598 70`
- Re-run pipeline after import: `python run.py --all`
- Notes: fat slows gastric emptying; fiber reduces appearance; protein drives insulin pulse; sugars fraction → faster absorption. All figures/tables regenerate from one command.
