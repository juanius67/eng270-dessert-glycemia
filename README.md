
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
  ├─ interface.py  # ctypes bridge + Python RK4 fallback
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
* `numpy`, `matplotlib`, `pyyaml` (installed below)
* Optional C build: **gcc** (MSYS2 MinGW64 on Windows) or any C11 compiler

---

## Quick start (Codex environment)

**Container image:** `universal`
**Setup script (Manual):**

```bash
pip install -r requirements.txt
make -C C
```

Then in the terminal:

```bash
python py/run_all.py
```

Results: `figures/*.png`, `tables/summary.csv`.

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
py py\run_all.py
```

Outputs appear in `figures\` and `tables\`.

> No C compiler? The Python RK4 fallback runs automatically (slower but fine).

---

## Reproduce everything (single command)

From repo root (venv active):

```bash
python py/run_all.py
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
dt: 0.01
t_end: 180
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
* `tables/summary.csv` with:

  ```
  dessert,peak_G(t0 units),t_peak_G(min),AUC_G,peak_I,t_peak_I(min),AUC_I
  ...
  ```

---

## How it works (brief)

* C exports `simulate(G,I,nsteps,dt,Params*)` (fixed-step RK4).
* Python loads `C/model.dll` or `C/libmodel.so` via `ctypes`; if not found, uses a NumPy RK4 fallback.
* `run_all.py` iterates desserts → runs sim → saves plots & metrics.

---

## Troubleshooting

**Figures/CSV not appearing**
Ensure `py/plotting.py` and `py/metrics.py` use repo-relative paths (already fixed). Run from repo root:

```bash
python py/run_all.py
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
