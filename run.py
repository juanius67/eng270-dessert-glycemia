#!/usr/bin/env python3
"""Compact entrypoint for dessert-driven glycemia experiments."""
from __future__ import annotations

import argparse
import csv
import ctypes
import math
import platform
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import yaml

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
CONFIG_PATH = ROOT / "configs" / "params.yaml"
DESSERTS_CONFIG_PATH = ROOT / "configs" / "desserts.yaml"
C_DIR = ROOT / "C"
FIG_DIR = ROOT / "figures"
TABLE_DIR = ROOT / "tables"

LEGACY_DESSERTS: Dict[str, Dict[str, float]] = {
    "chocotorta": {"dose_mgdL": 60.0, "k": 0.08},
    "brigadeiro": {"dose_mgdL": 85.0, "k": 0.50},
    "alfajor": {"dose_mgdL": 70.0, "k": 0.10},
    "acai": {"dose_mgdL": 90.0, "k": 0.15},
}

Vd_dL = 100.0
f_app_base = 0.35
beta_fiber = 0.04
beta_fat = 0.03
kfast_base = 0.50
kslow_base = 0.08

TARGET_PEAK = 50.0
TARGET_TOL = 5.0
MAX_CAL_STEPS = 4
A_MIN, A_MAX = 0.1, 200.0


@dataclass
class SimulationResult:
    times: List[float]
    glucose: List[float]
    insulin: List[float]
    A: float
    k: float


class Simulator:
    def __init__(self) -> None:
        self.backend, self._lib = load_backend()
        if self._lib is not None:
            class Params(ctypes.Structure):
                _fields_ = [
                    ("p1", ctypes.c_double),
                    ("p2", ctypes.c_double),
                    ("p3", ctypes.c_double),
                    ("p4", ctypes.c_double),
                    ("p5", ctypes.c_double),
                    ("p6", ctypes.c_double),
                    ("Gb", ctypes.c_double),
                    ("Ib", ctypes.c_double),
                    ("A", ctypes.c_double),
                    ("k", ctypes.c_double),
                ]

            self._params_struct = Params
            self._lib.simulate.argtypes = [
                ctypes.POINTER(ctypes.c_double),
                ctypes.POINTER(ctypes.c_double),
                ctypes.c_int,
                ctypes.c_double,
                ctypes.POINTER(Params),
            ]
            self._lib.simulate.restype = None
        else:
            self._params_struct = None

    @property
    def label(self) -> str:
        return "c" if self._lib is not None else "python"

    def run(self, params: Dict[str, float], dt: float, t_end: float, A: float, k: float) -> SimulationResult:
        if self._lib is not None:
            return self._run_c(params, dt, t_end, A, k)
        return self._run_python(params, dt, t_end, A, k)

    # --- backends ---------------------------------------------------------
    def _run_c(self, params: Dict[str, float], dt: float, t_end: float, A: float, k: float) -> SimulationResult:
        assert self._lib is not None and self._params_struct is not None
        steps = int(round(t_end / dt))
        nsteps = steps + 1
        times = [i * dt for i in range(nsteps)]

        g_arr = (ctypes.c_double * nsteps)()
        i_arr = (ctypes.c_double * nsteps)()
        params_struct = self._params_struct(
            params["p1"],
            params["p2"],
            params["p3"],
            params["p4"],
            params["p5"],
            params["p6"],
            params["Gb"],
            params["Ib"],
            A,
            k,
        )
        self._lib.simulate(g_arr, i_arr, nsteps, ctypes.c_double(dt), ctypes.byref(params_struct))
        glucose = [g_arr[i] for i in range(nsteps)]
        insulin = [i_arr[i] for i in range(nsteps)]
        return SimulationResult(times, glucose, insulin, A, k)

    def _run_python(self, params: Dict[str, float], dt: float, t_end: float, A: float, k: float) -> SimulationResult:
        steps = int(round(t_end / dt))
        nsteps = steps + 1
        times = [i * dt for i in range(nsteps)]

        G = params["Gb"]
        X = 0.0
        I = params["Ib"]
        glucose: List[float] = []
        insulin: List[float] = []

        def deriv(t: float, g: float, x: float, ins: float) -> Tuple[float, float, float]:
            Dt = A * math.exp(-k * t)
            dG = -(params["p1"] + x) * g + params["p1"] * params["Gb"] + Dt
            dX = -params["p2"] * x + params["p3"] * (ins - params["Ib"])
            secretion = params["p4"] * max(0.0, g - params["p5"])
            dI = -params["p6"] * (ins - params["Ib"]) + secretion
            return dG, dX, dI

        for idx, current_time in enumerate(times):
            glucose.append(G)
            insulin.append(I)
            if idx == nsteps - 1:
                break
            k1 = deriv(current_time, G, X, I)
            G1 = G + 0.5 * dt * k1[0]
            X1 = X + 0.5 * dt * k1[1]
            I1 = I + 0.5 * dt * k1[2]

            k2 = deriv(current_time + 0.5 * dt, G1, X1, I1)
            G2 = G + 0.5 * dt * k2[0]
            X2 = X + 0.5 * dt * k2[1]
            I2 = I + 0.5 * dt * k2[2]

            k3 = deriv(current_time + 0.5 * dt, G2, X2, I2)
            G3 = G + dt * k3[0]
            X3 = X + dt * k3[1]
            I3 = I + dt * k3[2]

            k4 = deriv(current_time + dt, G3, X3, I3)

            G += (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
            X += (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
            I += (dt / 6.0) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2])

        return SimulationResult(times, glucose, insulin, A, k)


def load_backend() -> Tuple[str, ctypes.CDLL | None]:
    system = platform.system().lower()
    if system.startswith("win"):
        lib_path = C_DIR / "model.dll"
        build_cmds = [
            ["cl", "/nologo", "/LD", "C\\model.c", "/Fe:C\\model.dll"],
            ["gcc", "-O3", "-shared", "-o", "C\\model.dll", "C\\model.c"],
        ]
    elif system == "darwin":
        lib_path = C_DIR / "libmodel.dylib"
        build_cmds = [["clang", "-O3", "-fPIC", "-shared", "-o", "C/libmodel.dylib", "C/model.c"]]
    else:
        lib_path = C_DIR / "libmodel.so"
        build_cmds = [["gcc", "-O3", "-fPIC", "-shared", "-o", "C/libmodel.so", "C/model.c"]]

    lib = try_load_library(lib_path)
    if lib is None and not lib_path.exists():
        for cmd in build_cmds:
            print("Build attempt:", " ".join(cmd))
            try:
                subprocess.run(cmd, check=True, cwd=ROOT)
            except (subprocess.CalledProcessError, FileNotFoundError, OSError) as exc:
                print(f"  build failed: {exc}")
                continue
            lib = try_load_library(lib_path)
            if lib is not None:
                break
    if lib is None:
        print("Using pure-Python RK4 backend.")
        return "python", None
    print(f"Using C backend: {lib_path}")
    return "c", lib


def try_load_library(path: Path) -> ctypes.CDLL | None:
    if not path.exists():
        return None
    try:
        return ctypes.CDLL(str(path.resolve()))
    except OSError as exc:
        print(f"  load failed for {path}: {exc}")
        return None


def load_config(path: Path) -> Dict[str, float]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    required = {"p1", "p2", "p3", "p4", "p5", "p6", "Gb", "Ib", "dt", "t_end"}
    missing = required.difference(data)
    if missing:
        raise KeyError(f"Missing config keys: {sorted(missing)}")
    return {key: float(data[key]) for key in required}


def load_desserts(use_nutrition: bool) -> Dict[str, Dict[str, float]]:
    if not use_nutrition:
        desserts: Dict[str, Dict[str, float]] = {}
        for name, specs in LEGACY_DESSERTS.items():
            desserts[name] = {
                "dose_mgdL": specs["dose_mgdL"],
                "k": specs["k"],
                "carbs_g": None,
                "sugars_g": None,
                "fiber_g": None,
                "fat_g": None,
                "protein_g": None,
                "f_fast": None,
                "f_app": None,
                "k_mod": None,
                "keff": specs["k"],
            }
        return desserts

    with DESSERTS_CONFIG_PATH.open("r", encoding="utf-8") as handle:
        raw_data = yaml.safe_load(handle) or {}

    desserts: Dict[str, Dict[str, float]] = {}
    for name, data in raw_data.items():
        carbs_g = float(data.get("carbs_g", 0.0))
        sugars_g = float(data.get("sugars_g", 0.0))
        fiber_g = float(data.get("fiber_g", 0.0))
        fat_g = float(data.get("fat_g", 0.0))
        protein_g = float(data.get("protein_g", 0.0))

        avail_carbs_g = max(0.0, carbs_g - 0.5 * fiber_g)
        f_fast = 0.0 if carbs_g <= 0.0 else min(1.0, max(0.0, sugars_g / carbs_g))
        f_app = f_app_base * (1.0 - beta_fiber * (fiber_g / 10.0))
        f_app = max(0.05, min(0.60, f_app))
        k_mod = 1.0 / (1.0 + beta_fat * (fat_g / 10.0) + beta_fiber * (fiber_g / 10.0))
        keff = (f_fast * kfast_base + (1.0 - f_fast) * kslow_base) * k_mod
        dose_total_mgdL = (avail_carbs_g * f_app * 1000.0) / Vd_dL

        desserts[name] = {
            "dose_mgdL": dose_total_mgdL,
            "k": keff,
            "carbs_g": carbs_g,
            "sugars_g": sugars_g,
            "fiber_g": fiber_g,
            "fat_g": fat_g,
            "protein_g": protein_g,
            "f_fast": f_fast,
            "f_app": f_app,
            "k_mod": k_mod,
            "keff": keff,
        }
    return desserts


def clamp(value: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, value))


def calibrate_amplitude(
    simulator: Simulator,
    params: Dict[str, float],
    dt: float,
    t_end: float,
    base_A: float,
    k: float,
) -> Tuple[SimulationResult, bool]:
    A = clamp(base_A, A_MIN, A_MAX)
    result: SimulationResult | None = None
    for step in range(1, MAX_CAL_STEPS + 1):
        result = simulator.run(params, dt, t_end, A, k)
        peak_delta = max(result.glucose) - params["Gb"]
        print(f"  calibration iter {step}: A={A:.3f}, peakΔG={peak_delta:.3f}")
        if abs(peak_delta - TARGET_PEAK) <= TARGET_TOL:
            return result, True
        if peak_delta <= 1e-8:
            A = clamp(A * 1.5, A_MIN, A_MAX)
        else:
            A = clamp(A * (TARGET_PEAK / peak_delta), A_MIN, A_MAX)
    assert result is not None
    return result, False


def integrate_trapezoid(
    times: Sequence[float],
    values: Sequence[float],
    limit: float,
    subtract: float | None = None,
    clamp_zero: bool = False,
) -> float:
    total = 0.0
    for idx in range(1, len(times)):
        t0 = times[idx - 1]
        t1 = times[idx]
        if t0 >= limit:
            break
        span_end = min(t1, limit)
        if span_end <= t0:
            continue
        v0 = values[idx - 1]
        v1 = values[idx]
        if subtract is not None:
            v0 -= subtract
            v1 -= subtract
        if clamp_zero:
            v0 = max(0.0, v0)
            v1 = max(0.0, v1)
        total += 0.5 * (v0 + v1) * (span_end - t0)
        if t1 >= limit:
            break
    return total


def compute_metrics(result: SimulationResult, params: Dict[str, float]) -> Dict[str, float]:
    peakG = max(result.glucose)
    t_peakG = result.times[result.glucose.index(peakG)]
    peakI = max(result.insulin)
    t_peakI = result.times[result.insulin.index(peakI)]
    peak_delta = peakG - params["Gb"]
    aucg = integrate_trapezoid(result.times, result.glucose, 120.0, subtract=params["Gb"], clamp_zero=True)
    auci = integrate_trapezoid(result.times, result.insulin, 120.0)
    return {
        "peakG": peakG,
        "t_peakG": t_peakG,
        "peak_delta": peak_delta,
        "peakI": peakI,
        "t_peakI": t_peakI,
        "AUCG": aucg,
        "AUCI": auci,
    }


def ensure_directories() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)


def save_line_plot(x: Sequence[float], y: Sequence[float], path: Path, title: str, ylabel: str) -> None:
    plt.figure(figsize=(6.0, 3.4))
    plt.plot(x, y, lw=1.8)
    plt.title(title)
    plt.xlabel("Time (min)")
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def save_overlay(curves: Iterable[Tuple[str, Sequence[float], Sequence[float]]], path: Path, title: str, ylabel: str, xlim: Tuple[float, float] | None = None) -> None:
    plt.figure(figsize=(6.0, 3.4))
    for name, x, y in curves:
        plt.plot(x, y, lw=1.5, label=name)
    plt.title(title)
    plt.xlabel("Time (min)")
    plt.ylabel(ylabel)
    if xlim is not None:
        plt.xlim(*xlim)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def make_dose_curve(A: float, k: float, duration: float = 60.0, step: float = 0.1) -> Tuple[List[float], List[float]]:
    n = int(round(duration / step)) + 1
    times = [i * step for i in range(n)]
    values = [A * math.exp(-k * t) for t in times]
    return times, values


def write_summary(rows: List[Dict[str, object]], latex: bool = False) -> None:
    fieldnames = [
        "name",
        "k",
        "backend",
        "mode",
        "dose_mgdL",
        "final_A",
        "peakG",
        "t_peakG",
        "AUCG_aboveGb_0_120",
        "peakI",
        "t_peakI",
        "AUCI_0_120",
    ]
    extra_fields = [
        "carbs_g",
        "sugars_g",
        "fiber_g",
        "fat_g",
        "protein_g",
        "f_fast",
        "f_app",
        "k_mod",
        "keff",
        "dose_mgdL",
        "final_A",
    ]
    for field in extra_fields:
        if field not in fieldnames:
            fieldnames.append(field)
    with (TABLE_DIR / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"Wrote {TABLE_DIR / 'summary.csv'}")

    if latex:
        latex_path = TABLE_DIR / "summary.tex"
        columns = [
            ("name", "name", "{value}"),
            ("k", "k", "{value:.3f}"),
            ("backend", "backend", "{value}"),
            ("mode", "mode", "{value}"),
            ("dose_mgdL", "dose", "{value:.2f}"),
            ("final_A", "A", "{value:.2f}"),
            ("peakG", "peakG", "{value:.2f}"),
            ("t_peakG", "t$_{peakG}$", "{value:.1f}"),
            ("AUCG_aboveGb_0_120", "AUCG", "{value:.1f}"),
            ("peakI", "peakI", "{value:.2f}"),
            ("t_peakI", "t$_{peakI}$", "{value:.1f}"),
            ("AUCI_0_120", "AUCI", "{value:.1f}"),
        ]
        with latex_path.open("w", encoding="utf-8") as handle:
            handle.write("% Auto-generated summary table\\n")
            handle.write(f"\\begin{{tabular}}{{l{'c' * (len(columns) - 1)}}}\\n")
            header = " & ".join(label for _, label, _ in columns)
            handle.write(f"{header} \\\\ \\hline\\n")
            for row in rows:
                formatted = []
                for key, _, fmt in columns:
                    value = row[key]
                    formatted.append(fmt.format(value=value))
                handle.write(" & ".join(formatted) + " \\\\ \n")
            handle.write("\\end{tabular}\\n")
        print(f"Wrote {latex_path}")


def write_manifest(entries: List[Dict[str, str]]) -> None:
    manifest_path = TABLE_DIR / "manifest.csv"
    fieldnames = ["path", "caption"]
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for entry in entries:
            writer.writerow(entry)
    print(f"Wrote {manifest_path}")


def run_pipeline(
    simulator: Simulator,
    params: Dict[str, float],
    dt: float,
    t_end: float,
    desserts: Dict[str, Dict[str, float]],
    calibrate: bool,
    latex: bool,
    nutrition_enabled: bool,
) -> None:
    ensure_directories()
    summary_rows: List[Dict[str, object]] = []
    glucose_curves: List[Tuple[str, List[float], List[float]]] = []
    insulin_curves: List[Tuple[str, List[float], List[float]]] = []
    dose_curves: List[Tuple[str, List[float], List[float]]] = []
    manifest_entries: List[Dict[str, str]] = []

    mode_label = "calibrated" if calibrate else "dose-driven"
    print(f"Running dessert pipeline ({mode_label}) with backend={simulator.label}...")

    nutrition_mode = "nutrition" if nutrition_enabled else "legacy"

    for name, specs in desserts.items():
        base_A = specs["k"] * specs["dose_mgdL"]
        if calibrate:
            result, converged = calibrate_amplitude(simulator, params, dt, t_end, base_A, specs["k"])
            mode = "calibrated" if converged else "calibrated*"
        else:
            result = simulator.run(params, dt, t_end, base_A, specs["k"])
            mode = "dose-driven"

        full_mode = f"{nutrition_mode}-{mode}"

        metrics = compute_metrics(result, params)
        glucose_curves.append((name, result.times, result.glucose))
        insulin_curves.append((name, result.times, result.insulin))
        dose_curves.append((name, *make_dose_curve(result.A, result.k)))

        summary_row: Dict[str, object] = {
            "name": name,
            "k": specs["k"],
            "backend": simulator.label,
            "mode": full_mode,
            "dose_mgdL": specs["dose_mgdL"],
            "final_A": result.A,
            "peakG": metrics["peakG"],
            "t_peakG": metrics["t_peakG"],
            "AUCG_aboveGb_0_120": metrics["AUCG"],
            "peakI": metrics["peakI"],
            "t_peakI": metrics["t_peakI"],
            "AUCI_0_120": metrics["AUCI"],
            "carbs_g": specs.get("carbs_g"),
            "sugars_g": specs.get("sugars_g"),
            "fiber_g": specs.get("fiber_g"),
            "fat_g": specs.get("fat_g"),
            "protein_g": specs.get("protein_g"),
            "f_fast": specs.get("f_fast"),
            "f_app": specs.get("f_app"),
            "k_mod": specs.get("k_mod"),
            "keff": specs.get("keff"),
        }
        summary_rows.append(summary_row)

        print(
            f"{name:12s} | backend={simulator.label:6s} | mode={full_mode:18s} | "
            f"keff={specs['k']:.3f} | dose={specs['dose_mgdL']:.2f} mg/dL | "
            f"peakΔG={metrics['peak_delta']:+6.2f} mg/dL @ {metrics['t_peakG']:.1f} min"
        )

        glucose_path = FIG_DIR / f"{name}_glucose.png"
        insulin_path = FIG_DIR / f"{name}_insulin.png"

        save_line_plot(result.times, result.glucose, glucose_path, f"{name.title()} glucose", "Glucose (mg/dL)")
        save_line_plot(result.times, result.insulin, insulin_path, f"{name.title()} insulin", "Insulin (mU/L)")

        manifest_entries.extend(
            [
                {
                    "path": str(glucose_path.relative_to(ROOT)),
                    "caption": f"{name.title()} glucose curve",
                },
                {
                    "path": str(insulin_path.relative_to(ROOT)),
                    "caption": f"{name.title()} insulin curve",
                },
            ]
        )

    overlay_glucose = FIG_DIR / "glucose_overlay.png"
    overlay_insulin = FIG_DIR / "insulin_overlay.png"
    overlay_dose = FIG_DIR / "D_overlay.png"

    save_overlay(glucose_curves, overlay_glucose, "Glucose overlay", "Glucose (mg/dL)")
    save_overlay(insulin_curves, overlay_insulin, "Insulin overlay", "Insulin (mU/L)")
    save_overlay(dose_curves, overlay_dose, "Dessert input D(t)=A e^{-k t}", "Dose (mg/dL)", xlim=(0.0, 60.0))

    manifest_entries.extend(
        [
            {
                "path": str(overlay_glucose.relative_to(ROOT)),
                "caption": "Dessert glucose overlay",
            },
            {
                "path": str(overlay_insulin.relative_to(ROOT)),
                "caption": "Dessert insulin overlay",
            },
            {
                "path": str(overlay_dose.relative_to(ROOT)),
                "caption": "Dessert input profiles",
            },
        ]
    )

    write_summary(summary_rows, latex=latex)
    write_manifest(manifest_entries)


def run_sanity(simulator: Simulator, params: Dict[str, float], dt: float, t_end: float) -> None:
    print(f"Running sanity checks with backend={simulator.label}...")
    baseline = simulator.run(params, dt, t_end, 0.0, 0.0)
    half = simulator.run(params, dt / 2.0, t_end, 0.0, 0.0)
    g_diff = abs(baseline.glucose[-1] - half.glucose[-1])
    i_diff = abs(baseline.insulin[-1] - half.insulin[-1])
    print(
        f"Baseline final G={baseline.glucose[-1]:.6f}, I={baseline.insulin[-1]:.6f}; "
        f"dt/2 final G={half.glucose[-1]:.6f}, I={half.insulin[-1]:.6f}"
    )
    print(f"Final-state diffs | ΔG={g_diff:.6e}, ΔI={i_diff:.6e}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Dessert glycemia simulator")
    parser.add_argument("--all", action="store_true", help="Run full dessert pipeline (default)")
    parser.add_argument("--sanity", action="store_true", help="Run baseline + dt-halving checks")
    parser.add_argument("--calibrate", action="store_true", help="Calibrate amplitudes to ~50 mg/dL peaks")
    parser.add_argument("--latex", action="store_true", help="Also emit LaTeX summary table")
    parser.add_argument(
        "--nutrition",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable nutrition-aware dessert mapping (default: on)",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    params = load_config(CONFIG_PATH)
    dt = params["dt"]
    t_end = params["t_end"]

    desserts = load_desserts(args.nutrition)

    simulator = Simulator()

    run_all = args.all or not args.sanity
    if run_all:
        run_pipeline(
            simulator,
            params,
            dt,
            t_end,
            desserts,
            calibrate=args.calibrate,
            latex=args.latex,
            nutrition_enabled=args.nutrition,
        )
    if args.sanity:
        run_sanity(simulator, params, dt, t_end)


if __name__ == "__main__":
    main(sys.argv[1:])
