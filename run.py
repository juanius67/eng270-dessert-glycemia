#!/usr/bin/env python3
"""One-stop simulation runner for dessert-driven glycemia experiments."""
from __future__ import annotations

import argparse
import csv
import ctypes
import math
import platform
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import yaml

# Use non-interactive backend for reproducibility
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CONFIG_PATH = Path("configs/params.yaml")
FIG_DIR = Path("figures")
TABLE_DIR = Path("tables")
C_DIR = Path("C")

DESSERTS: Dict[str, Dict[str, float]] = {
    "chocotorta": {"dose_mgdL": 60.0, "k": 0.08},
    "brigadeiro": {"dose_mgdL": 85.0, "k": 0.50},
    "alfajor": {"dose_mgdL": 70.0, "k": 0.10},
    "acai": {"dose_mgdL": 90.0, "k": 0.15},
}

TARGET_PEAK_DELTA = 50.0
TARGET_TOLERANCE = 5.0
MAX_CALIBRATION_STEPS = 4
A_MIN, A_MAX = 0.1, 200.0


def load_config(path: Path) -> Dict[str, float]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    required = {"p1", "p2", "p3", "p4", "p5", "p6", "Gb", "Ib", "dt", "t_end"}
    missing = required.difference(data)
    if missing:
        raise KeyError(f"Missing config keys: {sorted(missing)}")
    return {k: float(data[k]) for k in required}


def clamp(value: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, value))


@dataclass
class SimulationResult:
    times: List[float]
    glucose: List[float]
    insulin: List[float]
    A: float
    k: float


class Simulator:
    def __init__(self) -> None:
        self.backend, self._lib = self._load_backend()

    @property
    def backend_label(self) -> str:
        return "c" if self.backend == "c" else "python"

    def run(self, params: Dict[str, float], dt: float, t_end: float, A: float, k: float) -> SimulationResult:
        if self.backend == "c":
            return self._run_c(params, dt, t_end, A, k)
        return self._run_python(params, dt, t_end, A, k)

    # --- Backend management -------------------------------------------------
    def _load_backend(self) -> Tuple[str, Optional[ctypes.CDLL]]:
        system = platform.system().lower()
        if system.startswith("win"):
            lib_path = C_DIR / "model.dll"
            build_attempts = [
                ["cl", "/nologo", "/LD", str(C_DIR / "model.c"), f"/Fe:{C_DIR / 'model.dll'}"],
                ["gcc", "-O3", "-shared", "-o", str(C_DIR / "model.dll"), str(C_DIR / "model.c")],
            ]
        elif system == "darwin":
            lib_path = C_DIR / "libmodel.dylib"
            build_attempts = [["clang", "-O3", "-fPIC", "-shared", "-o", str(lib_path), str(C_DIR / "model.c")]]
        else:
            lib_path = C_DIR / "libmodel.so"
            build_attempts = [["gcc", "-O3", "-fPIC", "-shared", "-o", str(lib_path), str(C_DIR / "model.c")]]

        lib = self._try_load_library(lib_path)
        if lib is not None:
            print(f"Loaded C backend: {lib_path}")
            return "c", lib

        for cmd in build_attempts:
            try:
                print("Attempting build:", " ".join(cmd))
                subprocess.run(cmd, check=True, cwd=Path.cwd())
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                print(f"  Build failed: {exc}")
                continue
            lib = self._try_load_library(lib_path)
            if lib is not None:
                print(f"Built and loaded C backend: {lib_path}")
                return "c", lib

        print("Falling back to pure-Python RK4 backend.")
        return "python", None

    @staticmethod
    def _try_load_library(path: Path) -> Optional[ctypes.CDLL]:
        if not path.exists():
            return None
        try:
            return ctypes.CDLL(str(path.resolve()))
        except OSError as exc:
            print(f"  Failed to load {path}: {exc}")
            return None

    # --- C backend ----------------------------------------------------------
    def _run_c(self, params: Dict[str, float], dt: float, t_end: float, A: float, k: float) -> SimulationResult:
        assert self._lib is not None

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

        self._lib.simulate.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
            ctypes.c_double,
            ctypes.POINTER(Params),
        ]
        self._lib.simulate.restype = None

        steps = int(round(t_end / dt))
        nsteps = steps + 1
        times = [i * dt for i in range(nsteps)]

        g_arr = (ctypes.c_double * nsteps)()
        i_arr = (ctypes.c_double * nsteps)()
        params_struct = Params(
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

    # --- Python fallback ----------------------------------------------------
    def _run_python(self, params: Dict[str, float], dt: float, t_end: float, A: float, k: float) -> SimulationResult:
        steps = int(round(t_end / dt))
        nsteps = steps + 1
        times = [i * dt for i in range(nsteps)]

        G = params["Gb"]
        X = 0.0
        I = params["Ib"]

        glucose: List[float] = []
        insulin: List[float] = []

        def derivatives(t: float, g: float, x: float, ins: float) -> Tuple[float, float, float]:
            D_t = A * math.exp(-k * t)
            dG = - (params["p1"] + x) * g + params["p1"] * params["Gb"] + D_t
            dX = - params["p2"] * x + params["p3"] * (ins - params["Ib"])
            secretion = params["p4"] * max(0.0, g - params["p5"])
            dI = - params["p6"] * (ins - params["Ib"]) + secretion
            return dG, dX, dI

        for idx, current_time in enumerate(times):
            glucose.append(G)
            insulin.append(I)
            if idx == nsteps - 1:
                break
            k1 = derivatives(current_time, G, X, I)
            G1 = G + 0.5 * dt * k1[0]
            X1 = X + 0.5 * dt * k1[1]
            I1 = I + 0.5 * dt * k1[2]

            k2 = derivatives(current_time + 0.5 * dt, G1, X1, I1)
            G2 = G + 0.5 * dt * k2[0]
            X2 = X + 0.5 * dt * k2[1]
            I2 = I + 0.5 * dt * k2[2]

            k3 = derivatives(current_time + 0.5 * dt, G2, X2, I2)
            G3 = G + dt * k3[0]
            X3 = X + dt * k3[1]
            I3 = I + dt * k3[2]

            k4 = derivatives(current_time + dt, G3, X3, I3)

            G += (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
            X += (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
            I += (dt / 6.0) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2])

        return SimulationResult(times, glucose, insulin, A, k)


def compute_metrics(result: SimulationResult, params: Dict[str, float]) -> Dict[str, float]:
    Gb = params["Gb"]
    times = result.times
    glucose = result.glucose
    insulin = result.insulin

    peakG = max(glucose)
    t_peakG = times[glucose.index(peakG)]
    peakI = max(insulin)
    t_peakI = times[insulin.index(peakI)]

    delta_glucose = [max(g - Gb, 0.0) for g in glucose]
    aucG = trapezoidal_area(times, delta_glucose, upper_limit=120.0)
    aucI = trapezoidal_area(times, insulin, upper_limit=120.0)

    return {
        "peakG": peakG,
        "t_peakG": t_peakG,
        "peakI": peakI,
        "t_peakI": t_peakI,
        "AUCG": aucG,
        "AUCI": aucI,
        "peak_deltaG": max(delta_glucose),
    }


def trapezoidal_area(times: Iterable[float], values: Iterable[float], upper_limit: float) -> float:
    time_list = list(times)
    value_list = list(values)
    area = 0.0
    for idx in range(1, len(time_list)):
        t0 = time_list[idx - 1]
        t1 = time_list[idx]
        if t0 >= upper_limit:
            break
        v0 = value_list[idx - 1]
        v1 = value_list[idx]
        if t1 > upper_limit:
            if t1 == t0:
                continue
            frac = (upper_limit - t0) / (t1 - t0)
            v1 = v0 + (v1 - v0) * frac
            t1 = upper_limit
        area += 0.5 * (v0 + v1) * (t1 - t0)
        if t1 >= upper_limit:
            break
    return area


def ensure_directories() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)


def save_line_plot(times: List[float], values: List[float], path: Path, title: str, ylabel: str, color: Optional[str] = None) -> None:
    plt.figure(figsize=(6, 3))
    plt.plot(times, values, color=color)
    plt.title(title)
    plt.xlabel("Time (min)")
    plt.ylabel(ylabel)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def save_overlay(curves: List[Tuple[str, List[float], List[float]]], path: Path, title: str, ylabel: str, xlim: Optional[Tuple[float, float]] = None) -> None:
    plt.figure(figsize=(6, 3))
    for label, times, values in curves:
        plt.plot(times, values, label=label)
    plt.title(title)
    plt.xlabel("Time (min)")
    plt.ylabel(ylabel)
    plt.grid(True, alpha=0.3)
    if xlim:
        plt.xlim(*xlim)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def simulate_dessert(
    name: str,
    specs: Dict[str, float],
    simulator: Simulator,
    params: Dict[str, float],
    dt: float,
    t_end: float,
    calibrate: bool,
) -> Tuple[SimulationResult, Dict[str, float], str]:
    dose = specs["dose_mgdL"]
    k = specs["k"]
    initial_A = dose * k

    mode = "default"
    if calibrate:
        A, result = calibrate_amplitude(simulator, params, dt, t_end, initial_A, k, name)
        mode = "calibrated"
    else:
        A = initial_A
        result = simulator.run(params, dt, t_end, A, k)

    metrics = compute_metrics(result, params)
    print(
        f"{name:12s} | mode={mode:10s} | A={result.A:7.3f} | peak ΔG={metrics['peak_deltaG']:+6.2f} mg/dL @ {metrics['t_peakG']:.1f} min | "
        f"peak I={metrics['peakI']:.2f} @ {metrics['t_peakI']:.1f} min"
    )
    return result, metrics, mode


def calibrate_amplitude(
    simulator: Simulator,
    params: Dict[str, float],
    dt: float,
    t_end: float,
    initial_A: float,
    k: float,
    label: str,
) -> Tuple[float, SimulationResult]:
    A = clamp(initial_A, A_MIN, A_MAX)
    last_result = simulator.run(params, dt, t_end, A, k)
    metrics = compute_metrics(last_result, params)

    for iteration in range(1, MAX_CALIBRATION_STEPS + 1):
        peak_delta = metrics["peak_deltaG"]
        print(f"  [{label}] Iter {iteration}: A={A:.4f} -> peak ΔG={peak_delta:.3f}")
        if abs(peak_delta - TARGET_PEAK_DELTA) <= TARGET_TOLERANCE:
            break
        if peak_delta <= 1e-6:
            scale = 1.5
        else:
            scale = TARGET_PEAK_DELTA / peak_delta
        A = clamp(A * scale, A_MIN, A_MAX)
        last_result = simulator.run(params, dt, t_end, A, k)
        metrics = compute_metrics(last_result, params)
    return A, last_result


def run_all_desserts(simulator: Simulator, params: Dict[str, float], dt: float, t_end: float, calibrate: bool) -> None:
    ensure_directories()
    summary_rows: List[Dict[str, object]] = []
    glucose_curves: List[Tuple[str, List[float], List[float]]] = []
    insulin_curves: List[Tuple[str, List[float], List[float]]] = []
    dose_curves: List[Tuple[str, List[float], List[float]]] = []

    print("Running dessert simulations...")
    for name, specs in DESSERTS.items():
        result, metrics, mode = simulate_dessert(name, specs, simulator, params, dt, t_end, calibrate)

        save_line_plot(
            result.times,
            result.glucose,
            FIG_DIR / f"{name}_glucose.png",
            title=f"{name.title()} glucose",
            ylabel="Glucose (mg/dL)",
        )
        save_line_plot(
            result.times,
            result.insulin,
            FIG_DIR / f"{name}_insulin.png",
            title=f"{name.title()} insulin",
            ylabel="Insulin (mU/L)",
        )

        glucose_curves.append((name, result.times, result.glucose))
        insulin_curves.append((name, result.times, result.insulin))

        dose_times, dose_values = make_dose_curve(result.A, result.k, duration=60.0, step=0.1)
        dose_curves.append((name, dose_times, dose_values))

        summary_rows.append(
            {
                "name": name,
                "k": specs["k"],
                "backend": simulator.backend_label,
                "mode": mode,
                "dose_mgdL": specs["dose_mgdL"],
                "final_A": result.A,
                "peakG": metrics["peakG"],
                "t_peakG": metrics["t_peakG"],
                "AUCG_aboveGb_0_120": metrics["AUCG"],
                "peakI": metrics["peakI"],
                "t_peakI": metrics["t_peakI"],
                "AUCI_0_120": metrics["AUCI"],
            }
        )

    save_overlay(glucose_curves, FIG_DIR / "glucose_overlay.png", "Glucose overlay", "Glucose (mg/dL)")
    save_overlay(insulin_curves, FIG_DIR / "insulin_overlay.png", "Insulin overlay", "Insulin (mU/L)")
    save_overlay(dose_curves, FIG_DIR / "D_overlay.png", "Dessert input (A e^{-kt})", "Dose (mg/dL)", xlim=(0.0, 60.0))

    write_summary(summary_rows)


def make_dose_curve(A: float, k: float, duration: float, step: float) -> Tuple[List[float], List[float]]:
    nsteps = int(round(duration / step)) + 1
    times = [i * step for i in range(nsteps)]
    values = [A * math.exp(-k * t) for t in times]
    return times, values


def write_summary(rows: List[Dict[str, object]]) -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
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
    with (TABLE_DIR / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"Wrote {TABLE_DIR / 'summary.csv'}")


def run_sanity_checks(simulator: Simulator, params: Dict[str, float], dt: float, t_end: float) -> None:
    print("Running sanity checks (baseline + dt halving)...")
    baseline = simulator.run(params, dt, t_end, A=0.0, k=0.0)
    half_dt = simulator.run(params, dt / 2.0, t_end, A=0.0, k=0.0)
    g_diff = abs(baseline.glucose[-1] - half_dt.glucose[-1])
    i_diff = abs(baseline.insulin[-1] - half_dt.insulin[-1])
    print(
        f"Baseline final G={baseline.glucose[-1]:.4f}, I={baseline.insulin[-1]:.4f}; "
        f"dt/2 final G={half_dt.glucose[-1]:.4f}, I={half_dt.insulin[-1]:.4f}"
    )
    print(f"Final state differences | ΔG={g_diff:.6f}, ΔI={i_diff:.6f}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Dessert glycemia simulator")
    parser.add_argument("--all", action="store_true", help="Run full dessert pipeline (default)")
    parser.add_argument("--sanity", action="store_true", help="Run baseline and dt-halving checks")
    parser.add_argument("--calibrate", action="store_true", help="Calibrate dessert amplitudes")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(CONFIG_PATH)
    simulator = Simulator()

    dt = config["dt"]
    t_end = config["t_end"]

    run_all = args.all or not args.sanity
    if run_all:
        run_all_desserts(simulator, config, dt, t_end, calibrate=args.calibrate)
    if args.sanity:
        run_sanity_checks(simulator, config, dt, t_end)


if __name__ == "__main__":
    main()
