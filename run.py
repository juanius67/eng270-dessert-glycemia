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

import requests
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

DEFAULT_DT = 0.5
DEFAULT_T_END = 1440.0

Vd_dL = 100.0  # glucose distribution (dL)
f_app_base = 0.35  # baseline appearance fraction
beta_fiber = 0.04  # ↓ appearance per 10 g fiber
beta_fat = 0.03    # gastric emptying slowdown per 10 g fat
kfast_base = 0.50  # min^-1 (fast sugars)
kslow_base = 0.08  # min^-1 (slow starch)
alpha_prot = 0.12  # μU·mL^-1 per g protein (insulin drive)
kprot = 0.05       # min^-1

TARGET_PEAK = 50.0
TARGET_TOL = 5.0
MAX_CAL_STEPS = 4
FACTOR_MIN, FACTOR_MAX = 0.1, 200.0

MACRO_FIELDS = {
    "carbs_g": "carbohydrates_100g",
    "sugars_g": "sugars_100g",
    "fiber_g": "fiber_100g",
    "fat_g": "fat_100g",
    "protein_g": "proteins_100g",
}


@dataclass
class DoseProfile:
    Afast: float
    kfast: float
    Aslow: float
    kslow: float
    Aprot: float
    kprot: float

    def scaled(self, factor: float) -> "DoseProfile":
        return DoseProfile(
            self.Afast * factor,
            self.kfast,
            self.Aslow * factor,
            self.kslow,
            self.Aprot,
            self.kprot,
        )

    def total_dose(self) -> float:
        dose_fast = self.Afast / self.kfast if self.kfast > 0.0 else 0.0
        dose_slow = self.Aslow / self.kslow if self.kslow > 0.0 else 0.0
        return dose_fast + dose_slow

    def total_amplitude(self) -> float:
        return self.Afast + self.Aslow

    def harmonic_rate(self) -> float:
        dose_fast = self.Afast / self.kfast if self.kfast > 0.0 else 0.0
        dose_slow = self.Aslow / self.kslow if self.kslow > 0.0 else 0.0
        total_dose = dose_fast + dose_slow
        denom = 0.0
        if dose_fast > 0.0 and self.kfast > 0.0:
            denom += dose_fast / self.kfast
        if dose_slow > 0.0 and self.kslow > 0.0:
            denom += dose_slow / self.kslow
        if total_dose > 0.0 and denom > 0.0:
            return total_dose / denom
        if self.Afast > 0.0 and self.kfast > 0.0:
            return self.kfast
        if self.Aslow > 0.0 and self.kslow > 0.0:
            return self.kslow
        return 1.0


@dataclass
class SimulationResult:
    times: List[float]
    glucose: List[float]
    insulin: List[float]
    profile: DoseProfile


class Simulator:
    def __init__(self) -> None:
        _, lib = load_backend()
        self._lib = lib
        self._params_struct = None
        self._has_extended = False

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
            if hasattr(self._lib, "simulate_ex"):
                self._lib.simulate_ex.argtypes = [
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.c_int,
                    ctypes.c_double,
                    ctypes.POINTER(Params),
                    ctypes.c_double,
                    ctypes.c_double,
                    ctypes.c_double,
                    ctypes.c_double,
                    ctypes.c_double,
                    ctypes.c_double,
                ]
                self._lib.simulate_ex.restype = None
                self._has_extended = True

    @property
    def label(self) -> str:
        if self._lib is None:
            return "python"
        return "c-ex" if self._has_extended else "c"

    def run(self, params: Dict[str, float], dt: float, t_end: float, profile: DoseProfile) -> SimulationResult:
        if self._lib is not None and self._params_struct is not None:
            if self._has_extended:
                return self._run_c_extended(params, dt, t_end, profile)
            return self._run_c_legacy(params, dt, t_end, profile)
        return self._run_python(params, dt, t_end, profile)

    # --- backends ---------------------------------------------------------
    def _run_c_extended(self, params: Dict[str, float], dt: float, t_end: float, profile: DoseProfile) -> SimulationResult:
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
            0.0,
            0.0,
        )
        self._lib.simulate_ex(
            g_arr,
            i_arr,
            nsteps,
            ctypes.c_double(dt),
            ctypes.byref(params_struct),
            ctypes.c_double(profile.Afast),
            ctypes.c_double(profile.kfast),
            ctypes.c_double(profile.Aslow),
            ctypes.c_double(profile.kslow),
            ctypes.c_double(profile.Aprot),
            ctypes.c_double(profile.kprot),
        )
        glucose = [g_arr[i] for i in range(nsteps)]
        insulin = [i_arr[i] for i in range(nsteps)]
        return SimulationResult(times, glucose, insulin, profile)

    def _run_c_legacy(self, params: Dict[str, float], dt: float, t_end: float, profile: DoseProfile) -> SimulationResult:
        assert self._lib is not None and self._params_struct is not None
        steps = int(round(t_end / dt))
        nsteps = steps + 1
        times = [i * dt for i in range(nsteps)]

        g_arr = (ctypes.c_double * nsteps)()
        i_arr = (ctypes.c_double * nsteps)()
        keff = profile.harmonic_rate()
        params_struct = self._params_struct(
            params["p1"],
            params["p2"],
            params["p3"],
            params["p4"],
            params["p5"],
            params["p6"],
            params["Gb"],
            params["Ib"],
            profile.total_amplitude(),
            keff,
        )
        self._lib.simulate(
            g_arr,
            i_arr,
            nsteps,
            ctypes.c_double(dt),
            ctypes.byref(params_struct),
        )
        glucose = [g_arr[i] for i in range(nsteps)]
        insulin = [i_arr[i] for i in range(nsteps)]
        return SimulationResult(times, glucose, insulin, profile)

    def _run_python(self, params: Dict[str, float], dt: float, t_end: float, profile: DoseProfile) -> SimulationResult:
        steps = int(round(t_end / dt))
        nsteps = steps + 1
        times = [i * dt for i in range(nsteps)]

        G = params["Gb"]
        X = 0.0
        I = params["Ib"]
        glucose: List[float] = []
        insulin: List[float] = []

        def deriv(t: float, g: float, x: float, ins: float) -> Tuple[float, float, float]:
            fast = profile.Afast * math.exp(-profile.kfast * t)
            slow = profile.Aslow * math.exp(-profile.kslow * t)
            D = fast + slow
            secretion = params["p4"] * max(0.0, g - params["p5"])
            iprot = profile.Aprot * math.exp(-profile.kprot * t)
            dG = -(params["p1"] + x) * g + params["p1"] * params["Gb"] + D
            dX = -params["p2"] * x + params["p3"] * (ins - params["Ib"])
            dI = -params["p6"] * (ins - params["Ib"]) + secretion + iprot
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

        return SimulationResult(times, glucose, insulin, profile)


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
        data = yaml.safe_load(handle) or {}
    required = {"p1", "p2", "p3", "p4", "p5", "p6", "Gb", "Ib"}
    missing = required.difference(data)
    if missing:
        raise KeyError(f"Missing config keys: {sorted(missing)}")
    params = {key: float(data[key]) for key in required}
    params["dt"] = float(data.get("dt", DEFAULT_DT))
    params["t_end"] = float(data.get("t_end", DEFAULT_T_END))
    return params


def clip(value: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, value))


def load_desserts(use_nutrition: bool) -> Dict[str, Dict[str, float | str | None | DoseProfile]]:
    if not use_nutrition:
        desserts: Dict[str, Dict[str, float | str | None | DoseProfile]] = {}
        for name, specs in LEGACY_DESSERTS.items():
            k = specs["k"]
            dose = specs["dose_mgdL"]
            profile = DoseProfile(k * dose, k, 0.0, k, 0.0, kprot)
            desserts[name] = {
                "profile": profile,
                "dose_mgdL": profile.total_dose(),
                "carbs_g": None,
                "sugars_g": None,
                "fiber_g": None,
                "fat_g": None,
                "protein_g": None,
                "f_fast": None,
                "f_app": None,
                "k_mod": None,
                "portion_g": None,
                "barcode": None,
            }
        return desserts

    with DESSERTS_CONFIG_PATH.open("r", encoding="utf-8") as handle:
        raw_data = yaml.safe_load(handle) or {}

    desserts: Dict[str, Dict[str, float | str | None | DoseProfile]] = {}
    for name, data in raw_data.items():
        carbs_g = float(data.get("carbs_g", 0.0))
        sugars_g = float(data.get("sugars_g", 0.0))
        fiber_g = float(data.get("fiber_g", 0.0))
        fat_g = float(data.get("fat_g", 0.0))
        protein_g = float(data.get("protein_g", 0.0))

        avail_carbs_g = max(0.0, carbs_g - 0.5 * fiber_g)
        f_fast = 0.0 if carbs_g <= 0.0 else clip(sugars_g / carbs_g, 0.0, 1.0)
        f_app = clip(f_app_base * (1.0 - beta_fiber * (fiber_g / 10.0)), 0.05, 0.60)
        k_mod = 1.0 / (1.0 + beta_fat * (fat_g / 10.0) + beta_fiber * (fiber_g / 10.0))
        kfast = kfast_base * k_mod
        kslow = kslow_base * k_mod
        dose_total_mgdL = (avail_carbs_g * f_app * 1000.0) / Vd_dL
        dose_fast = dose_total_mgdL * f_fast
        dose_slow = dose_total_mgdL * (1.0 - f_fast)
        Afast = kfast * dose_fast
        Aslow = kslow * dose_slow
        Aprot = alpha_prot * protein_g

        profile = DoseProfile(Afast, kfast, Aslow, kslow, Aprot, kprot)
        desserts[name] = {
            "profile": profile,
            "dose_mgdL": profile.total_dose(),
            "carbs_g": carbs_g,
            "sugars_g": sugars_g,
            "fiber_g": fiber_g,
            "fat_g": fat_g,
            "protein_g": protein_g,
            "f_fast": f_fast,
            "f_app": f_app,
            "k_mod": k_mod,
            "portion_g": float(data.get("portion_g", 0.0)) if data.get("portion_g") is not None else None,
            "barcode": data.get("barcode"),
        }
    return desserts


def clamp(value: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, value))


def calibrate_profile(
    simulator: Simulator,
    params: Dict[str, float],
    dt: float,
    t_end: float,
    base_profile: DoseProfile,
) -> Tuple[SimulationResult, bool]:
    if base_profile.Afast <= 0.0 and base_profile.Aslow <= 0.0:
        result = simulator.run(params, dt, t_end, base_profile)
        return result, False

    factor = 1.0
    result: SimulationResult | None = None
    for step in range(1, MAX_CAL_STEPS + 1):
        profile = base_profile.scaled(factor)
        result = simulator.run(params, dt, t_end, profile)
        peak_delta = max(result.glucose) - params["Gb"]
        print(f"  calibration iter {step}: scale={factor:.3f}, peakΔG={peak_delta:.3f}")
        if abs(peak_delta - TARGET_PEAK) <= TARGET_TOL:
            return result, True
        if peak_delta <= 1e-8:
            factor = clamp(factor * 1.5, FACTOR_MIN, FACTOR_MAX)
        else:
            factor = clamp(factor * (TARGET_PEAK / peak_delta), FACTOR_MIN, FACTOR_MAX)
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


def save_overlay(
    curves: Iterable[Tuple[str, Sequence[float], Sequence[float]]],
    path: Path,
    title: str,
    ylabel: str,
    xlim: Tuple[float, float] | None = None,
) -> None:
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


def make_dose_curves(profile: DoseProfile, duration: float = 360.0, step: float = 0.5) -> Tuple[List[float], List[float], List[float], List[float]]:
    n = int(round(duration / step)) + 1
    times = [i * step for i in range(n)]
    fast = [profile.Afast * math.exp(-profile.kfast * t) for t in times]
    slow = [profile.Aslow * math.exp(-profile.kslow * t) for t in times]
    total = [fast[i] + slow[i] for i in range(n)]
    return times, fast, slow, total


def save_d_components(name: str, profile: DoseProfile) -> Path:
    times, fast, slow, total = make_dose_curves(profile)
    path = FIG_DIR / f"D_components_{name}.png"
    plt.figure(figsize=(6.0, 3.4))
    plt.plot(times, fast, label="fast", lw=1.6)
    plt.plot(times, slow, label="slow", lw=1.6)
    plt.plot(times, total, label="total", lw=2.0, linestyle="--")
    plt.title(f"{name.title()} appearance components")
    plt.xlabel("Time (min)")
    plt.ylabel("D(t) (mg/dL·min⁻¹)")
    plt.xlim(0.0, 360.0)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()
    return path


def write_summary(rows: List[Dict[str, object]], latex: bool = False) -> None:
    base_fields = [
        "name",
        "backend",
        "mode",
        "peak_delta",
        "t_peakG",
        "AUCG_0_120",
        "peakI",
        "t_peakI",
        "AUCI_0_120",
        "dose_mgdL",
        "Afast",
        "Aslow",
        "Aprot",
        "kfast",
        "kslow",
        "kprot",
        "carbs_g",
        "sugars_g",
        "fiber_g",
        "fat_g",
        "protein_g",
        "f_fast",
        "f_app",
        "k_mod",
        "portion_g",
        "barcode",
    ]
    fieldnames = list(base_fields)
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
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
            ("mode", "mode", "{value}"),
            ("peak_delta", "peakΔG", "{value:.2f}"),
            ("t_peakG", "t$_{peak}$", "{value:.1f}"),
            ("AUCG_0_120", "AUCG", "{value:.1f}"),
            ("peakI", "peakI", "{value:.2f}"),
            ("t_peakI", "t$_{peakI}$", "{value:.1f}"),
        ]
        with latex_path.open("w", encoding="utf-8") as handle:
            handle.write("% Auto-generated summary table\n")
            handle.write(f"\\begin{{tabular}}{{l{'c' * (len(columns) - 1)}}}\\n")
            header = " & ".join(label for _, label, _ in columns)
            handle.write(f"{header} \\ \\hline\\n")
            for row in rows:
                formatted = []
                for key, _, fmt in columns:
                    value = row[key]
                    formatted.append(fmt.format(value=value))
                handle.write(" & ".join(formatted) + " \\ \n")
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
    desserts: Dict[str, Dict[str, float | str | None | DoseProfile]],
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
        base_profile = specs["profile"]
        if calibrate:
            result, converged = calibrate_profile(simulator, params, dt, t_end, base_profile)
            mode = "calibrated" if converged else "calibrated*"
        else:
            result = simulator.run(params, dt, t_end, base_profile)
            mode = "dose-driven"

        final_profile = result.profile
        full_mode = f"{nutrition_mode}-{mode}"

        metrics = compute_metrics(result, params)
        glucose_curves.append((name, result.times, result.glucose))
        insulin_curves.append((name, result.times, result.insulin))
        dose_times, _, _, dose_total = make_dose_curves(final_profile)
        dose_curves.append((name, dose_times, dose_total))

        final_dose = final_profile.total_dose()
        summary_row: Dict[str, object] = {
            "name": name,
            "backend": simulator.label,
            "mode": full_mode,
            "peak_delta": metrics["peak_delta"],
            "t_peakG": metrics["t_peakG"],
            "AUCG_0_120": metrics["AUCG"],
            "peakI": metrics["peakI"],
            "t_peakI": metrics["t_peakI"],
            "AUCI_0_120": metrics["AUCI"],
            "dose_mgdL": final_dose,
            "Afast": final_profile.Afast,
            "Aslow": final_profile.Aslow,
            "Aprot": final_profile.Aprot,
            "kfast": final_profile.kfast,
            "kslow": final_profile.kslow,
            "kprot": final_profile.kprot,
            "carbs_g": specs.get("carbs_g"),
            "sugars_g": specs.get("sugars_g"),
            "fiber_g": specs.get("fiber_g"),
            "fat_g": specs.get("fat_g"),
            "protein_g": specs.get("protein_g"),
            "f_fast": specs.get("f_fast"),
            "f_app": specs.get("f_app"),
            "k_mod": specs.get("k_mod"),
            "portion_g": specs.get("portion_g"),
            "barcode": specs.get("barcode"),
            "peakG": metrics["peakG"],
        }
        summary_rows.append(summary_row)

        print(
            f"{name:12s} | backend={simulator.label:6s} | mode={full_mode:18s} | "
            f"peakΔG={metrics['peak_delta']:+6.2f} mg/dL @ {metrics['t_peakG']:.1f} min | "
            f"AUCG={metrics['AUCG']:.1f} | peakI={metrics['peakI']:.2f}"
        )

        glucose_path = FIG_DIR / f"{name}_glucose.png"
        insulin_path = FIG_DIR / f"{name}_insulin.png"
        components_path = save_d_components(name, final_profile)

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
                {
                    "path": str(components_path.relative_to(ROOT)),
                    "caption": f"{name.title()} appearance components",
                },
            ]
        )

    overlay_glucose = FIG_DIR / "glucose_overlay.png"
    overlay_insulin = FIG_DIR / "insulin_overlay.png"
    overlay_dose = FIG_DIR / "D_overlay.png"

    save_overlay(glucose_curves, overlay_glucose, "Glucose overlay", "Glucose (mg/dL)")
    save_overlay(insulin_curves, overlay_insulin, "Insulin overlay", "Insulin (mU/L)")
    save_overlay(dose_curves, overlay_dose, "Dessert appearance D(t)", "Dose (mg/dL·min⁻¹)", xlim=(0.0, 180.0))

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
                "caption": "Dessert appearance overlay",
            },
        ]
    )

    write_summary(summary_rows, latex=latex)
    write_manifest(manifest_entries)


def run_sanity(simulator: Simulator, params: Dict[str, float], dt: float, t_end: float) -> None:
    print(f"Running sanity checks with backend={simulator.label}...")
    zero_profile = DoseProfile(0.0, 1.0, 0.0, 1.0, 0.0, kprot)
    baseline = simulator.run(params, dt, t_end, zero_profile)
    half = simulator.run(params, dt / 2.0, t_end, zero_profile)
    g_diff = abs(baseline.glucose[-1] - half.glucose[-1])
    i_diff = abs(baseline.insulin[-1] - half.insulin[-1])
    print(
        f"Baseline final G={baseline.glucose[-1]:.6f}, I={baseline.insulin[-1]:.6f}; "
        f"dt/2 final G={half.glucose[-1]:.6f}, I={half.insulin[-1]:.6f}"
    )
    print(f"Final-state diffs | ΔG={g_diff:.6e}, ΔI={i_diff:.6e}")


def fetch_off_product(barcode: str) -> Dict[str, object]:
    url = f"https://world.openfoodfacts.org/api/v2/product/{barcode}.json"
    try:
        response = requests.get(url, timeout=15)
    except requests.RequestException as exc:
        raise SystemExit(f"OpenFoodFacts request failed: {exc}") from exc
    if response.status_code != 200:
        raise SystemExit(f"OpenFoodFacts responded with status {response.status_code}")
    try:
        payload = response.json()
    except ValueError as exc:
        raise SystemExit(f"Invalid JSON from OpenFoodFacts: {exc}") from exc
    if payload.get("status") != 1:
        raise SystemExit(f"Product {barcode} not found on OpenFoodFacts")
    return payload.get("product", {})


def import_off_entry(name: str, barcode: str, portion_g_raw: str) -> None:
    try:
        portion_g = float(portion_g_raw)
    except ValueError as exc:
        raise SystemExit(f"Invalid portion (grams): {portion_g_raw}") from exc
    product = fetch_off_product(barcode)
    nutriments = product.get("nutriments", {})

    entry: Dict[str, object] = {
        "portion_g": portion_g,
        "barcode": barcode,
    }
    for key, field in MACRO_FIELDS.items():
        value = nutriments.get(field)
        try:
            per100 = float(value) if value is not None else 0.0
        except (TypeError, ValueError):
            per100 = 0.0
        scaled = round(per100 * portion_g / 100.0, 1)
        entry[key] = scaled

    if DESSERTS_CONFIG_PATH.exists():
        with DESSERTS_CONFIG_PATH.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
    else:
        data = {}
    data[name] = entry
    with DESSERTS_CONFIG_PATH.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=True)
    print(f"Imported {name} from OpenFoodFacts (barcode {barcode}) into {DESSERTS_CONFIG_PATH}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Dessert glycemia simulator")
    parser.add_argument("--all", action="store_true", help="Run full dessert pipeline")
    parser.add_argument("--sanity", action="store_true", help="Run baseline + dt-halving checks")
    parser.add_argument("--calibrate", action="store_true", help="Calibrate amplitudes to ~50 mg/dL peaks")
    parser.add_argument("--latex", action="store_true", help="Also emit LaTeX summary table")
    parser.add_argument("--no-nutrition", action="store_true", help="Disable nutrition-aware mapping")
    parser.add_argument(
        "--import-off",
        nargs=3,
        metavar=("NAME", "BARCODE", "PORTION_G"),
        help="Import nutrition from OpenFoodFacts using barcode",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)

    if args.import_off is not None:
        name, barcode, portion = args.import_off
        import_off_entry(name, barcode, portion)

    params = load_config(CONFIG_PATH)
    dt = params.pop("dt")
    t_end = params.pop("t_end")

    use_nutrition = not args.no_nutrition
    desserts = load_desserts(use_nutrition)

    simulator = Simulator()

    run_all = args.all or (args.import_off is None and not args.sanity)
    if run_all:
        run_pipeline(
            simulator,
            params,
            dt,
            t_end,
            desserts,
            calibrate=args.calibrate,
            latex=args.latex,
            nutrition_enabled=use_nutrition,
        )
    if args.sanity:
        run_sanity(simulator, params, dt, t_end)


if __name__ == "__main__":
    main(sys.argv[1:])
