#!/usr/bin/env python3
"""Single CLI entry-point for nutrition-driven Bergman minimal model runs."""
from __future__ import annotations

import argparse
import csv
import json
import math
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import yaml

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ctypes import byref, c_double

from src import bindings
from src.bindings import get_build_metadata, set_params_from_dict

ROOT = Path(__file__).resolve().parent
CONFIGS_DIR = ROOT / "configs"
FROZEN_DIR = CONFIGS_DIR / "frozen"
DEFAULT_PARAMS_PATH = CONFIGS_DIR / "params.yaml"


@dataclass
class MealAppearance:
    A_fast: float
    k_fast: float
    A_slow: float
    k_slow: float
    A_prot: float
    k_prot: float
    dose_fast: float
    dose_slow: float
    dose_total: float


@dataclass
class SimulationResult:
    name: str
    times: np.ndarray
    glucose: np.ndarray
    insulin: np.ndarray
    d_fast: np.ndarray
    d_slow: np.ndarray
    u: np.ndarray


def load_params_yaml(path: Path = DEFAULT_PARAMS_PATH) -> Dict[str, float]:
    required_keys = {
        "Gb_mg_dL",
        "Ib_uU_mL",
        "S_G_min1",
        "p2_min1",
        "p3_min1_per_uU_mL",
        "n_min1",
        "Vd_dL",
        "f_hep",
        "f_app0",
        "beta_fiber_per10g",
        "beta_fat_per10g",
        "k_fast0_min1",
        "k_slow0_min1",
        "alpha_prot_uU_mL_per_g",
        "k_prot_min1",
        "dt_min",
        "t_end_min",
    }
    with Path(path).open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError("params.yaml must contain a mapping")

    legacy_map = {"p1": "S_G_min1", "p2": "p2_min1", "p3": "p3_min1_per_uU_mL", "p4": "n_min1"}
    used_legacy = False
    for legacy_key, new_key in legacy_map.items():
        if legacy_key in data and new_key not in data:
            data[new_key] = data[legacy_key]
            used_legacy = True
    if used_legacy:
        warnings.warn("Mapped legacy p1..p4 keys to unit-encoded names", RuntimeWarning, stacklevel=2)
    if "Gb_mg_dL" not in data:
        data["Gb_mg_dL"] = 90.0
    if "Ib_uU_mL" not in data:
        data["Ib_uU_mL"] = 7.0

    missing = sorted(required_keys - data.keys())
    if missing:
        raise KeyError(f"Missing parameters: {', '.join(missing)}")

    return {key: float(data[key]) for key in required_keys}


def _clip(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def compute_meal_appearance(dessert: Dict[str, float], params: Dict[str, float], calibrate: bool) -> MealAppearance:
    carbs = float(dessert.get("carbs_g", 0.0))
    sugars = float(dessert.get("sugars_g", 0.0))
    fiber = float(dessert.get("fiber_g", 0.0))
    fat = float(dessert.get("fat_g", 0.0))
    protein = float(dessert.get("protein_g", 0.0))

    C_avail = max(0.0, carbs - 0.5 * fiber)
    carbs_safe = carbs if carbs > 1e-9 else 1e-9
    f_fast = _clip(sugars / carbs_safe, 0.0, 1.0)
    f_app0 = params["f_app0"]
    beta_fiber = params["beta_fiber_per10g"] * fiber / 10.0
    f_app = _clip(f_app0 * (1.0 - beta_fiber), 0.05, 0.60)
    k_mod = 1.0 / (1.0 + params["beta_fat_per10g"] * fat / 10.0 + params["beta_fiber_per10g"] * fiber / 10.0)
    k_fast = params["k_fast0_min1"] * k_mod
    k_slow = params["k_slow0_min1"] * k_mod

    Vd = params["Vd_dL"]
    f_hep = params["f_hep"]
    dose_total = (1000.0 * C_avail / Vd) * f_app * (1.0 - f_hep)

    dose_fast = f_fast * dose_total
    dose_slow = (1.0 - f_fast) * dose_total
    if calibrate and dose_total > 0:
        dose_fast = dose_slow = 0.5 * dose_total

    A_fast = k_fast * dose_fast
    A_slow = k_slow * dose_slow

    A_prot = params["alpha_prot_uU_mL_per_g"] * protein
    k_prot = params["k_prot_min1"]

    return MealAppearance(A_fast, k_fast, A_slow, k_slow, A_prot, k_prot, dose_fast, dose_slow, dose_total)


def _exp_decay(amplitude: float, rate: float, t: float) -> float:
    if amplitude == 0.0:
        return 0.0
    if rate <= 0.0:
        return amplitude
    return amplitude * math.exp(-rate * t)


def simulate(name: str, params: Dict[str, float], kinetics: MealAppearance) -> SimulationResult:
    if bindings.lib is None:
        raise RuntimeError("C backend not loaded")
    dt = params["dt_min"]
    t_end = params["t_end_min"]
    steps = int(round(t_end / dt))
    times = np.linspace(0.0, t_end, steps + 1)

    glucose = np.zeros_like(times)
    insulin = np.zeros_like(times)
    d_fast = np.zeros_like(times)
    d_slow = np.zeros_like(times)
    u_series = np.zeros_like(times)

    G = c_double(params["Gb_mg_dL"])
    X = c_double(0.0)
    I = c_double(params["Ib_uU_mL"])

    Gb = c_double(params["Gb_mg_dL"])
    Ib = c_double(params["Ib_uU_mL"])
    dt_c = c_double(dt)

    for idx, t in enumerate(times):
        glucose[idx] = G.value
        insulin[idx] = I.value
        d_fast[idx] = _exp_decay(kinetics.A_fast, kinetics.k_fast, t)
        d_slow[idx] = _exp_decay(kinetics.A_slow, kinetics.k_slow, t)
        u_series[idx] = _exp_decay(kinetics.A_prot, kinetics.k_prot, t)
        if idx == len(times) - 1:
            break
        D0 = d_fast[idx] + d_slow[idx]
        u0 = u_series[idx]
        t_half = t + 0.5 * dt_c.value
        t_full = t + dt_c.value
        D_half = _exp_decay(kinetics.A_fast, kinetics.k_fast, t_half) + _exp_decay(kinetics.A_slow, kinetics.k_slow, t_half)
        u_half = _exp_decay(kinetics.A_prot, kinetics.k_prot, t_half)
        D_full = _exp_decay(kinetics.A_fast, kinetics.k_fast, t_full) + _exp_decay(kinetics.A_slow, kinetics.k_slow, t_full)
        u_full = _exp_decay(kinetics.A_prot, kinetics.k_prot, t_full)
        D_eff = (D0 + 4.0 * D_half + D_full) / 6.0
        u_eff = (u0 + 4.0 * u_half + u_full) / 6.0
        bindings.lib.step(byref(G), byref(X), byref(I), Gb.value, Ib.value, dt_c.value, D_eff, u_eff)

    return SimulationResult(name=name, times=times, glucose=glucose, insulin=insulin, d_fast=d_fast, d_slow=d_slow, u=u_series)


def _clip_series(times: np.ndarray, values: np.ndarray, t_limit: float) -> Tuple[np.ndarray, np.ndarray]:
    if t_limit >= times[-1]:
        return times, values
    idx = int(np.searchsorted(times, t_limit))
    if math.isclose(times[idx], t_limit, rel_tol=1e-9, abs_tol=1e-9):
        return times[: idx + 1], values[: idx + 1]
    t0 = times[idx - 1]
    t1 = times[idx]
    v0 = values[idx - 1]
    v1 = values[idx]
    frac = (t_limit - t0) / (t1 - t0)
    v_limit = v0 + frac * (v1 - v0)
    clipped_times = np.concatenate([times[:idx], np.array([t_limit])])
    clipped_values = np.concatenate([values[:idx], np.array([v_limit])])
    return clipped_times, clipped_values


def incremental_auc(times: np.ndarray, deltas: np.ndarray, window: float) -> float:
    clipped_times, clipped_deltas = _clip_series(times, deltas, window)
    positive = np.maximum(clipped_deltas, 0.0)
    return float(np.trapezoid(positive, clipped_times))


def time_to_baseline(times: np.ndarray, glucose: np.ndarray, baseline: float, tolerance: float = 5.0) -> float:
    deltas = np.abs(glucose - baseline)
    peak_idx = int(np.argmax(glucose))
    for idx in range(peak_idx, len(times)):
        if deltas[idx] <= tolerance and np.all(deltas[idx:] <= tolerance):
            return float(times[idx])
    return math.nan


def compute_metrics(result: SimulationResult, params: Dict[str, float]) -> Dict[str, float]:
    Gb = params["Gb_mg_dL"]
    delta = result.glucose - Gb
    peak_idx = int(delta.argmax())
    peak_delta = float(delta[peak_idx])
    t_peak = float(result.times[peak_idx])
    if 0 < peak_idx < len(delta) - 1:
        dt_left = result.times[peak_idx] - result.times[peak_idx - 1]
        dt_right = result.times[peak_idx + 1] - result.times[peak_idx]
        if dt_left > 0 and dt_right > 0 and math.isclose(dt_left, dt_right, rel_tol=1e-9, abs_tol=1e-9):
            y_prev = delta[peak_idx - 1]
            y_curr = delta[peak_idx]
            y_next = delta[peak_idx + 1]
            denom = y_prev - 2.0 * y_curr + y_next
            if abs(denom) > 1e-12:
                offset = 0.5 * (y_prev - y_next) / denom
                dt_avg = 0.5 * (dt_left + dt_right)
                t_peak = float(result.times[peak_idx] + offset * dt_avg)
                peak_delta = float(y_curr - 0.25 * (y_prev - y_next) * offset)
    iauc_120 = incremental_auc(result.times, delta, 120.0)
    iauc_240 = incremental_auc(result.times, delta, 240.0)
    t_baseline = time_to_baseline(result.times, result.glucose, Gb)
    insulin_peak_idx = int(np.argmax(result.insulin))
    insulin_peak_time = float(result.times[insulin_peak_idx])
    return {
        "peak_delta_G": peak_delta,
        "t_peak_min": t_peak,
        "iAUC_0_120": iauc_120,
        "iAUC_0_240": iauc_240,
        "time_to_baseline": t_baseline,
        "insulin_peak_time": insulin_peak_time,
    }


def plot_result(result: SimulationResult, params: Dict[str, float], zoom_end: float, full_end: float, out_root: Path, dpi: int) -> None:
    figures_root = out_root / "figures"
    zoom_dir = figures_root / "zoom"
    full_dir = figures_root / "full"
    zoom_dir.mkdir(parents=True, exist_ok=True)
    full_dir.mkdir(parents=True, exist_ok=True)

    def _plot(range_end: float, path: Path) -> None:
        mask = result.times <= range_end
        times = result.times[mask]
        glucose = result.glucose[mask]
        insulin = result.insulin[mask]
        D_total = (result.d_fast + result.d_slow)[mask]
        u_series = result.u[mask]

        fig, (ax_g, ax_i) = plt.subplots(2, 1, sharex=True, figsize=(9, 6))
        ax_g.plot(times, glucose, label="Glucose (mg/dL)", color="#1f77b4")
        ax_g.axhline(params["Gb_mg_dL"], linestyle="--", color="#444444", linewidth=1.0, label="Gb")
        ax_g.set_ylabel("Glucose (mg/dL)")
        ax_g.legend(loc="upper right")

        ax_i.plot(times, insulin, label="Insulin (μU/mL)", color="#d62728")
        ax_i.axhline(params["Ib_uU_mL"], linestyle="--", color="#555555", linewidth=1.0, label="Ib")
        ax_i.set_ylabel("Insulin (μU/mL)")
        ax_i.set_xlabel("Time (min)")
        ax_i.legend(loc="upper right")

        ax_i2 = ax_i.twinx()
        ax_i2.plot(times, D_total, color="#2ca02c", alpha=0.6, label="D(t)")
        ax_i2.plot(times, u_series, color="#9467bd", alpha=0.6, label="u(t)")
        ax_i2.set_ylabel("Appearance / stimulus")
        lines, labels = ax_i.get_legend_handles_labels()
        lines2, labels2 = ax_i2.get_legend_handles_labels()
        ax_i.legend(lines + lines2, labels + labels2, loc="upper left")

        fig.suptitle(result.name)
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi)
        plt.close(fig)

    _plot(zoom_end, zoom_dir / f"{result.name}.png")
    _plot(full_end, full_dir / f"{result.name}.png")


def write_summary(rows: Sequence[Tuple[str, Dict[str, float]]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "dessert",
        "peak_delta_G_mg_dL",
        "t_peak_min",
        "iAUC_0_120_mg_dL_min",
        "iAUC_0_240_mg_dL_min",
        "time_to_baseline_min",
        "insulin_peak_time_min",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for name, metrics in rows:
            writer.writerow(
                {
                    "dessert": name,
                    "peak_delta_G_mg_dL": f"{metrics['peak_delta_G']:.4f}",
                    "t_peak_min": f"{metrics['t_peak_min']:.2f}",
                    "iAUC_0_120_mg_dL_min": f"{metrics['iAUC_0_120']:.2f}",
                    "iAUC_0_240_mg_dL_min": f"{metrics['iAUC_0_240']:.2f}",
                    "time_to_baseline_min": "" if math.isnan(metrics["time_to_baseline"]) else f"{metrics['time_to_baseline']:.2f}",
                    "insulin_peak_time_min": f"{metrics['insulin_peak_time']:.2f}",
                }
            )


def write_json(path: Path, payload: Dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def collect_env_metadata() -> Dict[str, object]:
    import platform

    build = platform.python_build()
    return {
        "python": platform.python_version(),
        "python_build": " ".join(build),
        "numpy": np.__version__,
        "os": platform.platform(),
    }


def run_sensitivity(configs: Sequence[Tuple[str, Dict[str, float]]], params: Dict[str, float], calibrate: bool, out_dir: Path) -> None:
    rows: List[Dict[str, object]] = []
    perturbations = {"fat_g": "fat", "fiber_g": "fiber", "protein_g": "protein"}
    for name, dessert in configs:
        for key, label in perturbations.items():
            base_value = float(dessert.get(key, 0.0))
            for factor in (0.8, 1.2):
                modified = dict(dessert)
                modified[key] = max(0.0, base_value * factor)
                kinetics = compute_meal_appearance(modified, params, calibrate)
                result = simulate(name, params, kinetics)
                metrics = compute_metrics(result, params)
                rows.append(
                    {
                        "dessert": name,
                        "nutrient": label,
                        "factor": factor,
                        "peak_delta_G_mg_dL": metrics["peak_delta_G"],
                        "t_peak_min": metrics["t_peak_min"],
                        "iAUC_0_240_mg_dL_min": metrics["iAUC_0_240"],
                    }
                )
    path = out_dir / "tables" / "sensitivity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["dessert", "nutrient", "factor", "peak_delta_G_mg_dL", "t_peak_min", "iAUC_0_240_mg_dL_min"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_frozen_configs() -> List[Tuple[str, Dict[str, float]]]:
    configs: List[Tuple[str, Dict[str, float]]] = []
    for yaml_path in sorted(FROZEN_DIR.glob("*.yaml")):
        with yaml_path.open("r", encoding="utf-8") as handle:
            payload = yaml.safe_load(handle) or {}
        if not isinstance(payload, dict):
            raise ValueError(f"Frozen config {yaml_path} must be a mapping")
        name = str(payload.get("name", yaml_path.stem))
        configs.append((name, payload))
    return configs


def clean_outputs(base: Path) -> None:
    for rel in ("build", "figures", "tables"):
        target = base / rel
        if target.exists():
            for path in sorted(target.glob("**/*"), reverse=True):
                if path.is_file() or path.is_symlink():
                    path.unlink(missing_ok=True)
            for path in sorted(target.glob("**/*"), reverse=True):
                if path.is_dir():
                    try:
                        path.rmdir()
                    except OSError:
                        pass
            if target.exists():
                try:
                    target.rmdir()
                except OSError:
                    pass


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Dessert glycemia simulator")
    parser.add_argument("--reproduce", action="store_true", help="Clean, build, and run frozen configs")
    parser.add_argument("--calibrate", action="store_true", help="Equalize fast/slow appearance amplitudes")
    parser.add_argument("--sensitivity", action="store_true", help="Run ±20% nutrient perturbations")
    parser.add_argument("--window", choices=["0-120", "0-240", "0-1440"], default="0-240", help="Zoom window for plots")
    parser.add_argument("--no-plots", action="store_true", help="Skip figure generation")
    parser.add_argument("--summary-only", action="store_true", help="Only write summary table")
    parser.add_argument("--dpi", type=int, default=150, help="Figure DPI")
    parser.add_argument("--outdir", type=Path, default=None, help="Output directory (defaults to repository root)")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir).resolve() if args.outdir else ROOT

    if args.reproduce:
        clean_outputs(outdir)

    params = load_params_yaml(DEFAULT_PARAMS_PATH)
    set_params_from_dict(params)

    configs = load_frozen_configs()

    zoom_choice = args.window
    zoom_end = float(zoom_choice.split("-")[1])
    full_end = 1440.0

    generate_plots = not (args.no_plots or args.summary_only)

    summary_rows: List[Tuple[str, Dict[str, float]]] = []

    for name, dessert in configs:
        kinetics = compute_meal_appearance(dessert, params, args.calibrate)
        result = simulate(name, params, kinetics)
        metrics = compute_metrics(result, params)
        summary_rows.append((name, metrics))
        if generate_plots:
            plot_result(result, params, zoom_end=zoom_end, full_end=full_end, out_root=outdir, dpi=args.dpi)

    summary_path = outdir / "tables" / "summary.csv"
    write_summary(summary_rows, summary_path)

    env_meta = collect_env_metadata()
    write_json(outdir / "build" / "env.json", env_meta)
    write_json(outdir / "build" / "build.json", get_build_metadata())

    if args.sensitivity:
        run_sensitivity(configs, params, args.calibrate, outdir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
