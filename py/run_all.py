"""Run all dessert simulations and regenerate figures/tables."""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Tuple

import yaml
import numpy as np

from .desserts import desserts
from .interface import simulate_c
from .metrics import auc, peak_and_tpeak
from .plotting import plot_overlay, plot_single

_ROOT = Path(__file__).resolve().parents[1]
_PARAMS_PATH = _ROOT / "configs" / "params.yaml"
_FIG_DIR = _ROOT / "figures"
_TABLE_DIR = _ROOT / "tables"


def _load_params() -> Dict[str, float]:
    with _PARAMS_PATH.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _write_summary(rows: Tuple[Dict[str, float], ...]) -> None:
    _TABLE_DIR.mkdir(parents=True, exist_ok=True)
    path = _TABLE_DIR / "summary.csv"
    fieldnames = [
        "dessert",
        "peak_G",
        "t_peak_G",
        "AUC_G_aboveGb_0_120",
        "peak_I",
        "t_peak_I",
        "AUC_I_0_120",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    params = _load_params()

    glucose_curves: Dict[str, np.ndarray] = {}
    insulin_curves: Dict[str, np.ndarray] = {}
    summary_rows = []
    time_ref = None

    for name, (A, k) in desserts.items():
        t, G, I = simulate_c(params, A, k)
        if time_ref is None:
            time_ref = t

        glucose_curves[name] = G
        insulin_curves[name] = I

        plot_single(t, G, f"Glucose response: {name}", "Glucose (mg/dL)", _FIG_DIR / f"{name}_glucose.png")
        plot_single(t, I, f"Insulin response: {name}", "Insulin (µU/mL)", _FIG_DIR / f"{name}_insulin.png")

        peak_G, t_peak_G = peak_and_tpeak(G, t)
        peak_I, t_peak_I = peak_and_tpeak(I, t)
        auc_G = auc(G, t, baseline=float(params["Gb"]), tmax=120.0)
        auc_I = auc(I, t, tmax=120.0)

        summary_rows.append(
            {
                "dessert": name,
                "peak_G": round(peak_G, 6),
                "t_peak_G": round(t_peak_G, 6),
                "AUC_G_aboveGb_0_120": round(auc_G, 6),
                "peak_I": round(peak_I, 6),
                "t_peak_I": round(t_peak_I, 6),
                "AUC_I_0_120": round(auc_I, 6),
            }
        )

    if time_ref is None:
        raise RuntimeError("No desserts defined.")

    plot_overlay(time_ref, glucose_curves, "Glucose (mg/dL)", _FIG_DIR / "glucose_overlay.png")
    plot_overlay(time_ref, insulin_curves, "Insulin (µU/mL)", _FIG_DIR / "insulin_overlay.png")

    _write_summary(tuple(summary_rows))

    print("Done. See figures/ and tables/.")


if __name__ == "__main__":
    main()
