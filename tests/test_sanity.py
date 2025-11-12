"""Sanity checks for the dessert glycemia CLI."""
from __future__ import annotations

import math
from pathlib import Path

import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from run import compute_meal_appearance, compute_metrics, load_params_yaml, simulate
from src.bindings import set_params_from_dict


def test_steady_state_no_meal() -> None:
    params = load_params_yaml(ROOT / "configs" / "params.yaml")
    params = dict(params)
    params["t_end_min"] = 120.0
    params["dt_min"] = 1.0
    set_params_from_dict(params)

    fasting_kinetics = compute_meal_appearance({}, params, calibrate=False)
    result = simulate("fasting", params, fasting_kinetics)

    G_end = result.glucose[-1]
    I_end = result.insulin[-1]

    assert abs(G_end - params["Gb_mg_dL"]) < 0.5
    assert abs(I_end - params["Ib_uU_mL"]) < 0.5


def test_dt_halving_peak_delta_stability() -> None:
    params = load_params_yaml(ROOT / "configs" / "params.yaml")
    dessert = {
        "carbs_g": 30.0,
        "sugars_g": 18.0,
        "fiber_g": 2.0,
        "fat_g": 5.0,
        "protein_g": 3.0,
    }

    params_coarse = dict(params)
    params_coarse["dt_min"] = params["dt_min"]
    params_coarse["t_end_min"] = 300.0
    set_params_from_dict(params_coarse)
    kinetics_coarse = compute_meal_appearance(dessert, params_coarse, calibrate=False)
    result_coarse = simulate("test", params_coarse, kinetics_coarse)
    metrics_coarse = compute_metrics(result_coarse, params_coarse)

    params_fine = dict(params)
    params_fine["dt_min"] = params["dt_min"] / 2.0
    params_fine["t_end_min"] = 300.0
    set_params_from_dict(params_fine)
    kinetics_fine = compute_meal_appearance(dessert, params_fine, calibrate=False)
    result_fine = simulate("test", params_fine, kinetics_fine)
    metrics_fine = compute_metrics(result_fine, params_fine)

    peak_coarse = metrics_coarse["peak_delta_G"]
    peak_fine = metrics_fine["peak_delta_G"]
    rel_diff = abs(peak_coarse - peak_fine) / max(1.0, abs(peak_fine))
    assert rel_diff < 0.01

    assert not math.isnan(peak_coarse)
    assert not math.isnan(peak_fine)
