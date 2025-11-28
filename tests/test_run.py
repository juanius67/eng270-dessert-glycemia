from __future__ import annotations
import math
import numpy as np
import pytest
import yaml
from pathlib import Path
from unittest.mock import MagicMock, patch, mock_open

from run import (
    MealAppearance,
    SimulationResult,
    load_params_yaml,
    map_meal_from_label,
    _equalize_amplitudes,
    build_meal_appearance,
    incremental_auc,
    time_to_baseline,
    compute_metrics,
    _exp_decay,
    simulate_no_meal,
    simulate_peak_deltaG,
    _meal_from_mapping,
    _clip,
    _clip_series,
    _prepare_params
)

@pytest.fixture
def sample_params():
    return {
        "Gb_mg_dL": 90.0,
        "Ib_uU_mL": 10.0,
        "S_G_min1": 0.02,
        "p2_min1": 0.01,
        "p3_min1_per_uU_mL": 1e-5,
        "n_min1": 0.1,
        "Vd_dL": 150.0,
        "f_hep": 0.2,
        "f_app0": 0.5,
        "beta_fiber_per10g": 0.1,
        "beta_fat_per10g": 0.1,
        "k_fast0_min1": 0.05,
        "k_slow0_min1": 0.01,
        "alpha_prot_uU_mL_per_g": 0.1,
        "k_prot_min1": 0.02,
        "dt_min": 1.0,
        "t_end_min": 120.0,
    }

def test_load_params_yaml(tmp_path, sample_params):
    f = tmp_path / "params.yaml"
    with f.open("w") as fp:
        yaml.dump(sample_params, fp)

    loaded = load_params_yaml(f)
    assert loaded == sample_params

    # Test missing key
    incomplete = sample_params.copy()
    del incomplete["Gb_mg_dL"]
    with f.open("w") as fp:
        yaml.dump(incomplete, fp)

    with pytest.raises(ValueError, match="Missing parameters"):
        load_params_yaml(f)

    # Test invalid format
    with f.open("w") as fp:
        fp.write("not a dict")
    with pytest.raises(ValueError, match="must contain a mapping"):
        load_params_yaml(f)

def test_clip():
    assert _clip(0.5, 0.0, 1.0) == 0.5
    assert _clip(-1.0, 0.0, 1.0) == 0.0
    assert _clip(2.0, 0.0, 1.0) == 1.0

def test_map_meal_from_label(sample_params):
    dessert = {
        "carbs_g": 50.0,
        "sugars_g": 25.0,
        "fiber_g": 5.0,
        "fat_g": 10.0,
        "protein_g": 5.0,
    }

    mapping = map_meal_from_label(dessert, sample_params)

    assert "A_fast" in mapping
    assert "k_fast" in mapping
    assert "A_slow" in mapping
    assert "k_slow" in mapping
    assert "A_prot" in mapping
    assert "k_prot" in mapping

    assert mapping["k_prot"] == sample_params["k_prot_min1"]
    assert mapping["A_prot"] == sample_params["alpha_prot_uU_mL_per_g"] * dessert["protein_g"]

def test_equalize_amplitudes():
    meal = MealAppearance(
        A_fast=10.0, k_fast=0.1,
        A_slow=2.0, k_slow=0.01,
        A_prot=0.0, k_prot=0.01,
        dose_fast=100.0, dose_slow=200.0, dose_total=300.0
    )

    equalized = _equalize_amplitudes(meal)
    assert equalized.dose_total == meal.dose_total
    assert math.isclose(equalized.A_fast, equalized.A_slow)
    assert equalized.k_fast == meal.k_fast
    assert equalized.k_slow == meal.k_slow

    # Check invalid case
    invalid_meal = MealAppearance(
        A_fast=10.0, k_fast=-0.1,
        A_slow=2.0, k_slow=0.01,
        A_prot=0.0, k_prot=0.01,
        dose_fast=100.0, dose_slow=200.0, dose_total=300.0
    )
    assert _equalize_amplitudes(invalid_meal) == invalid_meal

def test_incremental_auc():
    times = np.array([0.0, 10.0, 20.0])
    deltas = np.array([0.0, 10.0, 5.0]) # Triangle-ish
    # Trapezoid rule:
    # 0-10: 0.5 * 10 * 10 = 50
    # 10-20: 0.5 * 10 * (10+5) = 75
    # Total = 125
    assert math.isclose(incremental_auc(times, deltas, 20.0), 125.0)

    # With negative values
    deltas_neg = np.array([0.0, 10.0, -5.0])
    # 0-10: 50
    # 10-20: crosses zero at some point.
    # But incremental_auc clamps negative to 0.
    # _clip_series will likely interpolate at 20.0 if window is 20.0.
    # Wait, positive = np.maximum(clipped_deltas, 0.0) handles negativity.
    # At 20.0 it is -5.0, so 0.0.

    val = incremental_auc(times, deltas_neg, 20.0)
    assert val > 0

def test_time_to_baseline():
    times = np.array([0.0, 10.0, 20.0, 30.0, 40.0])
    glucose = np.array([90.0, 120.0, 100.0, 92.0, 90.0])
    baseline = 90.0
    tolerance = 5.0

    # Peak is at 10.0 (120.0)
    # 20.0: |100-90| = 10 > 5
    # 30.0: |92-90| = 2 <= 5
    # 40.0: |90-90| = 0 <= 5
    # Should return 30.0

    assert time_to_baseline(times, glucose, baseline, tolerance) == 30.0

    # Never returns
    glucose_high = np.array([90.0, 120.0, 110.0, 110.0, 110.0])
    assert math.isnan(time_to_baseline(times, glucose_high, baseline, tolerance))

def test_compute_metrics(sample_params):
    times = np.linspace(0, 120, 13) # 0, 10, 20...
    glucose = np.full_like(times, 90.0)
    # Make a peak
    glucose[2] = 140.0 # at t=20
    glucose[1] = 115.0
    glucose[3] = 115.0

    insulin = np.full_like(times, 10.0)
    insulin[3] = 50.0 # Peak at t=30

    result = SimulationResult(
        name="test",
        times=times,
        glucose=glucose,
        insulin=insulin,
        d_fast=np.zeros_like(times),
        d_slow=np.zeros_like(times),
        u=np.zeros_like(times)
    )

    metrics = compute_metrics(result, sample_params)

    assert math.isclose(metrics["peak_delta_G"], 50.0, abs_tol=1.0) # 140 - 90
    assert math.isclose(metrics["t_peak_min"], 20.0, abs_tol=1.0)
    assert metrics["insulin_peak_time"] == 30.0
    assert metrics["iAUC_0_120"] > 0

def test_exp_decay():
    assert _exp_decay(10.0, 0.1, 0.0) == 10.0
    assert _exp_decay(0.0, 0.1, 10.0) == 0.0
    assert _exp_decay(10.0, -0.1, 10.0) == 10.0
    val = _exp_decay(10.0, 0.1, 10.0)
    assert math.isclose(val, 10.0 * math.exp(-1.0))

@patch("run.simulate")
def test_simulate_no_meal(mock_simulate, sample_params):
    mock_result = SimulationResult(
        name="fasting",
        times=np.array([0, 1]),
        glucose=np.array([90, 90]),
        insulin=np.array([10, 10]),
        d_fast=np.array([0, 0]),
        d_slow=np.array([0, 0]),
        u=np.array([0, 0])
    )
    mock_simulate.return_value = mock_result

    g, i = simulate_no_meal(sample_params, 120.0, 1.0)
    assert len(g) == 2
    assert len(i) == 2

    mock_simulate.assert_called_once()
    args = mock_simulate.call_args
    assert args[0][0] == "fasting"
    kinetics = args[0][2]
    assert kinetics.dose_total == 0.0

@patch("run.simulate")
def test_simulate_peak_deltaG(mock_simulate, sample_params, tmp_path):
    mock_result = SimulationResult(
        name="test",
        times=np.array([0, 1]),
        glucose=np.array([90, 140]), # Peak 140 -> delta 50
        insulin=np.array([10, 10]),
        d_fast=np.array([0, 0]),
        d_slow=np.array([0, 0]),
        u=np.array([0, 0])
    )
    mock_simulate.return_value = mock_result

    dessert_file = tmp_path / "dessert.yaml"
    with dessert_file.open("w") as f:
        yaml.dump({"name": "Cake", "carbs_g": 50}, f)

    peak = simulate_peak_deltaG(sample_params, dessert_file, 1.0)
    assert peak == 50.0

def test_prepare_params(sample_params):
    updated = _prepare_params(sample_params, dt_min=0.5, t_end_min=240.0)
    assert updated["dt_min"] == 0.5
    assert updated["t_end_min"] == 240.0
    assert updated["Gb_mg_dL"] == sample_params["Gb_mg_dL"]

def test_clip_series():
    times = np.array([0.0, 10.0, 20.0])
    values = np.array([0.0, 100.0, 200.0])

    # Clip at 15.0
    # idx=2 (20.0)
    # t0=10, t1=20
    # frac = 0.5
    # v = 100 + 0.5 * 100 = 150

    ct, cv = _clip_series(times, values, 15.0)
    assert ct[-1] == 15.0
    assert cv[-1] == 150.0
    assert len(ct) == 3

    # Clip exactly at point
    ct, cv = _clip_series(times, values, 10.0)
    assert ct[-1] == 10.0
    assert cv[-1] == 100.0
    assert len(ct) == 2

    # Clip beyond end
    ct, cv = _clip_series(times, values, 30.0)
    assert np.array_equal(ct, times)
    assert np.array_equal(cv, values)
