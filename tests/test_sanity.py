"""Sanity checks for dessert glycemia simulations."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable

import pytest

import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from run import (
    CONFIG_PATH,
    DESSERTS_CONFIG_PATH,
    DoseProfile,
    Simulator,
    compute_metrics,
    load_config,
    load_dessert_catalog,
)


@pytest.fixture(scope="module")
def params() -> Dict[str, float]:
    return load_config(CONFIG_PATH)


@pytest.fixture(scope="module")
def dessert_profiles() -> Dict[str, Dict[str, object]]:
    return load_dessert_catalog(DESSERTS_CONFIG_PATH)


@pytest.fixture()
def simulator() -> Simulator:
    return Simulator()


def _relative_diff(a: float, b: float) -> float:
    denom = max(1.0, abs(a))
    return abs(a - b) / denom


def test_steady_state(simulator: Simulator, params: Dict[str, float]) -> None:
    dt = params["dt"]
    t_end = params["t_end"]
    profile = DoseProfile(0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    result = simulator.run(params, dt, t_end, profile)

    gb = params["Gb_mg_dL"]
    ib = params["Ib_uU_mL"]

    max_glucose_dev = max(abs(value - gb) for value in result.glucose)
    max_insulin_dev = max(abs(value - ib) for value in result.insulin)

    assert max(max_glucose_dev, max_insulin_dev) < 1e-3


def test_dt_halving(params: Dict[str, float], dessert_profiles: Dict[str, Dict[str, object]]) -> None:
    simulator = Simulator()
    dessert_name = next(iter(dessert_profiles.keys()))
    profile = dessert_profiles[dessert_name]["profile"]  # type: ignore[index]

    dt = params["dt"]
    t_end = params["t_end"]

    base = simulator.run(params, dt, t_end, profile)
    refined = simulator.run(params, dt / 2.0, t_end, profile)

    metrics_base = compute_metrics(base, params)
    metrics_refined = compute_metrics(refined, params)

    keys: Iterable[str] = (
        "peak_G",
        "iAUC_0_120",
        "iAUC_0_240",
        "peak_I",
        "AUCI_0_240",
    )
    for key in keys:
        assert _relative_diff(metrics_base[key], metrics_refined[key]) < 0.01


@pytest.mark.skipif(Simulator().label == "python", reason="C backend missing")
def test_c_vs_python(params: Dict[str, float], dessert_profiles: Dict[str, Dict[str, object]]) -> None:
    simulator = Simulator()
    dessert_name = next(iter(dessert_profiles.keys()))
    profile = dessert_profiles[dessert_name]["profile"]  # type: ignore[index]

    dt = params["dt"]
    t_end = params["t_end"]

    if simulator._has_extended:  # type: ignore[attr-defined]
        result_c = simulator._run_c_extended(params, dt, t_end, profile)  # type: ignore[attr-defined]
    else:
        result_c = simulator._run_c_legacy(params, dt, t_end, profile)  # type: ignore[attr-defined]
    result_py = simulator._run_python(params, dt, t_end, profile)  # type: ignore[attr-defined]

    diffs = [abs(g_c - g_py) for g_c, g_py in zip(result_c.glucose, result_py.glucose)]
    assert max(diffs) < 1e-2
