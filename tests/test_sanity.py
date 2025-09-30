"""Sanity checks for dessert glycemia simulations."""

from __future__ import annotations

from typing import Dict

import pytest

from run import CONFIG_PATH, DESSERTS, Simulator, load_config


@pytest.fixture(scope="module")
def params() -> Dict[str, float]:
    return load_config(CONFIG_PATH)


@pytest.fixture()
def simulator() -> Simulator:
    return Simulator()


def test_baseline_remains_at_basals(simulator: Simulator, params: Dict[str, float]) -> None:
    dt = params["dt"]
    t_end = params["t_end"]
    result = simulator.run(params, dt, t_end, 0.0, 0.0)

    gb = params["Gb"]
    ib = params["Ib"]

    max_glucose_dev = max(abs(value - gb) for value in result.glucose)
    max_insulin_dev = max(abs(value - ib) for value in result.insulin)

    assert max(max_glucose_dev, max_insulin_dev) < 1e-3


def test_dt_halving_consistency(params: Dict[str, float]) -> None:
    simulator = Simulator()
    dessert_name = "chocotorta"
    specs = DESSERTS[dessert_name]

    dt = params["dt"]
    t_end = params["t_end"]
    A = specs["k"] * specs["dose_mgdL"]

    base = simulator.run(params, dt, t_end, A, specs["k"])
    refined = simulator.run(params, dt / 2.0, t_end, A, specs["k"])

    aligned_refined_glucose = refined.glucose[::2][: len(base.glucose)]

    max_glucose_diff = max(
        abs(g_base - g_ref)
        for g_base, g_ref in zip(base.glucose, aligned_refined_glucose)
    )

    assert max_glucose_diff < 1e-2


@pytest.mark.skipif(Simulator().label != "c", reason="C backend missing")
def test_c_backend_matches_python(params: Dict[str, float]) -> None:
    simulator = Simulator()
    dessert_name = "brigadeiro"
    specs = DESSERTS[dessert_name]

    dt = params["dt"]
    t_end = params["t_end"]
    A = specs["k"] * specs["dose_mgdL"]

    result_c = simulator._run_c(params, dt, t_end, A, specs["k"])
    result_py = simulator._run_python(params, dt, t_end, A, specs["k"])

    max_glucose_diff = max(
        abs(g_c - g_py) for g_c, g_py in zip(result_c.glucose, result_py.glucose)
    )

    assert max_glucose_diff < 1e-2
