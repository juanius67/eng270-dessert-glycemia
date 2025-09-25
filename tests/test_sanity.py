from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path

import numpy as np
import pytest
import yaml
from scipy.integrate import odeint

ROOT = Path(__file__).resolve().parents[1]
C_DIR = ROOT / "C"
PARAMS_PATH = ROOT / "configs" / "params.yaml"

_spec = importlib.util.spec_from_file_location("eng270_interface", ROOT / "py" / "interface.py")
_interface = importlib.util.module_from_spec(_spec)
assert _spec is not None and _spec.loader is not None
_spec.loader.exec_module(_interface)
simulate_c = _interface.simulate_c


@pytest.fixture(scope="session", autouse=True)
def build_library() -> None:
    subprocess.run(["make"], cwd=C_DIR, check=True)


def _load_params() -> dict[str, float]:
    with PARAMS_PATH.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def test_baseline_equilibrium() -> None:
    params = _load_params()
    t, G, I = simulate_c(params, 0.0, 0.1)
    assert np.max(np.abs(G - params["Gb"])) < 1e-3
    assert np.max(np.abs(I - params["Ib"])) < 1e-3


def test_no_insulin_matches_odeint() -> None:
    params = _load_params()
    params["p3"] = 0.0
    params["p4"] = 0.0

    A_test = 180.0
    k_test = 0.12

    t, G_c, I_c = simulate_c(params, A_test, k_test)

    def rhs(y, t_val):
        G_val, X_val, I_val = y
        Dt = A_test * np.exp(-k_test * t_val)
        dG = -(params["p1"] + X_val) * G_val + params["p1"] * params["Gb"] + Dt
        dX = -params["p2"] * X_val
        dI = -params["p6"] * (I_val - params["Ib"])
        return [dG, dX, dI]

    y0 = [params["Gb"], 0.0, params["Ib"]]
    sol = odeint(rhs, y0, t)

    assert np.max(np.abs(sol[:, 0] - G_c)) < 1e-5
    assert np.max(np.abs(sol[:, 2] - I_c)) < 1e-5


def test_dt_refinement_stability() -> None:
    params_coarse = _load_params()
    params_coarse["dt"] = 0.02
    params_coarse["t_end"] = 60.0

    params_fine = params_coarse.copy()
    params_fine["dt"] = params_coarse["dt"] / 2.0

    A_test = 200.0
    k_test = 0.1

    t_coarse, G_coarse, _ = simulate_c(params_coarse, A_test, k_test)
    t_fine, G_fine, _ = simulate_c(params_fine, A_test, k_test)

    np.testing.assert_allclose(t_coarse, t_fine[::2], rtol=0, atol=1e-12)
    max_dev = np.max(np.abs(G_coarse - G_fine[::2]))
    assert max_dev < 5e-4
