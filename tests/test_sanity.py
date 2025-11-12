import yaml
from pathlib import Path

from run import simulate_no_meal, simulate_peak_deltaG


def test_steady_state_no_meal():
    params = yaml.safe_load(Path("configs/params.yaml").read_text())
    G, I = simulate_no_meal(params, t_end_min=120.0, dt_min=0.5)
    assert abs(G[-1] - params["Gb_mg_dL"]) < 0.5
    assert abs(I[-1] - params["Ib_uU_mL"]) < 0.5


def test_dt_halving_peak_deltaG():
    params = yaml.safe_load(Path("configs/params.yaml").read_text())
    fp = sorted(Path("configs/frozen").glob("*.yaml"))[0]
    peak1 = simulate_peak_deltaG(params, fp, dt_min=0.5)
    peak2 = simulate_peak_deltaG(params, fp, dt_min=0.25)
    rel = abs(peak2 - peak1) / max(1e-6, abs(peak1))
    assert rel < 0.01
