#!/usr/bin/env python3
"""Compact entrypoint for dessert-driven glycemia experiments."""
from __future__ import annotations

import argparse
import csv
import ctypes
import hashlib
import json
import math
import os
import platform
import shutil
import struct
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

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
BUILD_DIR = ROOT / "build"
LIB_WIN = ROOT / "C" / "model.dll"
LIB_LIN = ROOT / "C" / "libmodel.so"
LIB_MAC = ROOT / "C" / "libmodel.dylib"
BUILD_META = BUILD_DIR / "build.json"
ENV_JSON = BUILD_DIR / "env.json"
SANITY_JSON = BUILD_DIR / "sanity.json"


def sha256(path: os.PathLike[str] | str) -> str:
    digest = hashlib.sha256()
    file_path = Path(path)
    with file_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def probe_arch() -> Tuple[str, str, int]:
    system = platform.system()
    machine = platform.machine()
    bits = struct.calcsize("P") * 8
    return system, machine, bits


def write_json(path: Path, obj: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2, sort_keys=True)


def read_json(path: Path) -> Dict[str, Any]:
    try:
        with Path(path).open("r", encoding="utf-8") as handle:
            return json.load(handle)
    except FileNotFoundError:
        return {}
    except json.JSONDecodeError:
        return {}


def compute_source_hash() -> str:
    model_c = C_DIR / "model.c"
    model_h = C_DIR / "model.h"
    return sha256(model_c) + sha256(model_h)


def get_backend_spec() -> Tuple[Path, List[List[str]]]:
    system = platform.system().lower()
    if system.startswith("win"):
        return (
            LIB_WIN,
            [
                ["cl", "/nologo", "/LD", "C\\model.c", "/Fe:C\\model.dll"],
                ["gcc", "-O3", "-shared", "-fPIC", "C/model.c", "-o", "C/model.dll"],
            ],
        )
    if system == "darwin":
        return (
            LIB_MAC,
            [["clang", "-O3", "-dynamiclib", "C/model.c", "-o", "C/libmodel.dylib"]],
        )
    return (
        LIB_LIN,
        [["gcc", "-O3", "-shared", "-fPIC", "C/model.c", "-o", "C/libmodel.so"]],
    )


def build_backend() -> Path:
    lib_path, commands = get_backend_spec()
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    source_hash = compute_source_hash()
    system, machine, bits = probe_arch()
    meta = read_json(BUILD_META)
    if (
        lib_path.exists()
        and meta.get("hash") == source_hash
        and meta.get("system") == system
        and meta.get("machine") == machine
        and meta.get("bits") == bits
        and meta.get("compiler") in {cmd[0] for cmd in commands}
    ):
        print(f"C backend already built ({lib_path})")
        return lib_path

    errors: List[Dict[str, str]] = []
    for cmd in commands:
        try:
            completed = subprocess.run(
                cmd,
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
        except FileNotFoundError as exc:
            errors.append({"cmd": " ".join(cmd), "error": str(exc), "stdout": "", "stderr": ""})
            continue
        except subprocess.CalledProcessError as exc:
            errors.append(
                {
                    "cmd": " ".join(cmd),
                    "error": str(exc),
                    "stdout": exc.stdout or "",
                    "stderr": exc.stderr or "",
                }
            )
            continue

        compiler_id = cmd[0]
        meta = {
            "hash": source_hash,
            "system": system,
            "machine": machine,
            "bits": bits,
            "compiler": compiler_id,
            "built_at": datetime.now(timezone.utc).isoformat(),
        }
        write_json(BUILD_META, meta)
        print(f"Built C backend using {compiler_id}: {lib_path}")
        if completed.stdout:
            print(completed.stdout.strip())
        if completed.stderr:
            print(completed.stderr.strip())
        return lib_path

    if errors:
        last = errors[-1]
        message_lines = [
            "Failed to build C backend.",
            f"Last command: {last['cmd']}",
            f"error: {last['error']}",
        ]
        if last["stdout"]:
            message_lines.append("stdout:\n" + last["stdout"])
        if last["stderr"]:
            message_lines.append("stderr:\n" + last["stderr"])
        raise RuntimeError("\n".join(message_lines))
    raise RuntimeError("No suitable compiler command found for building backend.")


def ensure_backend() -> Path | None:
    lib_path, commands = get_backend_spec()
    source_hash = compute_source_hash()
    system, machine, bits = probe_arch()
    meta = read_json(BUILD_META)
    expected_compilers = {cmd[0] for cmd in commands}
    if not lib_path.exists():
        return build_backend()
    if (
        meta.get("hash") != source_hash
        or meta.get("system") != system
        or meta.get("machine") != machine
        or meta.get("bits") != bits
        or meta.get("compiler") not in expected_compilers
    ):
        return build_backend()
    return lib_path


def save_env() -> None:
    env: Dict[str, Any] = {
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "python_version": sys.version,
        "pointer_bits": struct.calcsize("P") * 8,
    }

    try:
        import numpy  # type: ignore

        env["numpy_version"] = numpy.__version__
    except ImportError:
        pass

    try:
        import matplotlib as _matplotlib  # type: ignore

        env["matplotlib_version"] = _matplotlib.__version__
    except ImportError:
        pass

    meta = read_json(BUILD_META)
    if meta:
        env["compiler"] = meta.get("compiler")
        env["build"] = meta

    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        commit = result.stdout.strip()
        if commit:
            env["git_commit"] = commit
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    write_json(ENV_JSON, env)
    print(f"Wrote {ENV_JSON}")


def clean_outputs() -> None:
    for path in (FIG_DIR, TABLE_DIR, BUILD_DIR):
        if path.exists():
            shutil.rmtree(path)
            print(f"Removed {path}")

DEFAULT_DT = 0.5
DEFAULT_T_END = 1440.0
DEFAULT_PLOT_WINDOW = 240.0

NUTRITION_DEFAULTS: Dict[str, float] = {
    "Vd_dL": 120.0,
    "hepatic_first_pass": 0.25,
    "f_app_base": 0.22,
    "beta_fiber": 0.06,
    "beta_fat": 0.05,
    "kfast_base": 0.45,
    "kslow_base": 0.07,
    "alpha_prot": 0.10,
    "kprot": 0.05,
}

nutrition_settings: Dict[str, float] = NUTRITION_DEFAULTS.copy()

TARGET_PEAK = 50.0
TARGET_TOL = 5.0
MAX_CAL_STEPS = 4
FACTOR_MIN, FACTOR_MAX = 0.1, 200.0

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


@dataclass
class DessertRun:
    name: str
    result: SimulationResult
    mode: str
    converged: bool
    specs: Dict[str, float | str | None | DoseProfile]


class Simulator:
    def __init__(self) -> None:
        _, lib = load_backend()
        self._lib = lib
        self._params_struct = None
        self._has_extended = False

        if self._lib is not None:
            class Params(ctypes.Structure):
                _fields_ = [
                    ("S_G_min1", ctypes.c_double),
                    ("p2_min1", ctypes.c_double),
                    ("p3_min1_per_uU_per_mL", ctypes.c_double),
                    ("phi_G_uU_mL_min1_per_mg_dL", ctypes.c_double),
                    ("G_thr_mg_dL", ctypes.c_double),
                    ("n_min1", ctypes.c_double),
                    ("Gb_mg_dL", ctypes.c_double),
                    ("Ib_uU_mL", ctypes.c_double),
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
            params["S_G_min1"],
            params["p2_min1"],
            params["p3_min1_per_uU_per_mL"],
            params["phi_G_uU_mL_min1_per_mg_dL"],
            params["G_thr_mg_dL"],
            params["n_min1"],
            params["Gb_mg_dL"],
            params["Ib_uU_mL"],
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
            params["S_G_min1"],
            params["p2_min1"],
            params["p3_min1_per_uU_per_mL"],
            params["phi_G_uU_mL_min1_per_mg_dL"],
            params["G_thr_mg_dL"],
            params["n_min1"],
            params["Gb_mg_dL"],
            params["Ib_uU_mL"],
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

        G = params["Gb_mg_dL"]
        X = 0.0
        I = params["Ib_uU_mL"]
        glucose: List[float] = []
        insulin: List[float] = []

        def deriv(t: float, g: float, x: float, ins: float) -> Tuple[float, float, float]:
            fast = profile.Afast * math.exp(-profile.kfast * t)
            slow = profile.Aslow * math.exp(-profile.kslow * t)
            D = fast + slow
            secretion = params["phi_G_uU_mL_min1_per_mg_dL"] * max(0.0, g - params["G_thr_mg_dL"])
            iprot = profile.Aprot * math.exp(-profile.kprot * t)
            dG = -(params["S_G_min1"] + x) * g + params["S_G_min1"] * params["Gb_mg_dL"] + D
            dX = -params["p2_min1"] * x + params["p3_min1_per_uU_per_mL"] * (ins - params["Ib_uU_mL"])
            dI = -params["n_min1"] * (ins - params["Ib_uU_mL"]) + secretion + iprot
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
    lib_path, _ = get_backend_spec()
    lib = try_load_library(lib_path)
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
    global nutrition_settings
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    required = {
        "S_G_min1",
        "p2_min1",
        "p3_min1_per_uU_per_mL",
        "phi_G_uU_mL_min1_per_mg_dL",
        "G_thr_mg_dL",
        "n_min1",
        "Gb_mg_dL",
        "Ib_uU_mL",
    }
    missing = required.difference(data)
    if missing:
        raise KeyError(f"Missing config keys: {sorted(missing)}")
    params = {key: float(data[key]) for key in required}
    params["dt"] = float(data.get("dt", DEFAULT_DT))
    params["t_end"] = float(data.get("t_end", DEFAULT_T_END))
    params["plot_window_min"] = float(data.get("plot_window_min", DEFAULT_PLOT_WINDOW))
    for key, default in NUTRITION_DEFAULTS.items():
        nutrition_settings[key] = float(data.get(key, default))
    return params


def clip(value: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, value))


def nutrition_to_profile(
    carbs_g: float,
    sugars_g: float,
    fiber_g: float,
    fat_g: float,
    protein_g: float,
    settings: Dict[str, float],
) -> Tuple[DoseProfile, Dict[str, float]]:
    carbs_g = max(0.0, carbs_g)
    sugars_g = max(0.0, sugars_g)
    fiber_g = max(0.0, fiber_g)
    fat_g = max(0.0, fat_g)
    protein_g = max(0.0, protein_g)

    avail_carbs_g = max(0.0, carbs_g - 0.5 * fiber_g)
    f_fast = 0.0 if carbs_g <= 0.0 else clip(sugars_g / carbs_g, 0.0, 1.0)
    f_app = clip(
        settings["f_app_base"] * (1.0 - settings["beta_fiber"] * (fiber_g / 10.0)),
        0.05,
        0.60,
    )
    k_mod = 1.0 / (
        1.0
        + settings["beta_fat"] * (fat_g / 10.0)
        + settings["beta_fiber"] * (fiber_g / 10.0)
    )
    kfast = settings["kfast_base"] * k_mod
    kslow = settings["kslow_base"] * k_mod
    dose_total_mgdL = (
        (avail_carbs_g * f_app * 1000.0) / settings["Vd_dL"]
    ) * (1.0 - settings["hepatic_first_pass"])
    dose_fast = dose_total_mgdL * f_fast
    dose_slow = dose_total_mgdL * (1.0 - f_fast)
    Afast = kfast * dose_fast
    Aslow = kslow * dose_slow
    Aprot = settings["alpha_prot"] * protein_g

    profile = DoseProfile(Afast, kfast, Aslow, kslow, Aprot, settings["kprot"])
    extras = {
        "dose_mgdL": profile.total_dose(),
        "f_fast": f_fast,
        "f_app": f_app,
        "k_mod": k_mod,
        "kfast": kfast,
        "kslow": kslow,
    }
    return profile, extras


def build_dessert_entry(
    name: str,
    data: Dict[str, Any],
    settings: Dict[str, float],
) -> Dict[str, float | str | None | DoseProfile]:
    carbs_g = float(data.get("carbs_g", 0.0) or 0.0)
    sugars_g = float(data.get("sugars_g", 0.0) or 0.0)
    fiber_g = float(data.get("fiber_g", 0.0) or 0.0)
    fat_g = float(data.get("fat_g", 0.0) or 0.0)
    protein_g = float(data.get("protein_g", 0.0) or 0.0)
    profile, extras = nutrition_to_profile(
        carbs_g, sugars_g, fiber_g, fat_g, protein_g, settings
    )
    entry: Dict[str, float | str | None | DoseProfile] = {
        "profile": profile,
        "carbs_g": carbs_g,
        "sugars_g": sugars_g,
        "fiber_g": fiber_g,
        "fat_g": fat_g,
        "protein_g": protein_g,
        "portion_g": (
            float(data.get("portion_g", 0.0))
            if data.get("portion_g") is not None
            else None
        ),
        "barcode": data.get("barcode"),
    }
    entry.update(extras)
    return entry


def load_dessert_catalog(path: Path) -> Dict[str, Dict[str, float | str | None | DoseProfile]]:
    with path.open("r", encoding="utf-8") as handle:
        raw_data = yaml.safe_load(handle) or {}
    if not isinstance(raw_data, dict):
        raise ValueError(f"Dessert catalog {path} must be a mapping")
    desserts: Dict[str, Dict[str, float | str | None | DoseProfile]] = {}
    for name, data in sorted(raw_data.items()):
        if not isinstance(data, dict):
            raise ValueError(f"Dessert entry {name} must be a mapping")
        desserts[name] = build_dessert_entry(name, data, nutrition_settings)
    return desserts


def load_frozen_catalog(directory: Path) -> Dict[str, Dict[str, float | str | None | DoseProfile]]:
    desserts: Dict[str, Dict[str, float | str | None | DoseProfile]] = {}
    for yaml_path in sorted(directory.glob("*.yaml")):
        with yaml_path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
        if not isinstance(data, dict):
            raise ValueError(f"Frozen dessert file {yaml_path} must be a mapping")
        name = str(data.get("name") or yaml_path.stem)
        desserts[name] = build_dessert_entry(name, data, nutrition_settings)
    return desserts


def clamp(value: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, value))


def simulate_desserts(
    simulator: Simulator,
    params: Dict[str, float],
    dt: float,
    t_end: float,
    desserts: Dict[str, Dict[str, float | str | None | DoseProfile]],
    calibrate: bool,
) -> List[DessertRun]:
    runs: List[DessertRun] = []
    for name, specs in sorted(desserts.items()):
        base_profile = specs["profile"]
        if calibrate:
            result, converged = calibrate_profile(
                simulator, params, dt, t_end, base_profile
            )
            mode = "calibrated"
        else:
            result = simulator.run(params, dt, t_end, base_profile)
            converged = False
            mode = "dose-driven"
        runs.append(
            DessertRun(
                name=name,
                result=result,
                mode=mode,
                converged=converged,
                specs=specs,
            )
        )
    return runs


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
        peak_delta = max(result.glucose) - params["Gb_mg_dL"]
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


def compute_metrics(
    result: SimulationResult, params: Dict[str, float]
) -> Dict[str, float]:
    times = result.times
    glucose = result.glucose
    insulin = result.insulin
    peak_g = max(glucose)
    t_peak_g = times[glucose.index(peak_g)] if glucose else float("nan")
    peak_i = max(insulin)
    t_peak_i = times[insulin.index(peak_i)] if insulin else float("nan")
    gb = params["Gb_mg_dL"]
    nadir_window = [g for t, g in zip(times, glucose) if t <= 240.0]
    nadir_g = min(nadir_window) if nadir_window else float("nan")

    baseline_return = float("nan")
    for t, g in zip(times, glucose):
        if abs(g - gb) <= 5.0:
            baseline_return = t
            break

    iauc_120 = integrate_trapezoid(times, glucose, 120.0, subtract=gb, clamp_zero=True)
    iauc_240 = integrate_trapezoid(times, glucose, 240.0, subtract=gb, clamp_zero=True)
    aucg_240 = integrate_trapezoid(times, glucose, 240.0)
    auci_240 = integrate_trapezoid(times, insulin, 240.0)

    return {
        "peak_G": peak_g,
        "t_peak_G": t_peak_g,
        "peak_delta": peak_g - gb,
        "nadir_G_0_240": nadir_g,
        "baseline_return_min": baseline_return,
        "iAUC_0_120": iauc_120,
        "iAUC_0_240": iauc_240,
        "AUCG_0_240": aucg_240,
        "peak_I": peak_i,
        "t_peak_I": t_peak_i,
        "AUCI_0_240": auci_240,
    }


def ensure_dirs(*paths: os.PathLike[str] | str) -> None:
    for path in paths:
        os.makedirs(path, exist_ok=True)


def ensure_directories() -> None:
    ensure_dirs(FIG_DIR, TABLE_DIR)


def relative_to_root(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return os.path.relpath(path, ROOT)


def save_dual(
    fig: matplotlib.figure.Figure,
    base_name: str,
    t: Sequence[float],
    ys: Sequence[Sequence[float]],
    window: float,
    zoom_dir: os.PathLike[str] | str,
    full_dir: os.PathLike[str] | str,
    dpi: int,
) -> None:
    import numpy as np

    ensure_dirs(zoom_dir, full_dir)
    t_arr = np.asarray(t, dtype=float)
    y_arrays = [np.asarray(y, dtype=float) for y in ys]
    if t_arr.size == 0:
        t_arr = np.array([0.0, window], dtype=float)
        if not y_arrays:
            y_arrays = [np.zeros_like(t_arr)]
        else:
            y_arrays = [np.zeros_like(t_arr) for _ in y_arrays]
    mask = t_arr <= window
    if not mask.any():
        mask = np.ones_like(t_arr, dtype=bool)
    ystack = np.column_stack([arr[mask] for arr in y_arrays])
    ymin = float(np.nanmin(ystack))
    ymax = float(np.nanmax(ystack))
    if not (np.isfinite(ymin) and np.isfinite(ymax)):
        ymin, ymax = 0.0, 1.0
    if ymax <= ymin:
        ymax = ymin + 1.0
    pad = 0.02 * (ymax - ymin)
    ymin -= pad
    ymax += pad
    ax = fig.axes[0]
    ax.set_xlim(0.0, window)
    ax.set_ylim(ymin, ymax)
    fig.tight_layout()
    fig.savefig(os.path.join(str(zoom_dir), base_name), dpi=dpi)

    ax.set_xlim(0.0, t_arr[-1])
    ax.relim()
    ax.autoscale(axis="y", tight=False)
    fig.tight_layout()
    fig.savefig(
        os.path.join(str(full_dir), base_name.replace(".png", "_full.png")),
        dpi=dpi,
    )


def save_line_plot(
    x: Sequence[float],
    y: Sequence[float],
    base_name: str,
    title: str,
    ylabel: str,
    window: float,
    zoom_dir: os.PathLike[str] | str,
    full_dir: os.PathLike[str] | str,
    dpi: int,
) -> Tuple[Path, Path]:
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    ax.plot(x, y, lw=1.8)
    ax.set_title(title)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    save_dual(fig, base_name, x, [y], window, zoom_dir, full_dir, dpi)
    plt.close(fig)
    return Path(zoom_dir) / base_name, Path(full_dir) / base_name.replace(".png", "_full.png")


def save_overlay(
    curves: Iterable[Tuple[str, Sequence[float], Sequence[float]]],
    base_name: str,
    title: str,
    ylabel: str,
    window: float,
    zoom_dir: os.PathLike[str] | str,
    full_dir: os.PathLike[str] | str,
    dpi: int,
) -> Tuple[Path, Path]:
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    t_ref: Sequence[float] | None = None
    y_values: List[Sequence[float]] = []
    has_curves = False
    for name, x, y in curves:
        ax.plot(x, y, lw=1.5, label=name)
        if t_ref is None:
            t_ref = x
        y_values.append(y)
        has_curves = True
    if t_ref is None:
        t_ref = [0.0, window]
        y_values = [[0.0, 0.0]]
    ax.set_title(title)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    if has_curves:
        ax.legend()
    save_dual(fig, base_name, t_ref, y_values, window, zoom_dir, full_dir, dpi)
    plt.close(fig)
    return Path(zoom_dir) / base_name, Path(full_dir) / base_name.replace(".png", "_full.png")


def make_dose_curves(profile: DoseProfile, duration: float = 360.0, step: float = 0.5) -> Tuple[List[float], List[float], List[float], List[float]]:
    n = int(round(duration / step)) + 1
    times = [i * step for i in range(n)]
    fast = [profile.Afast * math.exp(-profile.kfast * t) for t in times]
    slow = [profile.Aslow * math.exp(-profile.kslow * t) for t in times]
    total = [fast[i] + slow[i] for i in range(n)]
    return times, fast, slow, total


def save_d_components(
    name: str,
    profile: DoseProfile,
    plot_window: float,
    t_end: float,
    zoom_dir: os.PathLike[str] | str,
    full_dir: os.PathLike[str] | str,
    dpi: int,
) -> Tuple[Path, Path]:
    times, fast, slow, total = make_dose_curves(profile, duration=t_end)
    base_name = f"D_components_{name}.png"
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    ax.plot(times, fast, label="fast", lw=1.6)
    ax.plot(times, slow, label="slow", lw=1.6)
    ax.plot(times, total, label="total", lw=2.0, linestyle="--")
    ax.set_title(f"{name.title()} appearance components")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("D(t) (mg/dL·min⁻¹)")
    ax.grid(True)
    ax.legend()
    save_dual(fig, base_name, times, [fast, slow, total], plot_window, zoom_dir, full_dir, dpi)
    plt.close(fig)
    return Path(zoom_dir) / base_name, Path(full_dir) / base_name.replace(".png", "_full.png")


def write_summary(rows: List[Dict[str, object]]) -> None:
    base_fields = [
        "name",
        "mode",
        "backend",
        "input_source",
        "Gb_mg_dL",
        "Ib_uU_mL",
        "peak_G",
        "t_peak_G",
        "nadir_G_0_240",
        "baseline_return_min",
        "iAUC_0_120",
        "iAUC_0_240",
        "AUCG_0_240",
        "peak_I",
        "t_peak_I",
        "AUCI_0_240",
        "carbs_g",
        "sugars_g",
        "fiber_g",
        "fat_g",
        "protein_g",
        "f_fast",
        "f_app",
        "k_mod",
        "kfast",
        "kslow",
        "dose_mgdL",
        "Afast",
        "Aslow",
        "Aprot",
        "peak_delta",
        "kprot",
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
    plot_window_min: float,
    plot_zoom_dir: Path,
    plot_full_dir: Path,
    dpi: int,
    desserts: Dict[str, Dict[str, float | str | None | DoseProfile]],
    calibrate: bool,
    input_source: str,
    emit_plots: bool,
    summary_only: bool,
) -> List[DessertRun]:
    ensure_directories()
    if emit_plots:
        ensure_dirs(plot_zoom_dir, plot_full_dir)
    summary_rows: List[Dict[str, object]] = []
    glucose_curves: List[Tuple[str, List[float], List[float]]] = []
    insulin_curves: List[Tuple[str, List[float], List[float]]] = []
    dose_curves: List[Tuple[str, List[float], List[float]]] = []
    manifest_entries: List[Dict[str, str]] = []

    mode_label = "calibrated" if calibrate else "dose-driven"
    print(f"Running dessert pipeline ({mode_label}) with backend={simulator.label}...")

    runs = simulate_desserts(simulator, params, dt, t_end, desserts, calibrate)

    for run in runs:
        name = run.name
        specs = run.specs
        result = run.result
        mode = run.mode
        final_profile = result.profile

        metrics = compute_metrics(result, params)
        if emit_plots:
            glucose_curves.append((name, result.times, result.glucose))
            insulin_curves.append((name, result.times, result.insulin))
            dose_times, _, _, dose_total = make_dose_curves(final_profile, duration=t_end)
            dose_curves.append((name, dose_times, dose_total))

        final_dose = final_profile.total_dose()
        summary_row: Dict[str, object] = {
            "name": name,
            "backend": simulator.label,
            "mode": mode,
            "input_source": input_source,
            "Gb_mg_dL": params["Gb_mg_dL"],
            "Ib_uU_mL": params["Ib_uU_mL"],
            "peak_G": metrics["peak_G"],
            "t_peak_G": metrics["t_peak_G"],
            "nadir_G_0_240": metrics["nadir_G_0_240"],
            "baseline_return_min": metrics["baseline_return_min"],
            "iAUC_0_120": metrics["iAUC_0_120"],
            "iAUC_0_240": metrics["iAUC_0_240"],
            "AUCG_0_240": metrics["AUCG_0_240"],
            "peak_I": metrics["peak_I"],
            "t_peak_I": metrics["t_peak_I"],
            "AUCI_0_240": metrics["AUCI_0_240"],
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
            "peak_delta": metrics["peak_delta"],
        }
        if calibrate:
            summary_row["calibration_converged"] = run.converged
        summary_rows.append(summary_row)

        print(
            f"{name:12s} | backend={simulator.label:6s} | mode={mode:11s} | "
            f"source={input_source:8s} | peakΔG={metrics['peak_delta']:+6.2f} mg/dL @ {metrics['t_peak_G']:.1f} min | "
            f"iAUC₀₋₁₂₀={metrics['iAUC_0_120']:.1f} | peakI={metrics['peak_I']:.2f}"
        )

        if emit_plots:
            glucose_base = f"{name}_glucose.png"
            insulin_base = f"{name}_insulin.png"
            glucose_zoom, glucose_full = save_line_plot(
                result.times,
                result.glucose,
                glucose_base,
                f"{name.title()} glucose",
                "Glucose (mg/dL)",
                plot_window_min,
                plot_zoom_dir,
                plot_full_dir,
                dpi,
            )
            insulin_zoom, insulin_full = save_line_plot(
                result.times,
                result.insulin,
                insulin_base,
                f"{name.title()} insulin",
                "Insulin (mU/L)",
                plot_window_min,
                plot_zoom_dir,
                plot_full_dir,
                dpi,
            )
            components_zoom, components_full = save_d_components(
                name,
                final_profile,
                plot_window_min,
                t_end,
                plot_zoom_dir,
                plot_full_dir,
                dpi,
            )

            manifest_entries.extend(
                [
                    {
                        "path": relative_to_root(glucose_zoom),
                        "caption": f"{name.title()} glucose curve",
                    },
                    {
                        "path": relative_to_root(glucose_full),
                        "caption": f"{name.title()} glucose curve (full)",
                    },
                    {
                        "path": relative_to_root(insulin_zoom),
                        "caption": f"{name.title()} insulin curve",
                    },
                    {
                        "path": relative_to_root(insulin_full),
                        "caption": f"{name.title()} insulin curve (full)",
                    },
                    {
                        "path": relative_to_root(components_zoom),
                        "caption": f"{name.title()} appearance components",
                    },
                    {
                        "path": relative_to_root(components_full),
                        "caption": f"{name.title()} appearance components (full)",
                    },
                ]
            )

    if emit_plots:
        overlay_glucose_base = "glucose_overlay.png"
        overlay_insulin_base = "insulin_overlay.png"
        overlay_dose_base = "D_overlay.png"

        overlay_glucose, overlay_glucose_full = save_overlay(
            glucose_curves,
            overlay_glucose_base,
            "Glucose overlay",
            "Glucose (mg/dL)",
            plot_window_min,
            plot_zoom_dir,
            plot_full_dir,
            dpi,
        )
        overlay_insulin, overlay_insulin_full = save_overlay(
            insulin_curves,
            overlay_insulin_base,
            "Insulin overlay",
            "Insulin (mU/L)",
            plot_window_min,
            plot_zoom_dir,
            plot_full_dir,
            dpi,
        )
        overlay_dose, overlay_dose_full = save_overlay(
            dose_curves,
            overlay_dose_base,
            "Dessert appearance D(t)",
            "Dose (mg/dL·min⁻¹)",
            plot_window_min,
            plot_zoom_dir,
            plot_full_dir,
            dpi,
        )

        manifest_entries.extend(
            [
                {
                    "path": relative_to_root(overlay_glucose),
                    "caption": "Dessert glucose overlay",
                },
                {
                    "path": relative_to_root(overlay_glucose_full),
                    "caption": "Dessert glucose overlay (full)",
                },
                {
                    "path": relative_to_root(overlay_insulin),
                    "caption": "Dessert insulin overlay",
                },
                {
                    "path": relative_to_root(overlay_insulin_full),
                    "caption": "Dessert insulin overlay (full)",
                },
                {
                    "path": relative_to_root(overlay_dose),
                    "caption": "Dessert appearance overlay",
                },
                {
                    "path": relative_to_root(overlay_dose_full),
                    "caption": "Dessert appearance overlay (full)",
                },
            ]
        )

    write_summary(summary_rows)
    if emit_plots and manifest_entries and not summary_only:
        write_manifest(manifest_entries)

    return runs


def run_sanity(
    simulator: Simulator,
    params: Dict[str, float],
    dt: float,
    t_end: float,
    desserts: Dict[str, Dict[str, float | str | None | DoseProfile]],
    runs: List[DessertRun] | None,
    calibrate: bool,
) -> None:
    print(f"Running sanity checks with backend={simulator.label}...")

    if runs is None:
        runs = simulate_desserts(simulator, params, dt, t_end, desserts, calibrate)

    run_map = {run.name: run for run in runs}
    ordered_names = sorted(desserts.keys())
    dt_check: Dict[str, Any]

    if ordered_names:
        first_name = ordered_names[0]
        base_profile = desserts[first_name]["profile"]
        base_result = simulator.run(params, dt, t_end, base_profile)
        half_result = simulator.run(params, dt / 2.0, t_end, base_profile)
        base_metrics = compute_metrics(base_result, params)
        half_metrics = compute_metrics(half_result, params)
        peak_diff = abs(base_metrics["peak_G"] - half_metrics["peak_G"])
        t_peak_diff = abs(base_metrics["t_peak_G"] - half_metrics["t_peak_G"])
        passed = peak_diff <= 1.0 and t_peak_diff <= 2.0
        dt_check = {
            "dessert": first_name,
            "peak_G_base": base_metrics["peak_G"],
            "peak_G_dt_half": half_metrics["peak_G"],
            "t_peak_base": base_metrics["t_peak_G"],
            "t_peak_dt_half": half_metrics["t_peak_G"],
            "peak_diff": peak_diff,
            "t_peak_diff": t_peak_diff,
            "passed": passed,
        }
        status = "PASS" if passed else "FAIL"
        print(
            f"dt-halving ({first_name}): peakΔ={peak_diff:.3f} mg/dL, tΔ={t_peak_diff:.3f} min -> {status}"
        )
        if not passed:
            raise AssertionError(
                "dt-halving check failed: Δpeak_G>1 mg/dL or Δt_peak>2 min"
            )
    else:
        dt_check = {"status": "skipped", "reason": "no desserts configured"}
        print("dt-halving check skipped (no desserts configured)")

    end_state_entries: List[Dict[str, Any]] = []
    gb = params["Gb_mg_dL"]
    ib = params["Ib_uU_mL"]
    for name, run in sorted(run_map.items()):
        final_g = run.result.glucose[-1]
        final_i = run.result.insulin[-1]
        g_diff = abs(final_g - gb)
        i_diff = abs(final_i - ib)
        passed = g_diff <= 5.0 and i_diff <= 0.5
        end_state_entries.append(
            {
                "name": name,
                "mode": run.mode,
                "final_G": final_g,
                "final_I": final_i,
                "delta_G": g_diff,
                "delta_I": i_diff,
                "passed": passed,
            }
        )
        status = "PASS" if passed else "FAIL"
        print(
            f"end-state ({name}): |ΔG|={g_diff:.3f} mg/dL, |ΔI|={i_diff:.3f} μU/mL -> {status}"
        )

    report = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "backend": simulator.label,
        "pipeline_mode": "calibrated" if calibrate else "dose-driven",
        "checks": {
            "dt_halving": dt_check,
            "end_state": end_state_entries,
        },
    }

    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    write_json(SANITY_JSON, report)
    print(f"Wrote {SANITY_JSON}")


def run_sensitivity_analysis(
    simulator: Simulator,
    params: Dict[str, float],
    dt: float,
    t_end: float,
    desserts: Dict[str, Dict[str, float | str | None | DoseProfile]],
    calibrate: bool,
) -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, object]] = []
    for name, specs in sorted(desserts.items()):
        base_macros = {
            "carbs_g": float(specs.get("carbs_g") or 0.0),
            "sugars_g": float(specs.get("sugars_g") or 0.0),
            "fiber_g": float(specs.get("fiber_g") or 0.0),
            "fat_g": float(specs.get("fat_g") or 0.0),
            "protein_g": float(specs.get("protein_g") or 0.0),
        }
        for macro in ("fat_g", "fiber_g", "protein_g"):
            baseline_value = base_macros[macro]
            for label, factor in (("-20%", 0.8), ("+20%", 1.2)):
                mutated = dict(base_macros)
                mutated[macro] = baseline_value * factor
                profile, extras = nutrition_to_profile(
                    mutated["carbs_g"],
                    mutated["sugars_g"],
                    mutated["fiber_g"],
                    mutated["fat_g"],
                    mutated["protein_g"],
                    nutrition_settings,
                )
                if calibrate:
                    result, converged = calibrate_profile(
                        simulator, params, dt, t_end, profile
                    )
                    mode = "calibrated"
                else:
                    result = simulator.run(params, dt, t_end, profile)
                    converged = False
                    mode = "dose-driven"
                metrics = compute_metrics(result, params)
                rows.append(
                    {
                        "name": name,
                        "mode": mode,
                        "macro": macro,
                        "delta": label,
                        "factor": factor,
                        "input_source": "sensitivity",
                        "Gb_mg_dL": params["Gb_mg_dL"],
                        "Ib_uU_mL": params["Ib_uU_mL"],
                        "carbs_g": mutated["carbs_g"],
                        "sugars_g": mutated["sugars_g"],
                        "fiber_g": mutated["fiber_g"],
                        "fat_g": mutated["fat_g"],
                        "protein_g": mutated["protein_g"],
                        "dose_mgdL": extras["dose_mgdL"],
                        "f_fast": extras["f_fast"],
                        "f_app": extras["f_app"],
                        "k_mod": extras["k_mod"],
                        "kfast": extras["kfast"],
                        "kslow": extras["kslow"],
                        "Afast": profile.Afast,
                        "Aslow": profile.Aslow,
                        "Aprot": profile.Aprot,
                        "kprot": profile.kprot,
                        "peak_G": metrics["peak_G"],
                        "peak_delta": metrics["peak_delta"],
                        "t_peak_G": metrics["t_peak_G"],
                        "nadir_G_0_240": metrics["nadir_G_0_240"],
                        "baseline_return_min": metrics["baseline_return_min"],
                        "iAUC_0_120": metrics["iAUC_0_120"],
                        "iAUC_0_240": metrics["iAUC_0_240"],
                        "AUCG_0_240": metrics["AUCG_0_240"],
                        "peak_I": metrics["peak_I"],
                        "t_peak_I": metrics["t_peak_I"],
                        "AUCI_0_240": metrics["AUCI_0_240"],
                        "calibration_converged": converged if calibrate else None,
                    }
                )
    path = TABLE_DIR / "sensitivity.csv"
    fieldnames: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"Wrote {path}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Dessert glycemia simulator")
    parser.add_argument(
        "--reproduce",
        action="store_true",
        help="Clean, rebuild, and run all frozen YAML configs",
    )
    parser.add_argument(
        "--calibrate",
        action="store_true",
        help="Equalize glucose peaks by rescaling appearance amplitudes",
    )
    parser.add_argument(
        "--sensitivity",
        action="store_true",
        help="Run ±20% fat/fiber/protein sensitivity analysis",
    )
    parser.add_argument("--window", type=int, help="Override plot window (minutes)")
    parser.add_argument("--dpi", type=int, default=150, help="Override plot DPI (default 150)")
    parser.add_argument("--no-plots", action="store_true", help="Skip plotting outputs")
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Only generate tables (implies --no-plots)",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)

    if args.summary_only:
        args.no_plots = True
    calibrate = args.calibrate
    if args.reproduce and calibrate:
        print("Reproduction mode enforces dose-driven simulations; ignoring --calibrate.")
        calibrate = False

    params = load_config(CONFIG_PATH)
    dt = params.pop("dt")
    t_end = params.pop("t_end")
    plot_window_min = float(params.pop("plot_window_min", DEFAULT_PLOT_WINDOW))
    plot_zoom_dir_raw = params.pop("plot_zoom_dir", "figures/zoom")
    plot_full_dir_raw = params.pop("plot_full_dir", "figures/full")
    if args.window is not None:
        plot_window_min = float(args.window)
    plot_window_min = max(1.0, plot_window_min)

    plot_zoom_dir = Path(plot_zoom_dir_raw)
    if not plot_zoom_dir.is_absolute():
        plot_zoom_dir = ROOT / plot_zoom_dir
    plot_full_dir = Path(plot_full_dir_raw)
    if not plot_full_dir.is_absolute():
        plot_full_dir = ROOT / plot_full_dir

    dpi = max(1, args.dpi)

    if args.reproduce:
        clean_outputs()
        try:
            build_backend()
        except RuntimeError as exc:
            raise SystemExit(str(exc)) from exc
        input_source = "frozen"
        desserts = load_frozen_catalog(ROOT / "configs" / "frozen")
    else:
        try:
            ensure_backend()
        except RuntimeError as exc:
            raise SystemExit(str(exc)) from exc
        input_source = "yaml"
        desserts = load_dessert_catalog(DESSERTS_CONFIG_PATH)

    save_env()

    simulator = Simulator()
    emit_plots = not args.no_plots and not args.summary_only

    pipeline_runs = run_pipeline(
        simulator,
        params,
        dt,
        t_end,
        plot_window_min,
        plot_zoom_dir,
        plot_full_dir,
        dpi,
        desserts,
        calibrate=calibrate,
        input_source=input_source,
        emit_plots=emit_plots,
        summary_only=args.summary_only,
    )
    run_sanity(simulator, params, dt, t_end, desserts, pipeline_runs, calibrate)

    if args.sensitivity:
        run_sensitivity_analysis(simulator, params, dt, t_end, desserts, calibrate)


if __name__ == "__main__":
    main(sys.argv[1:])
