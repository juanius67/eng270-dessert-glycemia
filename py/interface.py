"""ctypes bridge to the C minimal-model integrator."""
from __future__ import annotations

import ctypes as _ct
import sys
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
_C_DIR = _ROOT / "C"

_LIB_NAMES = {
    "win32": "model.dll",
    "cygwin": "model.dll",
}.get(sys.platform, "libmodel.so"), "model.dll"


class _Params(_ct.Structure):
    _fields_ = [
        ("p1", _ct.c_double),
        ("p2", _ct.c_double),
        ("p3", _ct.c_double),
        ("p4", _ct.c_double),
        ("p5", _ct.c_double),
        ("p6", _ct.c_double),
        ("Gb", _ct.c_double),
        ("Ib", _ct.c_double),
        ("A", _ct.c_double),
        ("k", _ct.c_double),
    ]


def _load_library() -> _ct.CDLL:
    preferred, fallback = _LIB_NAMES
    for candidate in (preferred, fallback):
        path = _C_DIR / candidate
        if path.exists():
            lib = _ct.CDLL(str(path))
            break
    else:
        raise FileNotFoundError(
            "Could not locate libmodel shared library. "
            "Build it first with `make -C C`."
        )

    lib.simulate.argtypes = [
        _ct.POINTER(_ct.c_double),
        _ct.POINTER(_ct.c_double),
        _ct.c_int,
        _ct.c_double,
        _ct.POINTER(_Params),
    ]
    lib.simulate.restype = None
    return lib


_LIB = None


def _get_library() -> _ct.CDLL:
    global _LIB
    if _LIB is None:
        _LIB = _load_library()
    return _LIB


def simulate_c(params: Dict[str, float], A: float, k: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run the C integrator and return ``(t, G, I)`` arrays."""

    dt = float(params["dt"])
    t_end = float(params["t_end"])
    nsteps = int(t_end / dt) + 1

    lib = _get_library()

    g_arr = np.empty(nsteps, dtype=np.float64)
    i_arr = np.empty(nsteps, dtype=np.float64)

    params_struct = _Params(
        float(params["p1"]),
        float(params["p2"]),
        float(params["p3"]),
        float(params["p4"]),
        float(params["p5"]),
        float(params["p6"]),
        float(params["Gb"]),
        float(params["Ib"]),
        float(A),
        float(k),
    )

    lib.simulate(
        g_arr.ctypes.data_as(_ct.POINTER(_ct.c_double)),
        i_arr.ctypes.data_as(_ct.POINTER(_ct.c_double)),
        nsteps,
        dt,
        _ct.byref(params_struct),
    )

    t_arr = np.linspace(0.0, t_end, nsteps, dtype=np.float64)
    return t_arr, g_arr, i_arr
