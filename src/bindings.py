"""ctypes bindings for the Bergman minimal model C backend."""
from __future__ import annotations

import ctypes
import sys
from ctypes import Structure, c_double
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parents[1]
C_DIR = ROOT / "C"


def _candidate_libraries() -> list[Path]:
    if sys.platform.startswith("win"):
        names = ["model.dll", "libmodel.so", "libmodel.dylib"]
    elif sys.platform == "darwin":
        names = ["libmodel.dylib", "libmodel.so", "model.dll"]
    else:
        names = ["libmodel.so", "libmodel.dylib", "model.dll"]
    return [C_DIR / name for name in names]


def _load_library() -> Optional[ctypes.CDLL]:
    for candidate in _candidate_libraries():
        if not candidate.exists():
            continue
        try:
            return ctypes.CDLL(str(candidate))
        except OSError:
            continue
    return None


lib: Optional[ctypes.CDLL] = _load_library()


class BergmanParams(Structure):
    _fields_ = [
        ("S_G", c_double),
        ("p2", c_double),
        ("p3", c_double),
        ("n", c_double),
        ("Gb", c_double),
        ("Ib", c_double),
    ]


if lib is not None:
    lib.set_params.argtypes = [BergmanParams]
    lib.set_params.restype = None
    if hasattr(lib, "simulate_dual"):
        lib.simulate_dual.argtypes = [
            ctypes.POINTER(c_double),
            ctypes.POINTER(c_double),
            ctypes.c_int,
            c_double,
            c_double,
            c_double,
            c_double,
            c_double,
            c_double,
            c_double,
        ]
        lib.simulate_dual.restype = None


def set_params_from_dict(params: dict[str, float]) -> None:
    """Push YAML parameters into the C core."""
    if lib is None:
        return
    bp = BergmanParams(
        float(params["S_G_min1"]),
        float(params["p2_min1"]),
        float(params["p3_min1_per_uU_mL"]),
        float(params["n_min1"]),
        float(params["Gb_mg_dL"]),
        float(params["Ib_uU_mL"]),
    )
    lib.set_params(bp)
