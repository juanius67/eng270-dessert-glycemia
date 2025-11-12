"""ctypes bindings for the Bergman minimal model RK4 core."""
from __future__ import annotations

import ctypes
import hashlib
import os
import platform
import subprocess
import sys
from ctypes import POINTER, Structure, c_double
from pathlib import Path
from typing import Dict, Optional

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "src"
BUILD_DIR = ROOT / "build"

BUILD_METADATA: Dict[str, object] = {}


def _shared_name() -> str:
    if sys.platform.startswith("win"):
        return "model.dll"
    if sys.platform == "darwin":
        return "libmodel.dylib"
    return "libmodel.so"


def _candidate_commands(output: Path) -> list[list[str]]:
    src = str(SRC_DIR / "model.c")
    include_flag = f"-I{SRC_DIR}" if not sys.platform.startswith("win") else f"/I{SRC_DIR}"
    if sys.platform.startswith("win"):
        return [
            [
                "cl",
                "/nologo",
                "/LD",
                "/DBUILDING_MODEL",
                include_flag,
                src,
                f"/Fe{output}",
            ],
            [
                "gcc",
                "-O3",
                "-shared",
                "-fPIC",
                "-DBUILDING_MODEL",
                include_flag,
                src,
                "-o",
                str(output),
            ],
        ]
    if sys.platform == "darwin":
        return [
            [
                "clang",
                "-O3",
                "-dynamiclib",
                "-fPIC",
                "-DBUILDING_MODEL",
                include_flag,
                src,
                "-o",
                str(output),
            ],
            [
                "gcc",
                "-O3",
                "-dynamiclib",
                "-fPIC",
                "-DBUILDING_MODEL",
                include_flag,
                src,
                "-o",
                str(output),
            ],
        ]
    return [
        [
            "gcc",
            "-O3",
            "-shared",
            "-fPIC",
            "-DBUILDING_MODEL",
            include_flag,
            src,
            "-o",
            str(output),
        ],
        [
            "clang",
            "-O3",
            "-shared",
            "-fPIC",
            "-DBUILDING_MODEL",
            include_flag,
            src,
            "-o",
            str(output),
        ],
    ]


def _build_library() -> Path:
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    output = BUILD_DIR / _shared_name()
    errors: list[str] = []
    for cmd in _candidate_commands(output):
        try:
            completed = subprocess.run(
                cmd,
                cwd=str(ROOT),
                check=True,
                capture_output=True,
                text=True,
            )
        except FileNotFoundError:
            errors.append("missing compiler: " + cmd[0])
            continue
        except subprocess.CalledProcessError as exc:
            message = exc.stderr or exc.stdout or str(exc)
            errors.append(message.strip())
            continue
        BUILD_METADATA.clear()
        BUILD_METADATA.update(
            {
                "compiler": cmd[0],
                "flags": cmd[1:],
                "command": cmd,
                "library": str(output),
            }
        )
        if completed.stdout:
            BUILD_METADATA["stdout"] = completed.stdout.strip()
        if completed.stderr:
            BUILD_METADATA["stderr"] = completed.stderr.strip()
        return output
    raise RuntimeError("; ".join(errors) if errors else "no compiler available")


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def model_source_hash() -> str:
    return _hash_file(SRC_DIR / "model.c")


def _load_library() -> Optional[ctypes.CDLL]:
    lib_path = BUILD_DIR / _shared_name()
    if lib_path.exists():
        try:
            lib = ctypes.CDLL(str(lib_path))
            BUILD_METADATA.setdefault("library", str(lib_path))
            return lib
        except OSError:
            pass
    try:
        lib_path = _build_library()
    except RuntimeError as exc:
        raise RuntimeError(f"Failed to build C backend: {exc}") from exc
    lib = ctypes.CDLL(str(lib_path))
    return lib


class BergmanParams(Structure):
    _fields_ = [
        ("S_G", c_double),
        ("p2", c_double),
        ("p3", c_double),
        ("n", c_double),
        ("Gb", c_double),
        ("Ib", c_double),
    ]


lib: Optional[ctypes.CDLL] = _load_library()

if lib is not None:
    lib.set_params.argtypes = [BergmanParams]
    lib.set_params.restype = None
    lib.step.argtypes = [
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
    ]
    lib.step.restype = None


def set_params_from_dict(params: Dict[str, float]) -> None:
    if lib is None:
        raise RuntimeError("C library not loaded; cannot set parameters")
    bp = BergmanParams(
        float(params["S_G_min1"]),
        float(params["p2_min1"]),
        float(params["p3_min1_per_uU_mL"]),
        float(params["n_min1"]),
        float(params["Gb_mg_dL"]),
        float(params["Ib_uU_mL"]),
    )
    lib.set_params(bp)


def get_build_metadata() -> Dict[str, object]:
    metadata = dict(BUILD_METADATA)
    metadata.setdefault("library", str(BUILD_DIR / _shared_name()))
    metadata.setdefault("compiler", None)
    metadata.setdefault("flags", [])
    metadata.setdefault("command", [])
    metadata["source_hash"] = model_source_hash()
    metadata["platform"] = platform.platform()
    metadata["python"] = sys.version
    return metadata
