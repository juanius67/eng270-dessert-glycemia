from __future__ import annotations  # WRITTEN BY AI
"""AI-assisted infrastructure helpers for environment metadata and JSON output."""  # WRITTEN BY AI
import json  # WRITTEN BY AI
import platform  # WRITTEN BY AI
from pathlib import Path  # WRITTEN BY AI
from typing import Any, Dict  # WRITTEN BY AI

import numpy as np  # WRITTEN BY AI


def collect_env_metadata() -> Dict[str, Any]:  # WRITTEN BY AI
    """Collect metadata about the current Python and OS environment."""  # WRITTEN BY AI
    build = platform.python_build()  # WRITTEN BY AI
    return {  # WRITTEN BY AI
        "python": platform.python_version(),  # WRITTEN BY AI
        "python_build": " ".join(build),  # WRITTEN BY AI
        "numpy": np.__version__,  # WRITTEN BY AI
        "os": platform.platform(),  # WRITTEN BY AI
    }  # WRITTEN BY AI


def write_json(path: Path, payload: Dict[str, object]) -> None:  # WRITTEN BY AI
    """Write a JSON payload to disk with stable formatting."""  # WRITTEN BY AI
    path.parent.mkdir(parents=True, exist_ok=True)  # WRITTEN BY AI
    with path.open("w", encoding="utf-8") as handle:  # WRITTEN BY AI
        json.dump(payload, handle, indent=2, sort_keys=True)  # WRITTEN BY AI
