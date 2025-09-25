"""Metrics helpers for glucose/insulin curves."""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np


def peak_and_tpeak(y: np.ndarray, t: np.ndarray) -> Tuple[float, float]:
    """Return the peak value and the time at which it occurs."""

    idx = int(np.argmax(y))
    return float(y[idx]), float(t[idx])


def auc(
    y: np.ndarray,
    t: np.ndarray,
    baseline: Optional[float] = None,
    tmax: Optional[float] = None,
) -> float:
    """Compute a trapezoidal AUC with optional baseline clipping and horizon."""

    if tmax is not None:
        mask = t <= tmax
        y = y[mask]
        t = t[mask]

    if baseline is not None:
        y = np.clip(y - baseline, 0.0, None)

    return float(np.trapz(y, t))
