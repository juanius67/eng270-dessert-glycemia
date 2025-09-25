"""Plotting helpers."""
from __future__ import annotations

from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def plot_single(t: np.ndarray, y: np.ndarray, title: str, ylabel: str, outpath: Path) -> None:
    outpath = Path(outpath)
    _ensure_parent(outpath)

    fig, ax = plt.subplots()
    ax.plot(t, y, label=title)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


def plot_overlay(
    t: np.ndarray,
    curves: Dict[str, np.ndarray],
    ylabel: str,
    outpath: Path,
) -> None:
    outpath = Path(outpath)
    _ensure_parent(outpath)

    fig, ax = plt.subplots()
    for name, values in curves.items():
        ax.plot(t, values, label=name)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{ylabel} overlay")
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
