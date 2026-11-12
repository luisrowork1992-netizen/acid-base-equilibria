"""Titration curves and equivalence-point detection."""

from __future__ import annotations

import numpy as np

from .equilibria import Solution, mix

__all__ = ["titration_curve", "equivalence_points"]


def titration_curve(analyte: Solution, v_analyte: float, titrant: Solution, volumes) -> np.ndarray:
    """pH after adding each titrant volume (same units as ``v_analyte``), with dilution."""
    volumes = np.asarray(volumes, dtype=float)
    out = np.empty_like(volumes)
    for i, v in enumerate(volumes):
        out[i] = analyte.pH() if v == 0 else mix((analyte, v_analyte), (titrant, v)).pH()
    return out


def equivalence_points(volumes, pH, min_jump: float = 1.0, window: float = 0.05) -> list:
    """Volumes where dpH/dV peaks. Each peak is refined with a parabola through its 3 points.

    A peak counts only if pH changes by at least ``min_jump`` across a window
    centred on it that spans ``window`` (a fraction) of the volume points. This
    filters out tiny bumps in flat regions.
    """
    v = np.asarray(volumes, dtype=float)
    p = np.asarray(pH, dtype=float)
    slope = np.diff(p) / np.diff(v)
    mid = 0.5 * (v[1:] + v[:-1])
    half = max(1, int(round(window * v.size / 2)))
    found = []
    for i in range(1, slope.size - 1):
        if slope[i] >= slope[i - 1] and slope[i] > slope[i + 1]:
            lo, hi = max(i - half, 0), min(i + 1 + half, p.size - 1)
            if abs(p[hi] - p[lo]) < min_jump:
                continue
            y0, y1, y2 = slope[i - 1:i + 2]
            denom = y0 - 2 * y1 + y2
            shift = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
            found.append(float(mid[i] + shift * (mid[i + 1] - mid[i])))
    return found
