"""Exact pH of aqueous acid-base mixtures from the charge balance.

A solution is described by the *analytical* (total) concentration of each
acid-base system plus the concentration of "strong" ions that never take part
in proton transfer (Na+, K+, Cl-, NO3-, ...). For example:

* 0.1 M acetic acid         -> ``Solution([Acid(4.76, 0.1)])``
* 0.1 M sodium acetate      -> ``Solution([Acid(4.76, 0.1)], cations=0.1)``
* 0.1 M ammonia             -> ``Solution([Acid(9.25, 0.1, charge=+1)])``
* 0.1 M ammonium chloride   -> ``Solution([Acid(9.25, 0.1, charge=+1)], anions=0.1)``
* 0.01 M HCl                -> ``Solution(anions=0.01)``

Activities are taken equal to concentrations (ideal dilute solution, 25 C).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

__all__ = ["Acid", "Solution", "alpha_fractions", "mix", "henderson_hasselbalch"]

LN10 = np.log(10.0)


@dataclass(frozen=True)
class Acid:
    """An acid-base system ``H_nA^charge``, given by its ``pKa`` values (one per proton).

    ``charge`` is the charge of the fully protonated form: 0 for acetic or
    phosphoric acid, +1 for ammonium, and so on.
    """

    pKa: tuple
    conc: float
    charge: int = 0
    name: str = ""

    def __init__(self, pKa, conc, charge=0, name=""):
        pk = tuple(float(p) for p in np.atleast_1d(pKa))
        if not pk:
            raise ValueError("need at least one pKa")
        if conc < 0:
            raise ValueError("concentration must be non-negative")
        object.__setattr__(self, "pKa", pk)
        object.__setattr__(self, "conc", float(conc))
        object.__setattr__(self, "charge", int(charge))
        object.__setattr__(self, "name", name)

    @property
    def n_protons(self) -> int:
        return len(self.pKa)

    def charges(self) -> np.ndarray:
        """Charge of each species, from fully protonated (index 0) to fully deprotonated."""
        return self.charge - np.arange(self.n_protons + 1)

    def alpha(self, pH):
        return alpha_fractions(self.pKa, pH)


def alpha_fractions(pKa, pH):
    """Mole fractions of ``H_nA, H_(n-1)A, ..., A`` at the given pH (scalar or array).

    Computed in log space, so very large or very small pH values do not overflow.
    """
    pKa = np.atleast_1d(np.asarray(pKa, dtype=float))
    pH = np.asarray(pH, dtype=float)
    n = pKa.size
    i = np.arange(n + 1)
    cum = np.concatenate([[0.0], np.cumsum(pKa)])  # sum of the first i pKa values
    # log10 term_i = -(n - i) pH - sum_{j<i} pKa_j   (relative to a common factor)
    logt = (-(n - i)[:, None] * pH.reshape(1, -1) - cum[:, None]) * LN10
    a = np.exp(logt - logsumexp(logt, axis=0))
    return a[:, 0] if pH.ndim == 0 else a.reshape((n + 1,) + pH.shape)


@dataclass
class Solution:
    acids: list = field(default_factory=list)
    cations: float = 0.0   # strong monovalent cations, mol/L (e.g. Na+)
    anions: float = 0.0    # strong monovalent anions, mol/L (e.g. Cl-)
    pKw: float = 14.0

    def charge_balance(self, pH: float) -> float:
        """Net charge concentration (mol/L) at a trial pH; zero at the true pH."""
        h = 10.0 ** -pH
        oh = 10.0 ** (pH - self.pKw)
        q = h - oh + self.cations - self.anions
        for a in self.acids:
            q += a.conc * float(np.dot(a.alpha(pH), a.charges()))
        return q

    def pH(self) -> float:
        """Solve the charge balance. The balance falls monotonically with pH, so the root is unique."""
        return float(brentq(self.charge_balance, -3.0, self.pKw + 3.0, xtol=1e-12, rtol=1e-14))

    def speciation(self, pH: float | None = None) -> dict:
        """Concentration of every species at the equilibrium (or a given) pH."""
        pH = self.pH() if pH is None else pH
        out = {"H+": 10.0 ** -pH, "OH-": 10.0 ** (pH - self.pKw)}
        for k, a in enumerate(self.acids):
            label = a.name or f"acid{k}"
            for i, frac in enumerate(a.alpha(pH)):
                out[f"{label}[{i}]"] = a.conc * float(frac)
        return out

    def buffer_capacity(self, pH: float | None = None) -> float:
        """Van Slyke buffer capacity beta = dC_base / dpH (mol/L per pH unit).

        ``beta = ln10 [h + Kw/h + sum_k C_k var_k]``, where ``var_k`` is the variance
        of the number of protons released by system k.
        """
        pH = self.pH() if pH is None else pH
        h = 10.0 ** -pH
        beta = h + 10.0 ** (pH - self.pKw)
        for a in self.acids:
            al = a.alpha(pH)
            i = np.arange(a.n_protons + 1)
            beta += a.conc * (np.dot(al, i**2) - np.dot(al, i) ** 2)
        return float(LN10 * beta)

    def scaled(self, factor: float) -> "Solution":
        return Solution([replace(a, conc=a.conc * factor) for a in self.acids],
                        self.cations * factor, self.anions * factor, self.pKw)


def mix(*parts) -> Solution:
    """Mix ``(Solution, volume)`` pairs. Acids with the same name, pKa and charge are merged."""
    total = sum(v for _, v in parts)
    if total <= 0:
        raise ValueError("total volume must be positive")
    pKw = parts[0][0].pKw
    acids: dict = {}
    cat = an = 0.0
    for sol, v in parts:
        f = v / total
        cat += sol.cations * f
        an += sol.anions * f
        for a in sol.acids:
            key = (a.name, a.pKa, a.charge)
            prev = acids.get(key)
            acids[key] = replace(a, conc=(prev.conc if prev else 0.0) + a.conc * f)
    return Solution(list(acids.values()), cat, an, pKw)


def henderson_hasselbalch(pKa: float, base: float, acid: float) -> float:
    """Approximate buffer pH, ``pKa + log10([A-]/[HA])``."""
    return pKa + np.log10(base / acid)
