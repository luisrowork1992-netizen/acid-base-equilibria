"""phcalc: exact pH, speciation, titration curves and buffer capacity of aqueous solutions."""

from .data import PKA, acid
from .equilibria import Acid, Solution, alpha_fractions, henderson_hasselbalch, mix
from .titration import equivalence_points, titration_curve

__version__ = "0.1.0"
