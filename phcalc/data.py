"""Representative pKa values at 25 C (rounded; tabulated values differ by about 0.05 between sources)."""

PKA = {
    "acetic": ((4.76,), 0),
    "formic": ((3.75,), 0),
    "benzoic": ((4.20,), 0),
    "hydrofluoric": ((3.17,), 0),
    "hypochlorous": ((7.53,), 0),
    "carbonic": ((6.35, 10.33), 0),
    "oxalic": ((1.25, 4.27), 0),
    "phosphoric": ((2.15, 7.20, 12.35), 0),
    "citric": ((3.13, 4.76, 6.40), 0),
    "ammonium": ((9.25,), 1),
    "tris": ((8.07,), 1),
}


def acid(name: str, conc: float):
    """Build an :class:`~phcalc.equilibria.Acid` from the table, e.g. ``acid("phosphoric", 0.1)``."""
    from .equilibria import Acid

    pka, charge = PKA[name]
    return Acid(pka, conc, charge=charge, name=name)
