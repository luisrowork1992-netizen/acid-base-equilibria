# acid-base-equilibria

![tests](https://github.com/luisrowork1992-netizen/acid-base-equilibria/actions/workflows/tests.yml/badge.svg)

`phcalc` computes the **exact pH and speciation of aqueous acid-base mixtures**. It solves the full charge balance numerically, so it needs no "x is small" approximations. It also simulates titration curves, finds equivalence points and computes buffer capacity.

## Features

- Monoprotic and polyprotic acids and bases, with any charge (for example NH₄⁺ or Tris-H⁺)
- Strong acids and bases, and spectator ions (Na⁺, Cl⁻, …)
- Exact pH from a monotone charge balance solved with Brent's method; the root is unique
- Alpha fractions computed in log space, so they stay stable at any pH
- Full speciation of every species
- `mix`: combines solutions by volume, with dilution
- Titration curves and automatic equivalence-point detection; peaks are refined with a parabola
- Van Slyke buffer capacity β, computed analytically
- A small table of common pKa values: acetic, formic, carbonic, phosphoric, citric, ammonium, Tris, and others

## How a solution is described

You give the total concentration of each acid-base system and the strong ions that came with it:

| Solution | `phcalc` |
|---|---|
| 0.1 M acetic acid | `Solution([acid("acetic", 0.1)])` |
| 0.1 M sodium acetate | `Solution([acid("acetic", 0.1)], cations=0.1)` |
| 0.1 M ammonia | `Solution([acid("ammonium", 0.1)])` |
| 0.1 M NH₄Cl | `Solution([acid("ammonium", 0.1)], anions=0.1)` |
| 0.01 M HCl | `Solution(anions=0.01)` |

## Install

```bash
git clone https://github.com/luisrowork1992-netizen/acid-base-equilibria.git
cd acid-base-equilibria
pip install -e .
```

## Quick start

```python
import numpy as np
import phcalc as ph

print(ph.Solution([ph.acid("acetic", 0.1)]).pH())        # 2.883
print(ph.Solution(anions=1e-8).pH())                     # 6.978, not 8

acetic = ph.Solution([ph.acid("acetic", 0.1)])
naoh = ph.Solution(cations=0.1)
v = np.linspace(0, 40, 801)
curve = ph.titration_curve(acetic, 25.0, naoh, v)
print(ph.equivalence_points(v, curve))                   # [25.0]
```

## Notebook

[`notebooks/ph_and_titrations.ipynb`](notebooks/ph_and_titrations.ipynb) covers five topics:

- exact pH compared with the textbook shortcut formulas
- speciation diagrams for phosphoric and citric acid
- a Sillén log-concentration diagram
- titration curves, including triprotic phosphoric acid
- buffer capacity, analytic compared with numerical

## Tests

```bash
pip install pytest
pytest
```

The tests check the following:

- pure water, and strong acids and bases
- the 10⁻⁸ M HCl trap
- weak acids and bases against the exact quadratic
- equimolar buffer pH = pKa
- the amphiprotic H₂PO₄⁻ formula
- alpha fractions sum to 1 and stay finite at extreme pH
- mass and charge balance
- buffer capacity against the numerical derivative
- the acetic-acid titration: initial pH, half-equivalence pH, equivalence pH 8.72 and the end point
- both phosphoric-acid end points

## Scope and limitations

- This is an ideal solution model: activity coefficients are 1. Above roughly 0.1 M ionic strength, results differ from measured pH by up to about 0.1 units.
- Strong ions are treated as monovalent. Represent a divalent ion as twice the concentration.
- The pKa values are for 25 °C and are rounded.

## References

- D. C. Harris, *Quantitative Chemical Analysis*, 10th ed., W. H. Freeman (2020).
- J. N. Butler, *Ionic Equilibrium: Solubility and pH Calculations*, Wiley (1998).
- D. D. Van Slyke, *J. Biol. Chem.* 52, 525 (1922).

## License

MIT
