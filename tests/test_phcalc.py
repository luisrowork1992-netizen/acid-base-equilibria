import math

import numpy as np
import pytest

import phcalc as ph

KW = 1e-14


def test_pure_water():
    assert ph.Solution().pH() == pytest.approx(7.0, abs=1e-9)


def test_strong_acid_and_base():
    assert ph.Solution(anions=0.01).pH() == pytest.approx(2.0, abs=1e-5)
    assert ph.Solution(cations=0.001).pH() == pytest.approx(11.0, abs=1e-5)


def test_very_dilute_strong_acid_includes_water():
    c = 1e-8
    h = (c + math.sqrt(c * c + 4 * KW)) / 2
    assert ph.Solution(anions=c).pH() == pytest.approx(-math.log10(h), abs=1e-9)
    assert ph.Solution(anions=c).pH() < 7.0


def test_weak_acid_matches_quadratic():
    ka, c = 10**-4.76, 0.1
    h = (-ka + math.sqrt(ka * ka + 4 * ka * c)) / 2
    assert ph.Solution([ph.acid("acetic", c)]).pH() == pytest.approx(-math.log10(h), abs=1e-4)


def test_weak_base_ammonia():
    kb, c = KW / 10**-9.25, 0.1
    oh = (-kb + math.sqrt(kb * kb + 4 * kb * c)) / 2
    assert ph.Solution([ph.acid("ammonium", c)]).pH() == pytest.approx(14 + math.log10(oh), abs=1e-4)


def test_ammonium_chloride_is_acidic():
    assert ph.Solution([ph.acid("ammonium", 0.1)], anions=0.1).pH() == pytest.approx(5.13, abs=0.01)


def test_equimolar_buffer_sits_at_pka():
    sol = ph.Solution([ph.acid("acetic", 0.2)], cations=0.1)  # 0.1 M HA + 0.1 M NaA
    assert sol.pH() == pytest.approx(4.76, abs=0.01)
    assert ph.henderson_hasselbalch(4.76, 0.1, 0.1) == pytest.approx(4.76)


def test_amphiprotic_dihydrogen_phosphate():
    ka1, ka2, c = 10**-2.15, 10**-7.20, 0.1
    approx = -0.5 * math.log10(ka1 * ka2 * c / (ka1 + c))
    assert ph.Solution([ph.acid("phosphoric", c)], cations=c).pH() == pytest.approx(approx, abs=0.02)


@pytest.mark.parametrize("pka", [[4.76], [2.15, 7.20, 12.35], [3.13, 4.76, 6.40]])
def test_alpha_fractions_sum_to_one(pka):
    pH = np.linspace(-2, 16, 181)
    a = ph.alpha_fractions(pka, pH)
    assert a.shape == (len(pka) + 1, pH.size)
    np.testing.assert_allclose(a.sum(axis=0), 1.0)
    for i, p in enumerate(pka):  # adjacent species are equal at pH = pKa_i
        ai = ph.alpha_fractions(pka, p)
        if len(pka) == 1 or min(abs(p - q) for q in pka if q != p) > 3:
            assert ai[i] == pytest.approx(ai[i + 1], rel=2e-3)


def test_alpha_fractions_extreme_ph_no_overflow():
    a = ph.alpha_fractions([2.15, 7.20, 12.35], [-50.0, 80.0])
    assert np.all(np.isfinite(a))
    assert a[0, 0] == pytest.approx(1.0) and a[-1, 1] == pytest.approx(1.0)


def test_speciation_mass_and_charge_balance():
    sol = ph.Solution([ph.acid("phosphoric", 0.05), ph.acid("ammonium", 0.02)], cations=0.06)
    p = sol.pH()
    sp = sol.speciation()
    assert sum(v for k, v in sp.items() if k.startswith("phosphoric")) == pytest.approx(0.05)
    assert abs(sol.charge_balance(p)) < 1e-13


def test_buffer_capacity_matches_numerical_derivative():
    acetate = ph.Solution([ph.acid("acetic", 0.1)], cations=0.04)
    p0 = acetate.pH()
    dc = 1e-6
    p1 = ph.Solution([ph.acid("acetic", 0.1)], cations=0.04 + dc).pH()
    assert acetate.buffer_capacity() == pytest.approx(dc / (p1 - p0), rel=1e-3)


def test_weak_acid_titration():
    acetic = ph.Solution([ph.acid("acetic", 0.1)])
    naoh = ph.Solution(cations=0.1)
    v = np.linspace(0, 40, 801)
    curve = ph.titration_curve(acetic, 25.0, naoh, v)
    assert curve[0] == pytest.approx(2.88, abs=0.01)
    assert np.all(np.diff(curve) > 0)
    # half-equivalence: pH = pKa
    assert curve[np.argmin(abs(v - 12.5))] == pytest.approx(4.76, abs=0.01)
    # equivalence: 0.05 M acetate, pH ~ 8.72
    assert curve[np.argmin(abs(v - 25.0))] == pytest.approx(8.72, abs=0.02)
    eq = ph.equivalence_points(v, curve)
    assert len(eq) == 1 and eq[0] == pytest.approx(25.0, abs=0.05)


def test_phosphoric_titration_shows_two_breaks():
    h3po4 = ph.Solution([ph.acid("phosphoric", 0.1)])
    v = np.linspace(0, 35, 1401)
    curve = ph.titration_curve(h3po4, 10.0, ph.Solution(cations=0.1), v)
    eq = ph.equivalence_points(v, curve)
    assert len(eq) == 2
    np.testing.assert_allclose(eq, [10.0, 20.0], atol=0.15)


def test_mix_merges_same_acid():
    a = ph.Solution([ph.acid("acetic", 0.2)])
    b = ph.Solution([ph.acid("acetic", 0.1)], cations=0.1)
    m = ph.mix((a, 1.0), (b, 1.0))
    assert len(m.acids) == 1
    assert m.acids[0].conc == pytest.approx(0.15)
    assert m.cations == pytest.approx(0.05)


def test_invalid_inputs():
    with pytest.raises(ValueError):
        ph.Acid([], 0.1)
    with pytest.raises(ValueError):
        ph.Acid(4.76, -1)
    with pytest.raises(ValueError):
        ph.mix((ph.Solution(), 0.0))
