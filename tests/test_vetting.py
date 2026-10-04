"""Vetting and crossmatch checks on synthetic signals with known truth."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tess_search import crossmatch, inject, lightcurve, search, vetting  # noqa: E402

R_STAR, M_STAR = 0.40, 0.40


def detect_and_vet(raw):
    lc = lightcurve.prepare(0, raw=raw)
    dets = search.search_star(lc, R_STAR, M_STAR, max_signals=1)[0]
    d = dets[0].as_dict()
    v = vetting.vet(lc, d, R_STAR, M_STAR)
    return d, v, vetting.classify(v)


def test_clean_planet_passes_vetting():
    raw = inject.synthetic_star(seed=11)
    raw = inject.inject(raw, 4.217, 1327.1, rp_rearth=1.8, r_star=R_STAR, m_star=M_STAR)
    d, v, (verdict, reasons) = detect_and_vet(raw)
    assert abs(d["period"] - 4.217) < 1e-3
    assert verdict == "candidate", (verdict, reasons)
    assert v["oddeven_sigma"] < 3
    assert abs(v["duration_ratio"] - 1) < 0.6


def test_eclipsing_binary_with_unequal_eclipses_is_rejected():
    # true binary period 6.0 d with alternating deep/shallow eclipses: the search
    # locks onto 3.0 d and the odd/even test must catch it
    raw = inject.synthetic_star(seed=12)
    P = 6.0
    a = inject.inject(raw, P, 1327.0, rp_rearth=3.0, r_star=R_STAR, m_star=M_STAR)
    raw2 = inject.inject(a, P, 1327.0 + P / 2, rp_rearth=1.8, r_star=R_STAR, m_star=M_STAR)
    d, v, (verdict, reasons) = detect_and_vet(raw2)
    assert verdict == "false positive", (d, reasons)
    assert any("odd/even" in r or "phase 0.5" in r or "second dip" in r for r in reasons), reasons


def test_secondary_eclipse_is_flagged():
    raw = inject.synthetic_star(seed=13)
    P = 5.3
    raw = inject.inject(raw, P, 1327.0, rp_rearth=3.5, r_star=R_STAR, m_star=M_STAR)
    raw = inject.inject(raw, P, 1327.0 + 0.37 * P, rp_rearth=1.6, r_star=R_STAR, m_star=M_STAR)
    d, v, (verdict, reasons) = detect_and_vet(raw)
    assert verdict == "false positive", reasons


def test_period_relations():
    assert crossmatch.period_relation(10.0, 10.02) == "same"
    assert crossmatch.period_relation(5.0, 10.0) == "half"
    assert crossmatch.period_relation(20.01, 10.0) == "double"
    assert crossmatch.period_relation(10.0, 3.3334) == "triple"
    assert crossmatch.period_relation(10.0, 10.5) is None


def test_harmonic_flags():
    assert vetting.harmonic_flags(13.72, np.nan, 0, 0, 1e-3)
    assert not vetting.harmonic_flags(9.3, np.nan, 0, 0, 1e-3)
    assert vetting.harmonic_flags(2.83, 2.81, 0.6, 0.02, 3e-4)
    # a quiet star: "rotation" peak that is really the planet must not be flagged
    assert not vetting.harmonic_flags(2.04, 2.01, 0.15, 0.001, 4e-4)


def test_timescale_noise_matches_white_noise():
    rng = np.random.default_rng(0)
    t = np.arange(0, 200, 2 / 1440)
    f = 1 + rng.normal(0, 1e-3, len(t))
    w = 1.0 / 24
    expected = 1e-3 / np.sqrt(w / (2 / 1440))
    got = vetting.timescale_noise(t, f, w)
    assert abs(got / expected - 1) < 0.15
