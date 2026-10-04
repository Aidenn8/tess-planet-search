"""End-to-end checks of cleaning + search on synthetic stars with known answers."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tess_search import inject, lightcurve, search  # noqa: E402

R_STAR, M_STAR = 0.40, 0.40


def run(raw):
    lc = lightcurve.prepare(0, raw=raw)
    return lc, search.search_star(lc, R_STAR, M_STAR)[0]


def test_period_grid_is_sorted_and_spans_range():
    p = search.period_grid(350, R_STAR, M_STAR)
    assert np.all(np.diff(p) > 0)
    assert p.min() < 0.41 and p.max() > 39
    # denser at short periods, as Ofir sampling requires
    assert np.diff(p)[0] < np.diff(p)[-1]


def test_transit_mask_and_count():
    t = np.arange(0, 100, 2 / 1440)
    m = search.transit_mask(t, 10.0, 5.0, 0.1, factor=1.0)
    assert abs(m.mean() - 0.1 / 10.0) < 1e-3
    assert search.count_transits(t, 10.0, 5.0, 0.1) == 10


def test_injected_planet_is_recovered():
    raw = inject.synthetic_star(seed=1)
    P, t0 = 7.3172, 1330.5
    raw = inject.inject(raw, P, t0, rp_rearth=1.6, r_star=R_STAR, m_star=M_STAR)
    lc, dets = run(raw)
    assert dets, "nothing detected"
    best = dets[0]
    assert abs(best.period - P) / P < 1e-4, best
    phase = ((best.t0 - t0) / P + 0.5) % 1 - 0.5
    assert abs(phase * P) < 0.02
    expected_depth = (1.6 * inject.R_EARTH / (R_STAR * inject.R_SUN)) ** 2
    assert 0.5 * expected_depth < best.depth < 1.3 * expected_depth
    assert best.snr > 10


def test_two_planets_are_both_found():
    raw = inject.synthetic_star(seed=2)
    raw = inject.inject(raw, 3.1416, 1326.2, rp_rearth=1.8, r_star=R_STAR, m_star=M_STAR)
    raw = inject.inject(raw, 11.789, 1331.0, rp_rearth=2.0, r_star=R_STAR, m_star=M_STAR)
    lc, dets = run(raw)
    found = sorted(round(d.period, 2) for d in dets if d.snr > 8)
    assert 3.14 in found and 11.79 in found, found


def test_pure_noise_gives_no_strong_detection():
    raw = inject.synthetic_star(seed=3)
    lc, dets = run(raw)
    assert all(d.snr < 8 for d in dets), [d.as_dict() for d in dets]


def test_flare_is_removed_but_transit_survives():
    raw = inject.synthetic_star(seed=4)
    raw = inject.inject(raw, 5.0, 1327.0, rp_rearth=2.0, r_star=R_STAR, m_star=M_STAR)
    # a big flare: fast rise, exponential decay, right next to a transit
    t = raw["time"]
    t_flare = 1327.0 + 5.0 * 20 + 0.2
    after = t >= t_flare
    shape = np.zeros_like(t)
    shape[after] = 0.05 * np.exp(-(t[after] - t_flare) / 0.01)
    raw["pdcsap"] = (raw["pdcsap"] * (1 + shape)).astype(np.float32)
    lc = lightcurve.prepare(0, raw=raw)
    peak = (lc.time >= t_flare) & (lc.time < t_flare + 0.006)  # first ~9 min of the flare
    assert peak.sum() == 0, "flare peak should be cut"
    in_transit = np.abs(lc.time - (1327.0 + 5.0 * 20)) < 0.02
    assert in_transit.sum() > 20, "transit points next to the flare should be kept"
