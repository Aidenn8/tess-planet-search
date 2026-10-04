"""Periodic transit search: seasonal box least squares (BLS), stacked.

Why seasons: these stars were observed over ~8 years with long gaps. A single
BLS over the whole span needs a period grid so fine (phase must stay aligned for
8 years) that it would take millions of trial periods per star. Instead we run
BLS separately on each season (a run of sectors without long gaps, usually up
to ~1 year), all on one common period grid, and add their log-likelihoods. A
real planet raises the likelihood at its period in every season, so the sum
grows; noise peaks land at different periods in each season and average out.
Each surviving peak is then refit coherently on all the data (`refine`), with
a narrow but very fine period grid, which pins the period down precisely.

Units: BLS "log likelihood" with per-point errors equals 0.5 * SNR^2 of the box
fit, so the stacked value converts to an (incoherent) SNR = sqrt(2 * sum).
"""
from dataclasses import asdict, dataclass

import numpy as np
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter

from .lightcurve import bin_lightcurve, robust_std

# Transit durations searched, in hours. M dwarfs are small, so transits are short:
# 0.5 h (ultra-short orbits around the smallest stars) to 4 h (40-day orbits).
DURATIONS_H = np.array([0.5, 0.75, 1.0, 1.4, 2.0, 2.8, 4.0])
SEARCH_BIN_MIN = 10.0
PERIOD_MIN, PERIOD_MAX = 0.4, 40.0


def period_grid(time_span, r_star, m_star, pmin=PERIOD_MIN, pmax=PERIOD_MAX, oversample=3):
    """Ofir (2014) optimal frequency sampling for a star of given radius/mass (solar units).

    Trial frequencies are spaced so that over `time_span` days the predicted
    transit time never drifts by more than a fraction of the expected duration.
    """
    G = 6.674e-11
    R = max(r_star, 0.1) * 6.957e8
    M = max(m_star, 0.08) * 1.989e30
    span = time_span * 86400.0
    fmin = 2.0 / span if 2.0 / span > 1 / (pmax * 86400) else 1 / (pmax * 86400)
    fmax = 1 / (pmin * 86400)
    A = (2 * np.pi) ** (2 / 3) / np.pi * R / (G * M) ** (1 / 3) / (span * oversample)
    C = fmin ** (1 / 3) - A / 3
    N = int(np.ceil((fmax ** (1 / 3) - fmin ** (1 / 3) + A / 3) * 3 / A))
    x = np.arange(N) + 1
    freqs = (A / 3 * x + C) ** 3  # Hz
    periods = 1 / freqs / 86400.0
    periods = periods[(periods >= pmin) & (periods <= pmax)]
    return np.sort(periods)


def sde(power, kernel=None):
    """Signal detection efficiency: peak height in units of the periodogram's scatter,
    after removing the slow rise of BLS power towards long periods."""
    kernel = kernel or max(101, (len(power) // 300) | 1)
    base = median_filter(power, size=kernel, mode="nearest")
    resid = power - base
    scatter = robust_std(resid)
    return resid / scatter if scatter > 0 else resid * 0


@dataclass
class Detection:
    period: float
    t0: float
    duration: float          # days
    depth: float             # fractional
    depth_err: float
    snr: float               # coherent box SNR on all data
    sde: float               # stacked-periodogram detection efficiency
    snr_stacked: float       # incoherent SNR from stacked seasons
    n_transits: int          # transits with data
    iteration: int

    def as_dict(self):
        return asdict(self)


def stacked_periodogram(lc, periods, mask=None, bin_min=SEARCH_BIN_MIN):
    """Sum of per-season BLS log-likelihoods on a common period grid."""
    durations = DURATIONS_H / 24.0
    total = np.zeros(len(periods))
    keep = np.ones(len(lc.time), bool) if mask is None else ~mask
    for s in np.unique(lc.season):
        m = keep & (lc.season == s)
        if m.sum() < 500:
            continue
        t, f, n = bin_lightcurve(lc.time[m], lc.flux[m], bin_min / 1440.0)
        sigma = robust_std(lc.flux[m] - 1) / np.sqrt(n)  # per-bin error from point scatter
        ok = n >= 2
        t, f, sigma = t[ok], f[ok], sigma[ok]
        span = t.max() - t.min()
        p_ok = periods < span / 2  # a season must cover at least ~2 cycles
        if p_ok.sum() == 0:
            continue
        bls = BoxLeastSquares(t, f, dy=sigma)
        res = bls.power(periods[p_ok], durations, objective="likelihood", oversample=5)
        ll = np.nan_to_num(np.asarray(res.log_likelihood), nan=0.0)
        ll[np.asarray(res.depth) <= 0] = 0.0  # only dips count
        total[p_ok] += ll
    return total


def refine(lc, period, mask=None, rel_width=None):
    """Coherent BLS on all data around `period`: returns best (P, t0, duration, stats)."""
    keep = np.ones(len(lc.time), bool) if mask is None else ~mask
    tb, fb = lc.time[keep], lc.flux[keep]
    sigma = np.full(len(tb), robust_std(fb - 1))
    span = tb.max() - tb.min()
    durations = DURATIONS_H / 24.0
    # coherent peak width in period is ~ P * duration / span; scan +-30 widths around the seasonal peak
    rel_width = rel_width or max(30 * (durations.min() / span), 2e-5)
    dP = period * durations.min() / span / 5
    grid = np.arange(period * (1 - rel_width), period * (1 + rel_width), dP)
    bls = BoxLeastSquares(tb, fb, dy=sigma)
    res = bls.power(grid, durations, objective="likelihood", oversample=10)
    i = int(np.nanargmax(res.log_likelihood))
    P, t0, dur = float(res.period[i]), float(res.transit_time[i]), float(res.duration[i])
    stats = bls.compute_stats(P, dur, t0)
    return P, t0, dur, float(res.depth[i]), float(res.depth_err[i]), stats


def transit_mask(time, period, t0, duration, factor=1.5):
    """True for points within factor*duration/2 of any transit centre."""
    phase = ((time - t0 + 0.5 * period) % period) - 0.5 * period
    return np.abs(phase) < factor * duration / 2


def count_transits(time, period, t0, duration, min_coverage=0.5):
    """Number of individual transits with at least min_coverage of their points observed."""
    epoch = np.round((time - t0) / period).astype(np.int64)
    in_tr = np.abs(time - (t0 + epoch * period)) < duration / 2
    expected = duration / (2.0 / 1440.0)
    counts = np.bincount(epoch[in_tr] - epoch.min())
    return int(np.sum(counts >= min_coverage * expected))


def search_star(lc, r_star, m_star, max_signals=3, sde_min=7.0, snr_min=6.0):
    """Iterative search: find the strongest signal, mask it, repeat.

    Thresholds here are deliberately loose (more false alarms, fewer missed
    planets); vetting decides what survives. Returns (detections, periods, last
    periodogram SDE array) so plots can show the periodogram.
    """
    span = max(np.ptp(lc.time[lc.season == s]) for s in np.unique(lc.season))
    periods = period_grid(span, r_star, m_star)
    mask = np.zeros(len(lc.time), bool)
    detections = []
    first_sde = None
    for it in range(max_signals):
        power = stacked_periodogram(lc, periods, mask=mask)
        s = sde(power)
        if first_sde is None:
            first_sde = s
        i = int(np.argmax(s))
        if s[i] < sde_min:
            break
        P, t0, dur, depth, depth_err, stats = refine(lc, periods[i], mask=mask)
        snr = depth / depth_err if depth_err > 0 else 0.0
        det = Detection(
            period=P, t0=t0, duration=dur, depth=depth, depth_err=depth_err, snr=snr,
            sde=float(s[i]), snr_stacked=float(np.sqrt(2 * max(power[i], 0))),
            n_transits=count_transits(lc.time[~mask], P, t0, dur), iteration=it,
        )
        detections.append(det)
        if snr < snr_min:
            break
        mask |= transit_mask(lc.time, P, t0, dur, factor=2.0)
    return detections, periods, first_sde
