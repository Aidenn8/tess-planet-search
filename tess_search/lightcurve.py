"""Turn a star's raw multi-sector data into a flat light curve ready for searching.

Steps (in `prepare`):
1. keep cadences SPOC marked as good (QUALITY == 0) with finite PDCSAP flux
2. normalise each sector to a median of 1
3. remove flares: M dwarfs flare often, and a flare's sharp rise and slow decay
   can drag the detrending curve upward and fake a dip next to it. Points far
   above the local level are cut, together with a short tail after each flare.
4. detrend with a robust sliding biweight filter (wotan) whose window (0.5 d)
   is several times longer than any transit we search for, so it removes
   starspot modulation and instrument drifts but leaves transits intact
5. split the time series into "seasons" (runs of sectors without long gaps);
   the search treats each season separately and then adds them up
"""
from dataclasses import dataclass, field

import numpy as np
from wotan import flatten

from .download import load_star

CADENCE_DAYS = 2.0 / 60 / 24


@dataclass
class LightCurve:
    tic: int
    time: np.ndarray       # BTJD (BJD - 2457000), days
    flux: np.ndarray       # detrended, normalised (1 = out of transit)
    flux_err: np.ndarray
    raw: np.ndarray        # normalised PDCSAP before detrending
    trend: np.ndarray
    sector: np.ndarray
    centr1: np.ndarray
    centr2: np.ndarray
    bkg: np.ndarray
    sap: np.ndarray        # normalised SAP flux (before PDC corrections)
    season: np.ndarray
    sector_meta: list = field(default_factory=list)
    star_meta: dict = field(default_factory=dict)

    def subset(self, mask):
        keep = {}
        for name, value in self.__dict__.items():
            keep[name] = value[mask] if isinstance(value, np.ndarray) and value.shape == self.time.shape else value
        return LightCurve(**keep)


def robust_std(x):
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def flare_mask(time, flux, window=0.25, nsigma=3.0, tail_cadences=10):
    """True for points to drop: upward outliers (flares) plus a short tail after each.

    The local level comes from a running median, so slow starspot variations
    are not mistaken for flares. Transits are dips, so they are never removed here.
    """
    trend = flatten(time, flux, method="median", window_length=window, return_trend=True)[1]
    resid = flux - trend
    sigma = robust_std(resid)
    high = resid > nsigma * sigma
    # a flare is several bright cadences in a row; a lone 3-sigma point is just noise
    pair = high & (np.roll(high, 1) | np.roll(high, -1))
    bad = pair | (resid > 5 * sigma)
    if bad.any():
        idx = np.flatnonzero(bad)
        for k in range(1, tail_cadences + 1):
            j = idx + k
            j = j[j < len(bad)]
            # only extend into cadences that are actually contiguous in time
            close = time[j] - time[j - k] < (k + 0.5) * CADENCE_DAYS
            bad[j[close]] = True
    return bad


def seasons_from_time(time, gap_days=40.0):
    """Label runs of data separated by gaps longer than gap_days."""
    breaks = np.flatnonzero(np.diff(time) > gap_days)
    season = np.zeros(len(time), dtype=np.int16)
    for b in breaks:
        season[b + 1:] += 1
    return season


def rotation(time, flux, sector):
    """Starspot rotation estimate: Lomb-Scargle on 30-min bins, sector by sector.

    Returns (period_days, median_power, amplitude). Per-sector because PDC can
    distort variability longer than ~13 days and sectors are normalised separately.
    """
    from astropy.timeseries import LombScargle

    periods, powers, amps = [], [], []
    for s in np.unique(sector):
        m = sector == s
        if m.sum() < 2000:
            continue
        t, f, _ = bin_lightcurve(time[m], flux[m], 30 / 1440)
        freq, pw = LombScargle(t, f).autopower(minimum_frequency=1 / 13, maximum_frequency=10,
                                                samples_per_peak=5)
        i = int(np.argmax(pw))
        periods.append(1 / freq[i])
        powers.append(pw[i])
        lo, hi = np.percentile(f, [5, 95])
        amps.append(hi - lo)
    if not periods:
        return np.nan, 0.0, 0.0
    j = int(np.argsort(powers)[len(powers) // 2])  # the sector with median peak power
    return float(periods[j]), float(np.median(powers)), float(np.median(amps))


def detrend_window(prot, power, amp, default=0.5, floor=0.3):
    """Shorter biweight window for fast, strongly spotted rotators, which otherwise
    leave rotation-period residuals that can mimic transits. Never below `floor`
    (0.3 d is still ~2x the longest transit we search)."""
    if np.isfinite(prot) and power > 0.3 and amp > 0.003:
        return float(np.clip(prot / 6, floor, default))
    return default


def prepare(tic, window=None, raw=None):
    """Load, clean and detrend all sectors of a star. `raw` lets tests inject signals."""
    d = load_star(tic) if raw is None else raw
    good = (d["quality"] == 0) & np.isfinite(d["pdcsap"]) & np.isfinite(d["time"])
    good &= np.isfinite(d["pdcsap_err"]) & (d["pdcsap"] > 0)
    time = d["time"][good]
    flux = d["pdcsap"][good].astype(np.float64)
    err = d["pdcsap_err"][good].astype(np.float64)
    sap = d["sap"][good].astype(np.float64)
    sector = d["sector"][good]

    for s in np.unique(sector):
        m = sector == s
        med = np.nanmedian(flux[m])
        flux[m] /= med
        err[m] /= med
        sap[m] /= np.nanmedian(sap[m])

    prot, prot_power, amp = rotation(time, flux, sector)
    if window is None:
        window = detrend_window(prot, prot_power, amp)

    drop = flare_mask(time, flux)
    flare_fraction = float(drop.mean())
    keep = ~drop
    time, flux, err, sap, sector = time[keep], flux[keep], err[keep], sap[keep], sector[keep]
    centr1 = d["centr1"][good][keep]
    centr2 = d["centr2"][good][keep]
    bkg = d["bkg"][good][keep]

    flat, trend = flatten(time, flux, method="biweight", window_length=window,
                          break_tolerance=0.3, edge_cutoff=0.0, return_trend=True)
    ok = np.isfinite(flat)
    lc = LightCurve(
        tic=int(tic), time=time, flux=flat, flux_err=err, raw=flux, trend=trend,
        sector=sector, centr1=centr1, centr2=centr2, bkg=bkg, sap=sap,
        season=seasons_from_time(time),
        sector_meta=d.get("sector_meta", []), star_meta=d.get("star_meta", {}),
    ).subset(ok)
    # final safety clip of wild single points far below (cosmic rays, glitches);
    # 15 sigma is far deeper than any planet we could detect per cadence here
    sigma = robust_std(lc.flux - 1)
    lc = lc.subset(lc.flux > 1 - 15 * sigma)
    lc.star_meta = dict(lc.star_meta, prot=prot, prot_power=prot_power, var_amp=amp,
                        window=window, flare_fraction=flare_fraction,
                        noise_ppm=float(robust_std(lc.flux - 1) * 1e6))
    return lc


def bin_lightcurve(time, flux, bin_days):
    """Average into fixed-width bins; returns bin-centre times, mean flux, counts."""
    idx = np.floor((time - time[0]) / bin_days).astype(np.int64)
    uniq, inv, counts = np.unique(idx, return_inverse=True, return_counts=True)
    fsum = np.bincount(inv, weights=flux)
    tsum = np.bincount(inv, weights=time)
    return tsum / counts, fsum / counts, counts
