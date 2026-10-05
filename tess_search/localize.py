"""Where on the sky does a transit come from? Pixel-level source localization.

TESS pixels are 21" across, so the light curve of a star also contains light from
its neighbours. If a neighbour is an eclipsing binary, its eclipses leak into the
target's light curve, diluted to planet-like depths. The light curve alone cannot
tell these apart; the pixels can, because each star spreads its light over the
pixels in a known pattern (the pixel response function, PRF) centred on its own
position.

Method (the same idea as SPOC's difference-image centroiding, made joint across
sectors and fitted with the PRF instead of a centroid):

1. For every sector, average the images taken during transit and the images taken
   just before and after each transit. "Out minus in" is the difference image: it
   shows only the light that disappeared during transit. Per-pixel uncertainties
   come from the scatter between individual transits, so they include red noise.
2. Fit the out-of-transit image with the PRF of every Gaia star (fluxes from
   Gaia photometry, positions moved to the sector's date with proper motions).
   This calibrates where the stars really fall on the pixels in that sector.
3. For trial sky positions on a grid around the target, model every sector's
   difference image as A x PRF(position) + constant, with one amplitude A shared
   by all sectors (a star loses the same light in every sector). The chi-square
   of that fit, summed over sectors, is a map of where the source can be.
4. Report the best position, its offset from the target with an uncertainty, and
   for every Gaia neighbour how strongly the data exclude it as the source.

The test cases (known planets on target, known nearby eclipsing binaries off
target, and synthetic eclipses injected on neighbours) are in scripts/12_localize.py.
"""
import subprocess
from dataclasses import dataclass, field

import numpy as np
from scipy import stats
from scipy.optimize import least_squares

from . import DATA

PRF_DIR = DATA / "prf"
PRF_URL = "https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/"
TESS_ZP_FLUX = 15000.0      # e-/s for a T = 10 star (TESS Instrument Handbook)
PIXEL_ARCSEC = 21.0
GAIA_EPOCH = 2016.0


# ---------------------------------------------------------------- stars


def gaia_stars(ra, dec, radius_arcsec=150.0, gmag_max=19.5):
    """Gaia DR3 stars around (ra, dec) with TESS magnitudes estimated from Gaia colours.

    Source IDs stay int64 end to end (a float conversion would round 19-digit IDs).
    """
    import astropy.units as u
    import pandas as pd
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    v = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "pmRA", "pmDE", "Gmag", "BPmag", "RPmag",
                        "RUWE", "Plx"], row_limit=-1, column_filters={"Gmag": f"<{gmag_max}"})
    res = v.query_region(SkyCoord(ra * u.deg, dec * u.deg), radius=radius_arcsec * u.arcsec,
                         catalog="I/355/gaiadr3")
    if not res:
        return pd.DataFrame()
    t = res[0]
    df = pd.DataFrame({
        "source_id": np.asarray(t["Source"], dtype=np.int64),
        "ra": np.asarray(t["RA_ICRS"], float), "dec": np.asarray(t["DE_ICRS"], float),
        "pmra": np.nan_to_num(np.ma.filled(t["pmRA"], np.nan)),
        "pmdec": np.nan_to_num(np.ma.filled(t["pmDE"], np.nan)),
        "gmag": np.ma.filled(t["Gmag"], np.nan).astype(float),
        "bp_rp": (np.ma.filled(t["BPmag"], np.nan) - np.ma.filled(t["RPmag"], np.nan)).astype(float),
        "ruwe": np.ma.filled(t["RUWE"], np.nan).astype(float),
        "plx": np.ma.filled(t["Plx"], np.nan).astype(float),
    })
    df["tmag"] = [tess_mag(g, c) for g, c in zip(df.gmag, df.bp_rp)]
    c0 = SkyCoord(ra * u.deg, dec * u.deg)
    df["sep_arcsec"] = c0.separation(SkyCoord(df.ra.values * u.deg, df.dec.values * u.deg)).arcsec
    return df.sort_values("tmag").reset_index(drop=True)


def tess_mag(g, bp_rp):
    """TESS magnitude from Gaia G and BP-RP (Stassun et al. 2019, TIC v8 relation)."""
    if not np.isfinite(g):
        return np.nan
    if not np.isfinite(bp_rp):
        return g - 0.430
    c = float(np.clip(bp_rp, -0.2, 3.5))
    return g - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473


def positions_at(stars, year):
    """Gaia positions moved from epoch 2016.0 to `year` with proper motions (mas/yr)."""
    dt = year - GAIA_EPOCH
    dec = stars.dec.values + stars.pmdec.values * dt / 3.6e6
    ra = stars.ra.values + stars.pmra.values * dt / 3.6e6 / np.cos(np.radians(stars.dec.values))
    return ra, dec


def btjd_to_year(btjd):
    return 2000.0 + (btjd + 2457000.0 - 2451545.0) / 365.25


# ---------------------------------------------------------------- PRF


def _prf_dir(sector, cam, ccd):
    start = "start_s0001" if sector < 4 else "start_s0004"
    return PRF_DIR / start, f"cam{cam}_ccd{ccd}"


def ensure_prf_files(sector, cam, ccd):
    """Mirror the SPOC PRF model files for one camera/CCD (25 small FITS files)."""
    import re

    base, sub = _prf_dir(sector, cam, ccd)
    d = base / sub
    if d.exists() and len(list(d.glob("*.fits"))) >= 25:
        return base
    d.mkdir(parents=True, exist_ok=True)
    url = f"{PRF_URL}{base.name}/{sub}/"
    page = subprocess.run(["curl", "-s", "-f", "-L", url], capture_output=True, text=True).stdout
    names = sorted(set(n for n in re.findall(r'href="([^"]+\.fits)"', page) if "phot" not in n))
    for n in names:
        if not (d / n).exists():
            subprocess.run(["curl", "-s", "-f", "-L", "-o", str(d / n), url + n], check=True)
    return base


class PRFModel:
    """SPOC PRF for one sector stamp: `image(x, y)` = normalised light of a point source
    at stamp pixel coordinates (x = column, y = row, integers at pixel centres)."""

    def __init__(self, stamp):
        from PRF import TESS_PRF

        base = ensure_prf_files(stamp.sector, stamp.camera, stamp.ccd)
        ny, nx = stamp.flux.shape[1:]
        col = stamp.col0 + nx / 2 + 1e-3
        row = stamp.row0 + ny / 2 + 1e-3
        self.prf = TESS_PRF(stamp.camera, stamp.ccd, stamp.sector, col, row, localdatadir=str(base))
        self.shape = (ny, nx)
        self._cache = {}

    def image(self, x, y):
        key = (float(x), float(y))  # exact key: rounding would hide the small steps a fit uses for slopes
        img = self._cache.get(key)
        if img is None:
            if not (-7 < x < self.shape[1] + 6 and -7 < y < self.shape[0] + 6):
                img = np.zeros(self.shape)
            else:
                img = self.prf.locate(float(x), float(y), self.shape)
            if len(self._cache) > 2000:
                self._cache.clear()
            self._cache[key] = img
        return img


# ---------------------------------------------------------------- calibration


@dataclass
class Calibrated:
    """One sector, with the star model that says where each Gaia star lands on its pixels."""
    stamp: object
    prf: PRFModel
    year: float
    dxy: tuple              # pixel shift (dx, dy) added to WCS positions, from the image fit
    fit: dict


def star_flux(tmag):
    return TESS_ZP_FLUX * 10 ** (-0.4 * (np.asarray(tmag, float) - 10))


def fit_direct(image, stamp, prf, stars, target_idx, year):
    """Fit a sector's median image with every Gaia star's PRF (fluxes from Gaia photometry):
    free shared pixel shift (dx, dy), flux scale and background."""
    ra, dec = positions_at(stars, year)
    x, y = stamp.sky_to_pix(ra, dec)
    flux = star_flux(stars.tmag.values)
    ny, nx = image.shape
    use = (x > -6) & (x < nx + 5) & (y > -6) & (y < ny + 5) & np.isfinite(flux)
    use &= flux > 1e-3 * flux[target_idx]
    idx = np.flatnonzero(use)
    w = 1 / np.sqrt(np.abs(image) + 25.0)

    def model(p):
        dx, dy, s, b = p
        img = np.full((ny, nx), b)
        for i in idx:
            img = img + s * flux[i] * prf.image(x[i] + dx, y[i] + dy)
        return img

    best = None
    for start in ((0.0, 0.0), (0.3, 0.3), (-0.3, -0.3)):
        r = least_squares(lambda p: ((model(p) - image) * w).ravel(), [start[0], start[1], 1.0, 0.0],
                          bounds=([-1.5, -1.5, 0.2, -np.inf], [1.5, 1.5, 5.0, np.inf]),
                          x_scale=[0.1, 0.1, 0.1, 10], diff_step=1e-3)
        if best is None or r.cost < best.cost:
            best = r
    dx, dy, s, b = best.x
    tgt = prf.image(x[target_idx] + dx, y[target_idx] + dy)
    return {"dx": float(dx), "dy": float(dy), "scale": float(s), "bkg": float(b),
            "target_flux_total": float(s * flux[target_idx]),
            "target_frac_in_aperture": float(tgt[stamp.aperture].sum()),
            "resid_frac": float(np.sum(np.abs(model(best.x) - image)) / np.sum(np.abs(image)))}


def remove_flares(stamp):
    """Drop flare cadences (same detector as the light-curve pipeline, run on the summed
    aperture flux). M dwarfs flare often; a flare inside an in-transit window makes a
    strongly negative difference image that no dimming source can explain."""
    import copy

    from .lightcurve import flare_mask

    lc = stamp.flux[:, stamp.aperture].sum(1)
    med = np.median(lc)
    if not np.isfinite(med) or med <= 0:
        return stamp
    bad = flare_mask(stamp.time, lc / med)
    out = copy.copy(stamp)
    out.time, out.flux = stamp.time[~bad], stamp.flux[~bad]
    return out


def calibrate(stamps, stars, target_idx, max_resid=0.5):
    cals = []
    for st in stamps:
        st = remove_flares(st)
        prf = PRFModel(st)
        year = btjd_to_year(np.median(st.time))
        fit = fit_direct(np.median(st.flux, 0), st, prf, stars, target_idx, year)
        fit["sector"] = st.sector
        if fit["resid_frac"] < max_resid:  # skip stamps the star model cannot describe
            cals.append(Calibrated(stamp=st, prf=prf, year=year, dxy=(fit["dx"], fit["dy"]), fit=fit))
    return cals


def pix(cal, ra, dec):
    x, y = cal.stamp.sky_to_pix(ra, dec)
    return x + cal.dxy[0], y + cal.dxy[1]


# ---------------------------------------------------------------- difference images


def _event_diff(stamp, tc, t14, w_out, n_core_exp):
    """Out-of-transit minus in-transit image for one event centred at tc, or None if poorly covered."""
    dt = stamp.time - tc
    mi = np.abs(dt) < 0.4 * t14
    mb = (dt < -0.75 * t14) & (dt > -0.75 * t14 - w_out)
    ma = (dt > 0.75 * t14) & (dt < 0.75 * t14 + w_out)
    if mi.sum() < max(3, 0.5 * n_core_exp) or mb.sum() < 3 or ma.sum() < 3:
        return None
    out = 0.5 * (stamp.flux[mb].mean(0) + stamp.flux[ma].mean(0))  # average cancels linear drifts
    return out - stamp.flux[mi].mean(0)


def _clip_events(diffs, aperture, nsig=5.0):
    """Drop events whose summed in-aperture change is a gross outlier among the events
    (a leftover flare or a glitch); robust to the signal itself since all events share it."""
    if len(diffs) < 5:
        return diffs
    tot = diffs[:, aperture].sum(1)
    med = np.median(tot)
    mad = 1.4826 * np.median(np.abs(tot - med))
    if mad <= 0:
        return diffs
    return diffs[np.abs(tot - med) < nsig * mad]


def difference_image(stamp, period, t0, t14, min_transits=1, n_null=60, seed=0):
    """Out-of-transit minus in-transit image for one sector (e-/s; positive where light was lost).

    Each transit is compared with windows just before and after it. Per-pixel
    errors come from "null events": the same measurement at ~60 random times in
    the sector away from any transit. Their scatter is the noise of a single
    event at the transit's own timescale (red noise included), and works for
    long periods with only one or two transits per sector.
    Returns (diff, err, n_transits) or None.
    """
    w_out = max(t14, 1.0 / 24)
    n_core_exp = 0.8 * t14 / (2 / 1440)
    t = stamp.time
    epochs = np.unique(np.round((t - t0) / period).astype(np.int64))
    diffs = [d for d in (_event_diff(stamp, t0 + e * period, t14, w_out, n_core_exp) for e in epochs)
             if d is not None]
    if len(diffs) < min_transits:
        return None
    diffs = _clip_events(np.array(diffs), stamp.aperture)
    if len(diffs) < min_transits:
        return None
    # null events: centres spaced by a full window; no part of a null window may touch a real
    # transit (|offset| > window half-width + half the transit duration)
    half = 0.75 * t14 + w_out
    rng = np.random.default_rng(seed + int(stamp.sector))
    centres = np.arange(t.min() + half, t.max() - half, 2 * half) + rng.uniform(0, half)
    ph = ((centres - t0 + 0.5 * period) % period) - 0.5 * period
    centres = centres[np.abs(ph) > half + 0.5 * t14]
    rng.shuffle(centres)
    nulls = []
    for c in centres:
        d = _event_diff(stamp, c, t14, w_out, n_core_exp)
        if d is not None:
            nulls.append(d)
        if len(nulls) >= n_null:
            break
    if len(nulls) < 10:
        return None
    nulls = _clip_events(np.array(nulls), stamp.aperture)
    sd1 = 1.4826 * np.median(np.abs(nulls - np.median(nulls, 0)), 0)
    sd1 = np.maximum(np.maximum(sd1, nulls.std(0, ddof=1) * 0.5), 0.3 * np.median(sd1))
    n = len(diffs)
    return np.mean(diffs, 0), sd1 / np.sqrt(n), n


# ---------------------------------------------------------------- localization


def _fit(items, xy):
    """chi2 of every sector's difference image modelled as A x PRF(x_s, y_s) + c_s, one A shared
    by all sectors (c_s absorbs a per-sector background change). xy = [(x_s, y_s), ...]."""
    num = den = 0.0
    terms = []
    for (cal, diff, err, _), (x, y) in zip(items, xy):
        p = cal.prf.image(x, y).ravel()
        d, w = diff.ravel(), 1 / err.ravel() ** 2
        W = w.sum()
        pc, dc = p - (w * p).sum() / W, d - (w * d).sum() / W
        num += (w * pc * dc).sum()
        den += (w * pc * pc).sum()
        terms.append((pc, dc, w))
    A = num / den if den > 0 else 0.0
    chi2 = sum(((dc - A * pc) ** 2 * w).sum() for pc, dc, w in terms)
    return float(chi2), float(A), float(1 / np.sqrt(den)) if den > 0 else np.inf


def sigma_from_dchi2(dchi2, dof=2):
    """Equivalent two-sided Gaussian sigma for a chi-square increase with `dof` parameters."""
    p = stats.chi2.sf(max(dchi2, 0.0), dof)
    return float(stats.norm.isf(p / 2)) if p > 1e-300 else 37.0


@dataclass
class Localization:
    tic: int
    offset_arcsec: float          # best-fitting source position minus target position
    offset_east: float            # east and north components (arcsec)
    offset_north: float
    offset_err_arcsec: float      # 1-sigma radius from the chi-square surface
    target_sigma: float           # how strongly the target is excluded as the source
    chi2_red: float
    amplitude: float              # e-/s lost at the best position
    amplitude_err: float
    expected_amplitude: float     # depth x target flux: what the target would lose if it hosts the signal
    n_sectors: int
    n_transits: int
    stars: list
    grid: dict


def localize(tic, cals, stars, target_idx, period, t0, t14, depth, half_width=100.0, step=5.0,
             fine_step=1.0):
    """Localize one periodic signal. Positions are offsets (arcsec east, north) from the target's
    position on each sector's date, so a fast-moving target stays at (0, 0)."""
    items = []
    for cal in cals:
        r = difference_image(cal.stamp, period, t0, t14)
        if r is not None:
            items.append((cal, *r))
    if not items:
        return None
    tgt = stars.iloc[[target_idx]]
    tpos = [positions_at(tgt, cal.year) for cal, *_ in items]
    cosd = np.cos(np.radians(stars.dec.values[target_idx]))

    def xy_at(east, north):
        out = []
        for (cal, *_), (ra, dec) in zip(items, tpos):
            out.append(tuple(v[0] for v in pix(cal, ra[0] + east / 3600 / cosd, dec[0] + north / 3600)))
        return out

    def scan(ce, cn, hw, st):
        offs = np.arange(-hw, hw + 1e-9, st)
        chi = np.empty((len(offs), len(offs)))
        for i, on in enumerate(offs):
            for j, oe in enumerate(offs):
                chi[i, j] = _fit(items, xy_at(ce + oe, cn + on))[0]
        return offs, chi

    offs, chi = scan(0.0, 0.0, half_width, step)
    i, j = np.unravel_index(np.argmin(chi), chi.shape)
    ce, cn = offs[j], offs[i]
    foffs, fchi = scan(ce, cn, 2 * step, fine_step)
    fi, fj = np.unravel_index(np.argmin(fchi), fchi.shape)
    be, bn = ce + foffs[fj], cn + foffs[fi]
    chi2_min, A, A_err = _fit(items, xy_at(be, bn))
    dof = sum(d.size for _, d, _, _ in items) - len(items) - 3
    chi2_red = chi2_min / dof
    scale = max(chi2_red, 1.0)  # inflate if the noise model underestimates the scatter
    inside = (fchi - chi2_min) / scale < 2.30   # 1-sigma region for 2 parameters
    err = max(float(np.sqrt(inside.sum() * fine_step ** 2 / np.pi)), fine_step / 2)

    rows = []
    med_scale = float(np.median([cal.fit["scale"] for cal, *_ in items]))
    for k in range(len(stars)):
        xy = [tuple(v[0] for v in pix(cal, *positions_at(stars.iloc[[k]], cal.year))) for cal, *_ in items]
        sep = float(stars.sep_arcsec.values[k])
        if sep > half_width and k != target_idx:
            continue
        c2, a, _ = _fit(items, xy)
        rows.append({"k": k, "c2": c2, "a": a, "sep": sep})
    # a star's own (proper-motion) track can fit better than any fixed grid offset
    chi2_ref = min([chi2_min] + [r["c2"] for r in rows])
    star_rows, rows = rows, []
    ymid = float(np.median([cal.year for cal, *_ in items]))
    sra, sdec = positions_at(stars, ymid)
    tra, tdec = positions_at(stars.iloc[[target_idx]], ymid)
    for r0 in star_rows:
        k, c2, a, sep = r0["k"], r0["c2"], r0["a"], r0["sep"]
        dchi2 = (c2 - chi2_ref) / scale
        f = star_flux(stars.tmag.values[k]) * med_scale
        rows.append({"source_id": int(stars.source_id.values[k]), "sep_arcsec": sep,
                     "east": float((sra[k] - tra[0]) * cosd * 3600), "north": float((sdec[k] - tdec[0]) * 3600),
                     "tmag": float(stars.tmag.values[k]), "is_target": k == target_idx,
                     "dchi2": float(dchi2), "excluded_sigma": sigma_from_dchi2(dchi2),
                     "amplitude": float(a), "implied_eclipse_depth": float(a / f) if f > 0 else np.nan})
    t_row = next(r for r in rows if r["is_target"])
    expected = depth * float(np.median([cal.fit["target_flux_total"] for cal, *_ in items]))
    loc = Localization(
        tic=int(tic), offset_arcsec=float(np.hypot(be, bn)), offset_east=float(be), offset_north=float(bn),
        offset_err_arcsec=err, target_sigma=t_row["excluded_sigma"], chi2_red=float(chi2_red),
        amplitude=A, amplitude_err=A_err * np.sqrt(scale), expected_amplitude=expected,
        n_sectors=len(items), n_transits=int(sum(n for *_, n in items)),
        stars=sorted(rows, key=lambda r: r["dchi2"]),
        grid={"offsets_arcsec": offs.tolist(), "dchi2": ((chi - chi2_ref) / scale).tolist(),
              "sectors": [cal.stamp.sector for cal, *_ in items]},
    )
    return loc, items


R68 = np.sqrt(stats.chi2.ppf(0.6827, 2))   # 68% containment radius of a 2-D Gaussian, in sigma (1.515)


def with_systematics(rec, sys_arcsec):
    """Fold a systematic position error (1-D sigma, arcsec) into a stored localization record.

    The chi-square map gives the statistical error only; the PRF model and the
    sky-to-pixel calibration add a floor of a few arcsec, measured on sources whose
    true position is known (scripts/12_localize.py). For a star at distance d from
    the best position, dchi2 = d^2 / s_stat^2, so with the floor it becomes
    dchi2 * s_stat^2 / (s_stat^2 + s_sys^2).
    """
    s_stat = rec["offset_err_arcsec"] / R68
    shrink = s_stat ** 2 / (s_stat ** 2 + sys_arcsec ** 2)
    out = dict(rec)
    out["offset_err_total_arcsec"] = float(R68 * np.hypot(s_stat, sys_arcsec))
    out["stars"] = [dict(r, excluded_sigma_total=sigma_from_dchi2(r["dchi2"] * shrink)) for r in rec["stars"]]
    t = next(r for r in out["stars"] if r["is_target"])
    out["target_sigma_total"] = t["excluded_sigma_total"]
    return out


def inject_eclipse(cals, stars, k, period, t0, t14, amplitude):
    """Copies of the calibrated sectors with a box eclipse of `amplitude` e-/s removed from
    star k's PRF during each event (a synthetic eclipsing binary at that star's position)."""
    import copy

    out = []
    for cal in cals:
        st = copy.copy(cal.stamp)
        st.flux = cal.stamp.flux.copy()
        x, y = pix(cal, *positions_at(stars.iloc[[k]], cal.year))
        p = cal.prf.image(x[0], y[0]).astype(np.float32)
        dt = ((st.time - t0 + 0.5 * period) % period) - 0.5 * period
        st.flux[np.abs(dt) < t14 / 2] -= amplitude * p
        out.append(Calibrated(stamp=st, prf=cal.prf, year=cal.year, dxy=cal.dxy, fit=cal.fit))
    return out


def aperture_fraction(cals, stars, k):
    """Median fraction of star k's light inside SPOC's photometric aperture."""
    fr = []
    for cal in cals:
        x, y = pix(cal, *positions_at(stars.iloc[[k]], cal.year))
        fr.append(float(cal.prf.image(x[0], y[0])[cal.stamp.aperture].sum()))
    return float(np.median(fr))
