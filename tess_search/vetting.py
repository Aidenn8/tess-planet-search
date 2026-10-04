"""Vetting: tests that separate planet-like signals from impostors.

Most periodic dips found by a transit search are not planets. The usual
impostors, and the test aimed at each:

  eclipsing binary (two stars)      odd/even depth mismatch, secondary eclipse,
                                    V-shaped dip, companion too big to be a planet
  signal from a neighbouring star   centroid shift during transit, duration too
                                    long for the target star
  instrument / data artefacts       one or two transits dominate, transits pile up at
                                    data gaps, period matches the spacecraft orbit,
                                    background flux changes in transit, signal not
                                    present in both halves of the data
  stellar variability               period matches the star's rotation, signal not
                                    distinct from other dips in the folded light curve

Noise is the crux. TESS light curves, centroids and backgrounds have correlated
("red") noise: slow wobbles from pointing jitter, thruster firings and scattered
light. Error bars computed as if each 2-minute point were independent come out
far too small, and every test then cries wolf (a first version of this module
flagged TOI-700 d as a 20-sigma centroid shift). So:
  * flux tests use the scatter measured at the transit's own timescale
    (`timescale_noise`), which includes the red noise;
  * centroid and background tests are calibrated empirically: the same statistic
    is computed at ~24 fake transit epochs (random phases), and the real one must
    stand out from that null distribution (`null_calibrated`).

The tests mirror ones used by the Kepler Robovetter (Thompson et al. 2018) and
TESS vetters such as LEO-Vetter (Kunimoto et al. 2025), simplified and documented.
"""
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

from .inject import expected_duration
from .lightcurve import robust_std

TESS_PIXEL_ARCSEC = 21.0
ORBIT_DAYS = 13.7  # TESS orbital period; scattered light and pointing repeat on it
N_NULL = 24


def phase_offset(time, period, t0):
    """Time from nearest transit centre, in days (between -P/2 and P/2)."""
    return ((time - t0 + 0.5 * period) % period) - 0.5 * period


def epochs(time, period, t0):
    return np.round((time - t0) / period).astype(np.int64)


def timescale_noise(time, flux, window, exclude=None, cadence=2 / 1440):
    """Robust scatter of `window`-day averages of the flux (includes red noise)."""
    keep = np.ones(len(time), bool) if exclude is None else ~exclude
    if keep.sum() < 1000:  # "transit" covers nearly the whole orbit (junk); use everything
        keep = np.ones(len(time), bool)
    t, f = time[keep], flux[keep]
    idx = np.floor((t - t[0]) / window).astype(np.int64)
    uniq, inv, counts = np.unique(idx, return_inverse=True, return_counts=True)
    means = np.bincount(inv, weights=f) / counts
    full = counts >= 0.6 * window / cadence
    if full.sum() < 20:
        return robust_std(f - 1) / np.sqrt(max(window / cadence, 1))
    return robust_std(means[full] - np.median(means[full]))


def weighted_mean(x, err):
    w = 1 / err ** 2
    m = np.sum(w * x) / np.sum(w)
    return m, 1 / np.sqrt(np.sum(w))


def per_transit_depths(lc, period, t0, duration, sigma_core, min_points=5):
    """Depth of each individual transit against its local out-of-transit level.

    sigma_core is the red-noise-aware scatter of a mean over the transit core
    (0.8 x duration); partial transits get proportionally larger errors.
    """
    dt = phase_offset(lc.time, period, t0)
    ep = epochs(lc.time, period, t0)
    core = np.abs(dt) < 0.4 * duration
    local = (np.abs(dt) > 0.75 * duration) & (np.abs(dt) < 3 * duration)
    n_expected = 0.8 * duration / (2 / 1440)
    rows = []
    for e in np.unique(ep[core]):
        mi = core & (ep == e)
        mo = local & (ep == e)
        if mi.sum() < min_points or mo.sum() < 2 * min_points:
            continue
        depth = np.mean(lc.flux[mo]) - np.mean(lc.flux[mi])
        s_in = sigma_core * np.sqrt(max(n_expected / mi.sum(), 1.0))
        s_out = sigma_core * np.sqrt(mi.sum() / mo.sum())
        rows.append((int(e), float(np.median(lc.time[mi])), depth, float(np.hypot(s_in, s_out)), int(mi.sum())))
    return np.array(rows, dtype=[("epoch", "i8"), ("time", "f8"), ("depth", "f8"),
                                  ("err", "f8"), ("n", "i8")])


def phase_snr_curve(lc, period, t0, duration, bins_per_duration=4):
    """Box SNR at every phase of the folded light curve, same width as the transit.

    Returns (phase_days, snr) using white-noise errors. Away from the transit the
    values should scatter with standard deviation 1; red noise makes it larger,
    and that measured scatter (`red_noise_factor`) rescales every phase-based test.
    """
    dt = phase_offset(lc.time, period, t0)
    width = duration / bins_per_duration
    nb = max(int(np.ceil(period / width)), 8)
    idx = np.clip(((dt + period / 2) / period * nb).astype(int), 0, nb - 1)
    s = np.bincount(idx, weights=lc.flux - 1, minlength=nb)
    n = np.bincount(idx, minlength=nb).astype(float)
    k = bins_per_duration
    s_box = np.convolve(np.concatenate([s, s[:k - 1]]), np.ones(k), "valid")[:nb]
    n_box = np.convolve(np.concatenate([n, n[:k - 1]]), np.ones(k), "valid")[:nb]
    sigma = robust_std(lc.flux - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        snr = -s_box / n_box / (sigma / np.sqrt(n_box))
    centres = (np.arange(nb) + k / 2) / nb * period - period / 2
    centres = (centres + period / 2) % period - period / 2
    return centres, np.nan_to_num(snr)


def trapezoid(t, depth, t14, frac, tc):
    """Trapezoid dip: total duration t14, ingress = frac * t14 (frac in 0..0.5)."""
    t12 = frac * t14
    x = np.abs(t - tc)
    y = np.zeros_like(t)
    flat = x <= t14 / 2 - t12
    slope = (x > t14 / 2 - t12) & (x < t14 / 2)
    y[flat] = depth
    y[slope] = depth * (t14 / 2 - x[slope]) / max(t12, 1e-6)
    return 1 - y


def fit_trapezoid(lc, period, t0, duration, depth):
    dt = phase_offset(lc.time, period, t0)
    m = np.abs(dt) < 2.5 * duration
    t, f = dt[m], lc.flux[m]
    if len(t) < 30:
        return None
    nb = 75
    edges = np.linspace(-2.5 * duration, 2.5 * duration, nb + 1)
    idx = np.clip(np.digitize(t, edges) - 1, 0, nb - 1)
    cnt = np.bincount(idx, minlength=nb)
    ok = cnt > 0
    fb = np.bincount(idx, weights=f, minlength=nb)[ok] / cnt[ok]
    tb = (0.5 * (edges[1:] + edges[:-1]))[ok]
    eb = robust_std(lc.flux - 1) / np.sqrt(cnt[ok])
    try:
        p, cov = curve_fit(trapezoid, tb, fb, p0=[max(depth, 1e-5), duration, 0.2, 0.0],
                           sigma=eb, bounds=([0, duration / 4, 0.0, -duration / 2],
                                             [1, 2.5 * duration, 0.5, duration / 2]),
                           maxfev=20000)
        perr = np.sqrt(np.diag(cov))
    except Exception:
        return None
    return {"depth": float(p[0]), "t14": float(p[1]), "ingress_frac": float(p[2]),
            "tc_offset": float(p[3]), "depth_err": float(perr[0]), "t14_err": float(perr[1]),
            "ingress_frac_err": float(perr[2])}


def _sector_slices(sector):
    """Data are time-sorted and sectors are contiguous; return [(start, stop), ...]."""
    change = np.flatnonzero(np.diff(sector)) + 1
    starts = np.concatenate([[0], change])
    stops = np.concatenate([change, [len(sector)]])
    return list(zip(starts, stops))


def inout_stats(lc, values_list, period, t0, duration, slices):
    """In-transit minus local out-of-transit means, sector by sector.

    Returns (chi2, dof, shift, signed): chi2 sums (difference / white error)^2 over
    sectors and series; shift is the median per-sector shift size, de-biased by
    subtracting the expected noise contribution (sizes of pure noise are always
    positive); signed sums difference/error, so its sign says up or down.
    Each sector is handled on its own (pixel axes rotate between cameras). The
    white-noise errors are wrong in absolute terms, but the same statistics are
    computed at fake epochs, so the calibration absorbs that.
    """
    chi2, dof, signed, shifts = 0.0, 0, 0.0, []
    for a, b in slices:
        dt = phase_offset(lc.time[a:b], period, t0)
        core = np.abs(dt) < 0.4 * duration
        local = (np.abs(dt) > 0.75 * duration) & (np.abs(dt) < 3 * duration)
        if core.sum() < 10 or local.sum() < 20:
            continue
        diffs = []
        for vals in values_list:
            v = vals[a:b]
            ci, co = v[core], v[local]
            ci, co = ci[np.isfinite(ci)], co[np.isfinite(co)]
            if len(ci) < 10 or len(co) < 20:
                diffs = None
                break
            noise = robust_std(co - np.median(co))
            if not np.isfinite(noise) or noise <= 0:
                diffs = None
                break
            d = np.mean(ci) - np.mean(co)
            err = noise * np.sqrt(1 / len(ci) + 1 / len(co))
            diffs.append((d, err))
        if not diffs:
            continue
        for d, err in diffs:
            chi2 += (d / err) ** 2
            signed += d / err
            dof += 1
        size2 = sum(d ** 2 - err ** 2 for d, err in diffs)
        shifts.append(np.sqrt(max(size2, 0.0)))
    shift = float(np.median(shifts)) if shifts else np.nan
    return chi2, dof, shift, signed


def null_calibrated(lc, values_list, period, t0, duration, seed=0):
    """Compare the real in/out statistics with the same statistics at fake epochs.

    Returns a dict:
      z, p        how far the real chi-square sits above the fake-epoch chi-squares
                  (robust standard deviations) and the fraction of fakes at least as large
      shift       real de-biased shift size; shift_z its significance vs the fakes
      shift_excess real shift minus the typical fake shift (the part not explained by noise)
      signed_z    signed version (positive = value goes up during transit)
    """
    nan = {"z": np.nan, "p": np.nan, "shift": np.nan, "shift_z": np.nan, "shift_excess": np.nan,
           "signed_z": np.nan}
    slices = _sector_slices(lc.sector)
    real, dof, shift, signed = inout_stats(lc, values_list, period, t0, duration, slices)
    if dof == 0:
        return nan
    rng = np.random.default_rng(seed)
    excl = min(3 * duration / period, 0.15)  # stay clear of the transit and of phase 0.5
    nulls = []
    tries = 0
    while len(nulls) < N_NULL and tries < 4 * N_NULL:
        tries += 1
        ph = rng.uniform(excl, 1 - excl)
        if abs(ph - 0.5) < excl:
            continue
        c, d, s, sg = inout_stats(lc, values_list, period, t0 + ph * period, duration, slices)
        if d == dof:
            nulls.append((c, s, sg))
    if len(nulls) < 8:
        return dict(nan, shift=shift)
    nulls = np.array(nulls)

    def z_of(x, arr):
        spread = robust_std(arr)
        return float((x - np.median(arr)) / spread) if spread > 0 else 0.0

    return {
        "z": z_of(real, nulls[:, 0]),
        "p": float((np.sum(nulls[:, 0] >= real) + 1) / (len(nulls) + 1)),
        "shift": shift,
        "shift_z": z_of(shift, nulls[:, 1]),
        "shift_excess": float(max(shift - np.median(nulls[:, 1]), 0.0)),
        "signed_z": z_of(signed, nulls[:, 2]),
    }


def sap_depth(lc, period, t0, duration):
    """Depth measured in SAP flux (before PDC), crowding-corrected, for comparison."""
    dt = phase_offset(lc.time, period, t0)
    core = np.abs(dt) < 0.4 * duration
    local = (np.abs(dt) > 0.75 * duration) & (np.abs(dt) < 3 * duration)
    ep = epochs(lc.time, period, t0)
    depths = []
    for e in np.unique(ep[core]):
        mi, mo = core & (ep == e), local & (ep == e)
        if mi.sum() >= 5 and mo.sum() >= 10:
            depths.append(np.median(lc.sap[mo]) - np.mean(lc.sap[mi]))
    if len(depths) < 2:
        return np.nan, np.nan
    crowd = np.nanmedian([m.get("crowdsap") or np.nan for m in lc.sector_meta]) if lc.sector_meta else 1.0
    crowd = crowd if np.isfinite(crowd) and crowd > 0 else 1.0
    depths = np.array(depths) / crowd
    return float(np.mean(depths)), float(np.std(depths) / np.sqrt(len(depths)))


def edge_fraction(lc, period, t0, duration, gap=0.5, near=0.5):
    """Fraction of observed transits within `near` days of a data gap > `gap` days."""
    dt = phase_offset(lc.time, period, t0)
    gaps = np.flatnonzero(np.diff(lc.time) > gap)
    edges = np.concatenate([lc.time[gaps], lc.time[gaps + 1], lc.time[[0, -1]]])
    tc = np.unique(np.round((lc.time[np.abs(dt) < 0.4 * duration] - t0) / period))
    tc = t0 + tc * period
    if len(tc) == 0:
        return np.nan
    dist = np.min(np.abs(tc[:, None] - edges[None, :]), axis=1)
    return float(np.mean(dist < near))


def harmonic_flags(period, prot, prot_power, var_amp, depth):
    """Near the spacecraft orbit/sector or the star's rotation (or simple multiples)?

    Rotation is only considered when the star is clearly spotted (periodogram
    power > 0.2 and variability several times the transit depth); otherwise the
    "rotation period" can just be the planet's own signal.
    """
    flags = []
    for base, name in ((ORBIT_DAYS, "orbit"), (2 * ORBIT_DAYS, "sector")):
        for k in (1 / 4, 1 / 3, 1 / 2, 1, 2, 3):
            if abs(period / (base * k) - 1) < 0.01:
                flags.append(f"period near {name} x {k:g}")
    if np.isfinite(prot) and prot_power > 0.2 and var_amp > max(0.002, 3 * depth):
        for k in (0.5, 1, 2):
            if abs(period / (prot * k) - 1) < 0.03:
                flags.append(f"period near rotation x {k:g} (Prot={prot:.2f} d)")
    return flags


def vet(lc, det, r_star, m_star):
    """Compute every metric for one detection. Returns a dict."""
    P, t0, dur = det["period"], det["t0"], det["duration"]
    out = {"period": P, "t0": t0, "duration_h": dur * 24, "depth_ppm": det["depth"] * 1e6,
           "snr": det["snr"], "sde": det["sde"], "n_transits": det["n_transits"]}

    trap = fit_trapezoid(lc, P, t0, dur, det["depth"])
    out["trap"] = trap
    dur_fit = trap["t14"] if trap and trap["t14"] > 0 else dur
    depth = trap["depth"] if trap and trap["depth"] > 0 else det["depth"]
    out["depth_ppm"] = depth * 1e6

    in_any = np.abs(phase_offset(lc.time, P, t0)) < dur_fit
    sigma_core = timescale_noise(lc.time, lc.flux, 0.8 * dur_fit, exclude=in_any)
    white_core = robust_std(lc.flux - 1) / np.sqrt(0.8 * dur_fit / (2 / 1440))
    out["red_noise_beta"] = float(sigma_core / white_core) if white_core > 0 else np.nan

    pt = per_transit_depths(lc, P, t0, dur_fit, sigma_core)
    out["n_transits_measured"] = int(len(pt))
    if len(pt) >= 2:
        mean, merr = weighted_mean(pt["depth"], pt["err"])
        out["snr_red"] = float(mean / merr)  # transit SNR with red-noise-aware errors
        snr_i = pt["depth"] / pt["err"]
        good_snr2 = np.clip(snr_i, 0, None) ** 2
        out["max_single_frac"] = float(good_snr2.max() / good_snr2.sum()) if good_snr2.sum() > 0 else 1.0
        out["chi2_depth"] = float(np.sum(((pt["depth"] - mean) / pt["err"]) ** 2) / (len(pt) - 1))
        excess_var = np.var(pt["depth"]) - np.mean(pt["err"] ** 2)
        out["depth_scatter_frac"] = float(np.sqrt(max(excess_var, 0)) / abs(mean)) if mean != 0 else np.nan
        out["frac_positive"] = float(np.mean(pt["depth"] > 0))
        odd, even = pt["epoch"] % 2 == 1, pt["epoch"] % 2 == 0
        if odd.sum() >= 1 and even.sum() >= 1:
            do, eo = weighted_mean(pt["depth"][odd], pt["err"][odd])
            de, ee = weighted_mean(pt["depth"][even], pt["err"][even])
            out["depth_odd_ppm"], out["depth_even_ppm"] = do * 1e6, de * 1e6
            out["oddeven_sigma"] = float(abs(do - de) / np.hypot(eo, ee))
        half = len(pt) // 2
        if half >= 2:
            d1, e1 = weighted_mean(pt["depth"][:half], pt["err"][:half])
            d2, e2 = weighted_mean(pt["depth"][half:], pt["err"][half:])
            out["halves_snr"] = [float(d1 / e1), float(d2 / e2)]
            out["halves_diff_sigma"] = float(abs(d1 - d2) / np.hypot(e1, e2))
        out["per_transit"] = {"time": pt["time"].tolist(), "depth": pt["depth"].tolist(),
                              "err": pt["err"].tolist()}

    phase, snr_curve = phase_snr_curve(lc, P, t0, dur_fit)
    near_primary = np.abs(phase) < 1.5 * dur_fit
    other = snr_curve[~near_primary]
    out["primary_snr_fold"] = float(snr_curve[near_primary].max()) if near_primary.any() else np.nan
    if other.size > 10:
        red = max(robust_std(other), 1.0)
        out["red_noise_factor"] = float(red)
        out["red_snr"] = float(out["primary_snr_fold"] / red)
        j = int(np.argmax(np.where(~near_primary, snr_curve, -np.inf)))
        out["sec_snr"] = float(snr_curve[j] / red)
        out["sec_phase"] = float(phase[j] / P % 1)
        dts = phase_offset(lc.time, P, t0 + phase[j])
        inbox = np.abs(dts) < dur_fit / 2
        out["sec_depth_ppm"] = float((1 - np.mean(lc.flux[inbox])) * 1e6) if inbox.sum() > 5 else np.nan
        excl = near_primary | (np.abs(((phase - phase[j]) + P / 2) % P - P / 2) < 1.5 * dur_fit)
        third = snr_curve[~excl]
        out["third_snr"] = float(third.max() / red) if third.size else np.nan
        out["uniqueness"] = float(out["red_snr"] - out["third_snr"])
        dt5 = phase_offset(lc.time, P, t0 + P / 2)
        in5 = np.abs(dt5) < dur_fit / 2
        out5 = (np.abs(dt5) > dur_fit) & (np.abs(dt5) < 3 * dur_fit)
        if in5.sum() > 5 and out5.sum() > 10:
            d5 = np.mean(lc.flux[out5]) - np.mean(lc.flux[in5])
            n_tr5 = max(len(np.unique(epochs(lc.time[in5], P, t0 + P / 2))), 1)
            e5 = sigma_core / np.sqrt(n_tr5)
            out["phase05_depth_ppm"] = float(d5 * 1e6)
            out["phase05_sigma"] = float(d5 / e5)

    out["rp_rearth"] = float(np.sqrt(max(depth, 0)) * r_star * 109.1)
    exp = expected_duration(P, r_star, m_star, b=0.0)
    out["expected_duration_h"] = float(exp * 24)
    out["duration_ratio"] = float(dur_fit / exp) if exp > 0 else np.nan

    c = null_calibrated(lc, [lc.centr1, lc.centr2], P, t0, dur_fit)
    out["centroid_z"], out["centroid_p"] = c["z"], c["p"]
    out["centroid_shift_px"], out["centroid_shift_z"] = c["shift"], c["shift_z"]
    # where the dimming source sits relative to the out-of-transit light centre:
    # a source at distance r carrying the whole dip moves the centroid by ~ depth * r
    out["centroid_offset_arcsec"] = (float(c["shift_excess"] / max(depth, 1e-6) * TESS_PIXEL_ARCSEC)
                                     if np.isfinite(c["shift_excess"]) else np.nan)
    b = null_calibrated(lc, [lc.bkg], P, t0, dur_fit, seed=1)
    out["bkg_z"], out["bkg_p"], out["bkg_signed_z"] = b["z"], b["p"], b["signed_z"]
    sd, se = sap_depth(lc, P, t0, dur_fit)
    out["sap_depth_ppm"], out["sap_depth_err_ppm"] = sd * 1e6, se * 1e6
    out["edge_frac"] = edge_fraction(lc, P, t0, dur_fit)
    meta = lc.star_meta or {}
    out["harmonic_flags"] = harmonic_flags(P, meta.get("prot", np.nan), meta.get("prot_power", 0),
                                           meta.get("var_amp", 0), depth)
    return out


SNR_CANDIDATE = 7.3   # same threshold Kepler used for its Threshold Crossing Events
SNR_WEAK = 6.0


def classify(v):
    """Turn metrics into (verdict, list of reasons)."""
    fp, flag = [], []

    def g(k, d=np.nan):
        x = v.get(k)
        return d if x is None or (isinstance(x, float) and not np.isfinite(x)) else x

    snr = g("snr_red", 0.0)

    # --- eclipsing binary signatures
    if g("oddeven_sigma", 0) > 3:
        fp.append(f"odd/even depths differ ({g('oddeven_sigma'):.1f} sigma): likely eclipsing binary at 2x period")
    if g("phase05_sigma", 0) > 4 and g("phase05_depth_ppm", 0) > 0.1 * v["depth_ppm"]:
        fp.append(f"dip at phase 0.5 ({g('phase05_sigma'):.1f} sigma): secondary eclipse, or true period is P/2")
    # small planets emit far too little light for a visible secondary eclipse, so any
    # clear second dip at a fixed phase points to an eclipsing binary (known planets in
    # the validation set never exceed 3.3 sigma here)
    if g("sec_snr", 0) > 5 and g("sec_depth_ppm", 0) > 0.1 * v["depth_ppm"]:
        fp.append(f"significant second dip at phase {g('sec_phase'):.2f} ({g('sec_snr'):.1f} sigma)")
    if g("rp_rearth", 0) > 20:
        fp.append(f"companion radius {g('rp_rearth'):.0f} R_earth is too large for a planet")
    trap = v.get("trap") or {}
    if trap and snr > 20 and trap.get("ingress_frac", 0) > 0.4:
        fp.append("V-shaped dip (grazing eclipse)")

    # --- off-target source (empirically calibrated). A real on-target transit also
    # nudges the light centre slightly (by depth x the star's distance from that
    # centre, a fraction of a pixel), so a "significant" but tiny shift is expected
    # for strong signals. Reject only when the shift size clearly beats the noise
    # AND implies a source more than ~1.5 pixels (30") from the light centre.
    off = g("centroid_offset_arcsec", 0)
    if g("centroid_shift_z", 0) > 5 and g("centroid_p", 1) <= 1 / (N_NULL + 1) + 1e-9 and off > 30:
        fp.append(f"centroid moves during transit: source ~{off:.0f}\" from the target "
                  f"(shift z={g('centroid_shift_z'):.1f} vs fake epochs)")
    elif g("centroid_shift_z", 0) > 4 and off > 15:
        flag.append(f"possible off-target source ~{off:.0f}\" away (shift z={g('centroid_shift_z'):.1f})")
    # a planet crossing this star cannot take much longer than the central-transit
    # duration (an eccentric orbit stretches it by at most ~1.5-2x for plausible e)
    if g("duration_ratio", 1) > 3.0:
        fp.append(f"dip lasts {g('duration_ratio'):.1f}x longer than a planet could take to cross "
                  f"this star (blend or variability)")
    elif g("duration_ratio", 1) > 2.0:
        flag.append(f"transit {g('duration_ratio'):.1f}x longer than expected for this star")

    # --- artefacts
    if v["n_transits"] < 3 or g("n_transits_measured", 0) < 3:
        fp.append("fewer than 3 observed transits")
    if g("max_single_frac", 0) > 0.5:
        fp.append(f"one transit carries {100 * g('max_single_frac'):.0f}% of the signal")
    if g("edge_frac", 0) > 0.5:
        fp.append(f"{100 * g('edge_frac'):.0f}% of transits sit next to data gaps")
    if g("red_snr", 99) < 5:
        fp.append(f"not distinct from other dips in folded light curve (red-noise SNR {g('red_snr'):.1f})")
    elif g("red_snr", 99) < 7:
        flag.append(f"modest red-noise SNR {g('red_snr'):.1f}")
    if g("uniqueness", 99) < 2:
        flag.append(f"another dip in the folded light curve is nearly as strong (uniqueness {g('uniqueness'):.1f})")
    # background RISING in transit means the dip could be over-subtracted sky; a
    # background dip is expected for bright stars whose wings leak into sky pixels
    if g("bkg_signed_z", 0) > 6 and g("bkg_p", 1) <= 1 / (N_NULL + 1) + 1e-9:
        flag.append(f"background rises during transit (z={g('bkg_signed_z'):.1f})")
    hs = v.get("halves_snr")
    if hs and min(hs) < 1.5 and snr > 9:
        flag.append(f"signal weak in one half of the data (half SNRs {hs[0]:.1f}, {hs[1]:.1f})")
    # depth consistency: excess scatter matters only if it is a sizeable fraction of
    # the depth (at high SNR, sector-to-sector dilution differences are significant but small)
    if g("chi2_depth", 1) > 3 and g("depth_scatter_frac", 0) > 0.3:
        flag.append(f"transit depths vary by ~{100 * g('depth_scatter_frac'):.0f}% (chi2/dof {g('chi2_depth'):.1f})")
    for f in v.get("harmonic_flags", []):
        (fp if ("orbit" in f or "sector" in f or snr < 15) else flag).append(f)

    if fp:
        return "false positive", fp + flag
    if snr >= SNR_CANDIDATE and g("red_snr", 0) >= 7 and not flag:
        return "candidate", []
    if snr >= SNR_WEAK:
        return "weak candidate", flag or [f"SNR {snr:.1f} below {SNR_CANDIDATE}"]
    return "below threshold", flag or [f"SNR {snr:.1f} below {SNR_WEAK}"]
