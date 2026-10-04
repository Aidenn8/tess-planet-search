"""Deeper checks for signals that survive vetting, using Gaia DR3.

TESS pixels are 21" wide, so light from neighbouring stars mixes into the
target's aperture. A dip could come from an eclipsing binary next door whose
eclipse is diluted down to a planet-like depth. For each Gaia neighbour within
~3 pixels we ask: how deep would its own eclipse have to be to produce the
observed dip? If that exceeds 100% the neighbour cannot be the source.

Also reports the target's Gaia RUWE (> 1.4 hints at an unresolved companion)
and Gaia's non-single-star flag.
"""
import numpy as np

SEARCH_RADIUS_ARCSEC = 63.0  # 3 TESS pixels


def fit_transit(lc, period, t0, duration, depth, r_star, m_star, window=3.0):
    """Least-squares fit of a physical transit model (batman, quadratic limb darkening
    fixed at M-dwarf values) to all data within `window` durations of each transit.

    Free parameters: period, mid-transit time, Rp/R*, impact parameter b, a/R*.
    Returns best values, 1-sigma errors (from the covariance, scaled by the reduced
    chi-square), and derived planet radius, transit duration and stellar density.
    """
    import batman
    from scipy.optimize import least_squares

    from .inject import LIMB_DARKENING, a_over_r

    ph = ((lc.time - t0 + 0.5 * period) % period) - 0.5 * period
    m = np.abs(ph) < window * duration
    t, f = lc.time[m], lc.flux[m]
    sigma = np.std(f[np.abs(ph[m]) > duration]) or 1e-3
    # reference epoch near the middle of the data keeps period and t0 uncorrelated
    n_mid = np.round((np.median(t) - t0) / period)
    tref = t0 + n_mid * period
    params = batman.TransitParams()
    params.ecc, params.w, params.u, params.limb_dark = 0.0, 90.0, list(LIMB_DARKENING), "quadratic"
    params.t0, params.per, params.rp, params.a, params.inc = tref, period, np.sqrt(max(depth, 1e-6)), 20.0, 89.0
    model = batman.TransitModel(params, t, supersample_factor=3, exp_time=2 / 1440)

    def curve(x):
        P, tc, k, b, aR = x
        params.per, params.t0, params.rp, params.a = P, tc, k, aR
        params.inc = float(np.degrees(np.arccos(np.clip(b / aR, 0, 1))))
        return model.light_curve(params)

    aR0 = a_over_r(period, r_star, m_star)
    # b and a/R* are degenerate at low SNR; the star's catalogue density pins a/R*
    # (rho ~ (a/R*)^3 / P^2), so a 20% density uncertainty is a ~7% prior on a/R*
    aR_sigma = aR0 * 0.2 / 3
    x0 = [period, tref, np.sqrt(max(depth, 1e-6)), 0.4, aR0]
    lo = [period * (1 - 1e-3), tref - duration, 0.002, 0.0, 1.5]
    hi = [period * (1 + 1e-3), tref + duration, 0.5, 1.2, 200.0]

    def resid(x):
        return np.append((curve(x) - f) / sigma, (x[4] - aR0) / aR_sigma)

    best, best_cost = None, np.inf
    for b0 in (0.2, 0.5, 0.8):  # a few starts in impact parameter
        x0[3] = b0
        r = least_squares(resid, x0, bounds=(lo, hi), x_scale="jac")
        if r.cost < best_cost:
            best, best_cost = r, r.cost
    J = best.jac
    dof = max(len(f) - len(best.x), 1)
    chi2r = 2 * best.cost / dof
    try:
        cov = np.linalg.inv(J.T @ J) * chi2r
        err = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        err = np.full(len(best.x), np.nan)
    P, tc, k, b, aR = best.x
    rho_sun = 1.41  # g/cm^3
    G = 6.674e-8
    rho_star = 3 * np.pi * aR ** 3 / (G * (P * 86400) ** 2)
    t14 = P / np.pi * np.arcsin(np.sqrt(max((1 + k) ** 2 - b ** 2, 0)) / aR / np.sin(np.arccos(min(b / aR, 1))))
    return {
        "period": float(P), "period_err": float(err[0]), "t0_btjd": float(tc), "t0_err": float(err[1]),
        "t0_bjd": float(tc + 2457000), "rp_rs": float(k), "rp_rs_err": float(err[2]), "b": float(b), "b_err": float(err[3]),
        "a_rs": float(aR), "a_rs_err": float(err[4]), "depth_ppm": float(k ** 2 * 1e6),
        "rp_rearth": float(k * r_star * 109.1), "rp_rearth_err": float(err[2] * r_star * 109.1),
        "t14_h": float(t14 * 24), "rho_star_fit_cgs": float(rho_star),
        "rho_star_tic_cgs": float(m_star / r_star ** 3 * rho_sun), "chi2_reduced": float(chi2r), "n_points": int(len(f)),
    }


def gaia_neighbours(ra, dec, radius_arcsec=SEARCH_RADIUS_ARCSEC):
    """Gaia DR3 sources around (ra, dec), via VizieR (I/355/gaiadr3).

    VizieR rather than the ESA archive, which was timing out (Oct 2026) while
    being prepared for Gaia DR4.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    v = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "Gmag", "RPmag", "RUWE", "NSS", "Plx"],
               row_limit=-1)
    centre = SkyCoord(ra * u.deg, dec * u.deg)
    res = v.query_region(centre, radius=radius_arcsec * u.arcsec, catalog="I/355/gaiadr3")
    if not res:
        return __import__("pandas").DataFrame()
    df = res[0].to_pandas()
    sep = centre.separation(SkyCoord(df.RA_ICRS.values * u.deg, df.DE_ICRS.values * u.deg)).arcsec
    return __import__("pandas").DataFrame({
        "source_id": df.Source.astype("int64"), "ra": df.RA_ICRS, "dec": df.DE_ICRS,
        "phot_g_mean_mag": df.Gmag, "phot_rp_mean_mag": df.RPmag, "ruwe": df.RUWE,
        "non_single_star": df.NSS, "parallax": df.Plx, "sep_arcsec": sep,
    }).sort_values("phot_g_mean_mag").reset_index(drop=True)


def neighbour_check(ra, dec, depth, target_gaia_id=None):
    """Which neighbours could produce a dip of fractional `depth` in the target's aperture?

    Uses Gaia RP magnitudes (closest Gaia band to TESS) when available, else G.
    Needed eclipse depth for neighbour i: depth * (F_target + sum F_others) / F_i,
    approximated with flux of target + that neighbour only (conservative: lower).
    """
    df = gaia_neighbours(ra, dec)
    if df.empty:
        return {"n_gaia": 0, "neighbours": [], "possible_sources": 0}
    mag = df.phot_rp_mean_mag.fillna(df.phot_g_mean_mag)
    if target_gaia_id is not None and int(target_gaia_id) in set(df.source_id.astype(np.int64)):
        is_target = df.source_id.astype(np.int64) == int(target_gaia_id)
    else:
        is_target = df.sep_arcsec == df.sep_arcsec.min()
    m_t = float(mag[is_target].iloc[0])
    target = df[is_target].iloc[0]
    out = []
    for (_, row), m in zip(df.iterrows(), mag):
        if row.source_id == target.source_id or not np.isfinite(m):
            continue
        flux_ratio = 10 ** (-0.4 * (m - m_t))           # neighbour / target
        needed = depth * (1 + flux_ratio) / flux_ratio   # eclipse depth the neighbour would need
        out.append({"source_id": int(row.source_id), "sep_arcsec": float(row.sep_arcsec),
                    "delta_mag": float(m - m_t), "needed_depth": float(needed),
                    "could_be_source": bool(needed < 0.8 and row.sep_arcsec < SEARCH_RADIUS_ARCSEC)})
    return {
        "n_gaia": int(len(df)),
        "target_ruwe": float(target.ruwe) if np.isfinite(target.ruwe) else None,
        "target_non_single_star": int(target.non_single_star) if np.isfinite(target.non_single_star) else None,
        "neighbours": out,
        "possible_sources": int(sum(n["could_be_source"] for n in out)),
    }
