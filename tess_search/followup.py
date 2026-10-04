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
