"""Injection-recovery: plant fake planets in real data and see if the search finds them.

This is how we measure what the search can and cannot detect. Each injected
transit is a physical model (batman, Kreidberg 2015) with limb darkening
typical of M dwarfs in the TESS band, multiplied into the raw PDCSAP and SAP
fluxes before any cleaning or detrending, so it suffers everything a real
planet would.
"""
import numpy as np

G = 6.674e-11
R_SUN, M_SUN, R_EARTH = 6.957e8, 1.989e30, 6.371e6
LIMB_DARKENING = (0.20, 0.40)  # quadratic u1, u2, roughly Claret (2017) TESS for Teff~3500 K


def a_over_r(period_days, r_star, m_star):
    P = period_days * 86400.0
    a = (G * m_star * M_SUN * P ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
    return a / (r_star * R_SUN)


def transit_model(time, period, t0, rp_rearth, r_star, m_star, b=0.3):
    import batman

    params = batman.TransitParams()
    params.t0, params.per = t0, period
    params.rp = rp_rearth * R_EARTH / (r_star * R_SUN)
    params.a = a_over_r(period, r_star, m_star)
    params.inc = float(np.degrees(np.arccos(min(b / params.a, 1.0))))
    params.ecc, params.w = 0.0, 90.0
    params.u, params.limb_dark = list(LIMB_DARKENING), "quadratic"
    m = batman.TransitModel(params, np.ascontiguousarray(time, dtype=np.float64))
    return m.light_curve(params)


def expected_duration(period, r_star, m_star, b=0.0):
    """Full transit duration (days) for a circular orbit, small planet."""
    aR = a_over_r(period, r_star, m_star)
    return period / np.pi * np.arcsin(np.sqrt(max(1 - b ** 2, 0)) / aR)


def inject(raw, period, t0, rp_rearth, r_star, m_star, b=0.3, dilution=1.0):
    """Return a copy of a load_star() dict with a transit multiplied into the fluxes.

    dilution = fraction of aperture flux from the target (CROWDSAP); PDCSAP is
    already corrected for crowding, SAP is not, so SAP gets the diluted signal.
    """
    out = dict(raw)
    model = transit_model(raw["time"], period, t0, rp_rearth, r_star, m_star, b)
    out["pdcsap"] = (raw["pdcsap"] * model).astype(np.float32)
    out["sap"] = (raw["sap"] * (1 - dilution * (1 - model))).astype(np.float32)
    return out


def synthetic_star(seed=0, noise_ppm=600.0, seasons=((1325, 1680), (2035, 2390), (2760, 3000)),
                   sector_days=27.4, gap_days=1.0):
    """A fake load_star() dict: white noise, 2-min cadence, sectors with mid-sector gaps."""
    rng = np.random.default_rng(seed)
    times, sectors = [], []
    sec = 1
    for start, stop in seasons:
        t = start
        while t + sector_days <= stop:
            half = (sector_days - gap_days) / 2
            for a, b in ((t, t + half), (t + half + gap_days, t + sector_days)):
                tt = np.arange(a, b, 2 / 1440)
                times.append(tt)
                sectors.append(np.full(len(tt), sec, np.int16))
            t += sector_days
            sec += 1
    time = np.concatenate(times)
    n = len(time)
    flux = (1 + rng.normal(0, noise_ppm * 1e-6, n)) * 5000.0
    err = np.full(n, noise_ppm * 1e-6 * 5000.0)
    return {
        "time": time, "pdcsap": flux.astype(np.float32), "pdcsap_err": err.astype(np.float32),
        "sap": flux.astype(np.float32), "bkg": np.zeros(n, np.float32),
        "quality": np.zeros(n, np.int32), "centr1": np.zeros(n, np.float32),
        "centr2": np.zeros(n, np.float32), "sector": np.concatenate(sectors),
        "sector_meta": [], "star_meta": {},
    }
