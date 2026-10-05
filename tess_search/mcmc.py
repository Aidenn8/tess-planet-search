"""Transit parameters with full uncertainties: an MCMC fit (emcee) of a physical transit model.

The earlier least-squares fit (followup.fit_transit) gives best values and
covariance errors, ignores the uncertainty in the star's radius, and treats the
limb darkening as known. Here:

  * free parameters: mid-transit offset, Rp/R*, impact parameter b, stellar
    density rho*, and the quadratic limb-darkening coefficients (Gaussian priors
    around typical M-dwarf values in the TESS band, sigma 0.1, kept physical)
  * rho* has a Gaussian prior from the TESS Input Catalog (M/R^3; 3% radius error
    from the TIC, 5% mass error assumed for the Mann et al. 2019 relation), which
    sets a/R* through Kepler's third law; a second fit with a wide prior on rho*
    measures the density from the transit shape alone, as a check that the
    transit is on a star like the target (a blend on a different star, or an
    eccentric orbit, would pull it away)
  * the light curve is re-detrended with the transits masked (the search's filter sags
    slightly into deep transits), folded at the period and binned to 1 minute; bin errors
    are the point scatter / sqrt(n), multiplied by the red-noise factor measured at
    the transit's own timescale
  * planet radius, semi-major axis, equilibrium temperature and insolation are
    computed per sample, with the stellar radius drawn from its TIC uncertainty
"""
import numpy as np

from .lightcurve import robust_std

G_CGS = 6.674e-8
RHO_SUN = 1.41          # g/cm^3
R_SUN_AU = 0.00465047
LD_PRIOR = ((0.20, 0.10), (0.40, 0.10))   # u1, u2 (mean, sigma)
MASS_FRAC_ERR = 0.05


def remask_detrend(lc, period, t0, t14, mask_factor=1.5):
    """Detrend again with the transits masked. The search's robust filter sees the
    in-transit points and, for deep transits, sags slightly into them, which makes the
    transit shallower and more V-shaped. Here the trend is fitted to out-of-transit
    points only and interpolated across each transit."""
    import copy

    from wotan import flatten

    dt = ((lc.time - t0 + 0.5 * period) % period) - 0.5 * period
    out = np.abs(dt) > mask_factor * t14
    window = (lc.star_meta or {}).get("window", 0.5)
    trend = np.full(len(lc.time), np.nan)
    _, tr = flatten(lc.time[out], lc.raw[out], method="biweight", window_length=window, break_tolerance=0.3,
                    edge_cutoff=0.0, return_trend=True)
    good = np.isfinite(tr)
    trend = np.interp(lc.time, lc.time[out][good], tr[good])
    new = copy.copy(lc)
    new.flux = lc.raw / trend
    return new


def fold_bin(lc, period, t0, t14, window=3.0, bin_min=1.0):
    dt = ((lc.time - t0 + 0.5 * period) % period) - 0.5 * period
    m = np.abs(dt) < window * t14
    edges = np.arange(-window * t14, window * t14 + 1e-12, bin_min / 1440)
    idx = np.clip(np.digitize(dt[m], edges) - 1, 0, len(edges) - 2)
    n = np.bincount(idx, minlength=len(edges) - 1)
    s = np.bincount(idx, weights=lc.flux[m], minlength=len(edges) - 1)
    ok = n > 0
    sigma_pt = robust_std(lc.flux[~m] - 1)
    # red-noise factor at the transit timescale: scatter of T14-long means vs white expectation
    tt = lc.time[~m]
    ff = lc.flux[~m]
    k = np.floor((tt - tt[0]) / t14).astype(np.int64)
    _, inv, cnt = np.unique(k, return_inverse=True, return_counts=True)
    means = np.bincount(inv, weights=ff) / cnt
    full = cnt >= 0.6 * t14 / (2 / 1440)
    beta = 1.0
    if full.sum() > 20:
        beta = max(1.0, robust_std(means[full] - 1) / (sigma_pt / np.sqrt(np.median(cnt[full]))))
    centres = 0.5 * (edges[1:] + edges[:-1])
    return centres[ok], s[ok] / n[ok], sigma_pt * beta / np.sqrt(n[ok]), beta


def a_over_r(rho, period_days):
    return (G_CGS * rho * (period_days * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)


class TransitPosterior:
    def __init__(self, t, f, e, period, t14, rho_prior):
        import batman

        self.t, self.f, self.ivar = t, f, 1 / e ** 2
        self.period, self.t14 = period, t14
        self.rho_prior = rho_prior            # (mean, sigma) g/cm^3, or None for a wide prior
        self.params = batman.TransitParams()
        p = self.params
        p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w = 0.0, period, 0.03, 20.0, 89.0, 0.0, 90.0
        p.u, p.limb_dark = [0.2, 0.4], "quadratic"
        # 1-minute bins of 2-minute cadences: smooth the model over ~2.5 minutes
        self.model = batman.TransitModel(p, t, supersample_factor=5, exp_time=2.5 / 1440)

    def curve(self, theta):
        dt0, k, b, log_rho, u1, u2 = theta
        p = self.params
        p.t0, p.rp = dt0, k
        p.a = a_over_r(10 ** log_rho, self.period)
        p.inc = float(np.degrees(np.arccos(np.clip(b / p.a, 0, 1))))
        p.u = [u1, u2]
        return self.model.light_curve(p)

    def log_prior(self, theta):
        dt0, k, b, log_rho, u1, u2 = theta
        if not (abs(dt0) < 0.5 * self.t14 and 0.001 < k < 0.3 and 0 <= b < 1 + k and -1.5 < log_rho < 2.5):
            return -np.inf
        # physical limb darkening (Kipping 2013): u1 + u2 < 1, u1 > 0, u1 + 2 u2 > 0
        if not (u1 > 0 and u1 + u2 < 1 and u1 + 2 * u2 > 0):
            return -np.inf
        lp = -0.5 * (((u1 - LD_PRIOR[0][0]) / LD_PRIOR[0][1]) ** 2 + ((u2 - LD_PRIOR[1][0]) / LD_PRIOR[1][1]) ** 2)
        if self.rho_prior is not None:
            rho = 10 ** log_rho
            lp += -0.5 * ((rho - self.rho_prior[0]) / self.rho_prior[1]) ** 2 + np.log(rho)  # prior on rho, sampled in log
        return lp

    def __call__(self, theta):
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        r = self.f - self.curve(theta)
        return lp - 0.5 * np.sum(r * r * self.ivar)


def fit(lc, period, t0, t14, depth, star, nwalkers=32, nsteps=6000, burn=2000, density_prior=True, seed=1,
        remask=True):
    """star: dict with rad, e_rad, mass, teff. Returns summary dict and the flat chain."""
    import emcee

    if remask:
        lc = remask_detrend(lc, period, t0, t14)
    t, f, e, beta = fold_bin(lc, period, t0, t14)
    rho_tic = star["mass"] / star["rad"] ** 3 * RHO_SUN
    rho_err = rho_tic * np.hypot(MASS_FRAC_ERR, 3 * star["e_rad"] / star["rad"])
    post = TransitPosterior(t, f, e, period, t14, (rho_tic, rho_err) if density_prior else None)
    rng = np.random.default_rng(seed)
    p0 = np.column_stack([
        rng.normal(0, 0.02 * t14, nwalkers),
        np.abs(rng.normal(np.sqrt(max(depth, 1e-5)), 0.002, nwalkers)),
        rng.uniform(0.1, 0.6, nwalkers),
        rng.normal(np.log10(rho_tic), 0.03, nwalkers),
        rng.normal(LD_PRIOR[0][0], 0.02, nwalkers),
        rng.normal(LD_PRIOR[1][0], 0.02, nwalkers),
    ])
    sampler = emcee.EnsembleSampler(nwalkers, 6, post)
    sampler.run_mcmc(p0, nsteps, progress=False)
    try:
        tau = float(np.max(sampler.get_autocorr_time(discard=burn, quiet=True)))
    except Exception:
        tau = np.nan
    chain = sampler.get_chain(discard=burn, thin=5, flat=True)
    lnp = sampler.get_log_prob(discard=burn, thin=5, flat=True)
    chain = chain[np.isfinite(lnp)]

    dt0, k, b, log_rho, u1, u2 = chain.T
    rho = 10 ** log_rho
    aR = a_over_r(rho, period)
    inc = np.degrees(np.arccos(np.clip(b / aR, 0, 1)))
    arg = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / aR / np.sin(np.radians(inc))
    t14_h = period / np.pi * np.arcsin(np.clip(arg, 0, 1)) * 24
    rs = rng.normal(star["rad"], star["e_rad"], len(chain))
    rp = k * rs * 109.1
    a_au = aR * rs * R_SUN_AU
    teff = rng.normal(star["teff"], star.get("e_teff", 157.0), len(chain))
    teq = teff * np.sqrt(1 / (2 * aR))                          # zero albedo, full redistribution
    insol = rs ** 2 * (teff / 5772.0) ** 4 / a_au ** 2          # Earth units
    depth_ppm = np.array([1 - post.curve(th).min() for th in chain[rng.choice(len(chain), 400, replace=False)]]) * 1e6

    def q(x):
        lo, med, hi = np.percentile(x, [15.865, 50, 84.135])
        return {"median": float(med), "minus": float(med - lo), "plus": float(hi - med)}

    summary = {
        "density_prior": density_prior, "rho_tic_cgs": float(rho_tic), "rho_tic_err_cgs": float(rho_err),
        "n_bins": int(len(t)), "red_noise_beta": float(beta), "n_samples": int(len(chain)),
        "autocorr_steps": tau, "converged": bool(np.isfinite(tau) and (nsteps - burn) > 30 * tau),
        "acceptance": float(np.mean(sampler.acceptance_fraction)),
        "t0_btjd": q(t0 + dt0), "rp_rs": q(k), "b": q(b), "rho_cgs": q(rho), "a_rs": q(aR), "inc_deg": q(inc),
        "t14_h": q(t14_h), "depth_ppm": q(depth_ppm), "rp_rearth": q(rp), "a_au": q(a_au), "teq_k": q(teq),
        "insolation_earth": q(insol), "u1": q(u1), "u2": q(u2),
        "rho_ratio_to_tic": q(rho / rho_tic),
    }
    return summary, chain, (t, f, e), post
