"""Step 14: MCMC transit fits for the candidates, checked first on two known planets.

    .venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/14_mcmc.py

For each signal two fits (tess_search/mcmc.py): with the TIC stellar-density prior
(the adopted parameters) and with a wide density prior (the density the transit
shape alone implies, compared with the TIC value). Validation: TOI-700 d and
L 98-59 c, whose published radii come from independent analyses.
Writes results/hardening/mcmc/<label>.json and <label>.png.
"""
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for var in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(var, "1")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from tess_search import DATA, RESULTS

OUT = RESULTS / "hardening" / "mcmc"
CANDIDATE_FILES = ["TIC32090583_5", "TIC229689348_1", "TIC198412174_1", "TIC149390648_1", "TIC294053492_1"]
# published radii for the validation planets (R_earth): TOI-700 d, Gilbert et al. 2023; L 98-59 c, Demangeon et al. 2021
VALIDATION = [("TOI-700d", 150428135, 3, 1.073, 0.059), ("L98-59c", 307210830, 1, 1.385, 0.085)]


def jobs():
    out = []
    for f in CANDIDATE_FILES:
        d = json.loads((RESULTS / "followup" / f"{f}.json").read_text())
        fit = d["fit"]
        out.append({"label": f"TIC{d['tic']}", "tic": d["tic"], "period": fit["period"], "t0": fit["t0_btjd"],
                    "t14": fit["t14_h"] / 24, "depth": fit["depth_ppm"] * 1e-6})
    sig = pd.read_csv(RESULTS / "all_signals.csv")
    for label, tic, k, rp, rp_err in VALIDATION:
        r = sig[(sig.tic == tic) & (sig.signal == k)].iloc[0]
        out.append({"label": label, "tic": tic, "period": float(r.period), "t0": float(r.t0_btjd),
                    "t14": float(r.duration_h) / 24, "depth": float(r.depth_ppm) * 1e-6,
                    "published_rp": rp, "published_rp_err": rp_err})
    return out


def plot(job, res, data, post, chain, path):
    import matplotlib
    matplotlib.use("Agg")
    import corner
    import matplotlib.pyplot as plt

    t, f, e = data
    fig = plt.figure(figsize=(13, 5.5))
    ax = fig.add_axes([0.06, 0.12, 0.42, 0.78])
    nb = 6
    tb = t[: len(t) // nb * nb].reshape(-1, nb).mean(1)
    fb = f[: len(f) // nb * nb].reshape(-1, nb).mean(1)
    eb = np.sqrt((e[: len(e) // nb * nb].reshape(-1, nb) ** 2).sum(1)) / nb
    ax.errorbar(tb * 24, (fb - 1) * 1e6, eb * 1e6, fmt="o", ms=3, color="k", lw=0.8, label="data, 6-min bins")
    rng = np.random.default_rng(0)
    for th in chain[rng.choice(len(chain), 60, replace=False)]:
        ax.plot(t * 24, (post.curve(th) - 1) * 1e6, color="tab:red", alpha=0.06, lw=1)
    ax.set_xlabel("hours from mid-transit")
    ax.set_ylabel("relative flux (ppm)")
    s = res["prior"]
    ax.set_title(f"{job['label']}: Rp = {s['rp_rearth']['median']:.2f} +{s['rp_rearth']['plus']:.2f} "
                 f"-{s['rp_rearth']['minus']:.2f} R_earth, b = {s['b']['median']:.2f}, "
                 f"T14 = {s['t14_h']['median']:.2f} h", fontsize=10, loc="left")
    ax.legend(fontsize=8)
    sub = fig.add_axes([0.55, 0.08, 0.43, 0.86])
    sub.axis("off")
    k, b, log_rho = chain[:, 1], chain[:, 2], chain[:, 3]
    cf = corner.corner(np.column_stack([k, b, 10 ** log_rho]), labels=["Rp/R*", "b", "rho* (g/cc)"],
                       quantiles=[0.16, 0.5, 0.84], fig=plt.figure(figsize=(5, 5)))
    cf.canvas.draw()
    img = np.asarray(cf.canvas.buffer_rgba())
    sub.imshow(img)
    plt.close(cf)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def run(job, scale=1):
    from tess_search import lightcurve, mcmc

    tic_row = pd.read_parquet(DATA / "tic_15plus.parquet").set_index("tic").loc[job["tic"]]
    star = {"rad": float(tic_row.rad), "e_rad": float(tic_row.e_rad), "mass": float(tic_row.mass),
            "teff": float(tic_row.Teff), "e_teff": float(tic_row.e_Teff)}
    lc = lightcurve.prepare(job["tic"])
    if "published_rp" in job:
        # validation planets start from the search's box-fit period, good to a few minutes over the
        # baseline; refine period and epoch with the least-squares transit fit the candidates had,
        # otherwise the folded ingress is blurred and the fit drifts towards grazing geometries
        from tess_search import followup
        f = followup.fit_transit(lc, job["period"], job["t0"], job["t14"], job["depth"], star["rad"], star["mass"])
        job = dict(job, period=f["period"], t0=f["t0_btjd"], t14=f["t14_h"] / 24)
    res = {"label": job["label"], "tic": job["tic"], "period": job["period"], "star": star}
    s1, chain, data, post = mcmc.fit(lc, job["period"], job["t0"], job["t14"], job["depth"], star, density_prior=True,
                                     nsteps=6000 * scale, burn=2000 * scale)
    # without the density prior b and rho* are strongly correlated: the chain mixes ~3-4x slower
    s2, _, _, _ = mcmc.fit(lc, job["period"], job["t0"], job["t14"], job["depth"], star, density_prior=False, seed=2,
                           nsteps=30000 * scale, burn=8000 * scale)
    res.update(prior=s1, free_density=s2)
    if "published_rp" in job:
        res["published_rp"] = job["published_rp"]
        res["published_rp_err"] = job["published_rp_err"]
    (OUT / f"{job['label']}.json").write_text(json.dumps(res, indent=1, default=float))
    plot(job, res, data, post, chain, OUT / f"{job['label']}.png")
    r = s1["rp_rearth"]
    rho = s2["rho_ratio_to_tic"]
    print(f"{job['label']}: Rp {r['median']:.2f} +{r['plus']:.2f} -{r['minus']:.2f}  b {s1['b']['median']:.2f}  "
          f"converged {s1['converged']}  free-density rho/rho_TIC {rho['median']:.2f} "
          f"+{rho['plus']:.2f} -{rho['minus']:.2f}" + (f"  published {job['published_rp']}" if "published_rp" in job else ""),
          flush=True)
    return res


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", help="labels to (re)fit")
    ap.add_argument("--scale", type=int, default=1, help="multiply chain lengths (high-SNR targets mix slower)")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    todo = [j for j in jobs() if not args.only or j["label"] in args.only]
    with ProcessPoolExecutor(min(3, len(todo))) as pool:
        list(pool.map(run, todo, [args.scale] * len(todo)))
