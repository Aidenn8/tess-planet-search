"""Step 18: do public ground-based light curves already cover a predicted transit?

TOI-218 has 22 TFOP time-series observations on ExoFOP (taken for TOI-218.01 and .02), with public
AstroImageJ measurement tables. This script downloads the tables for a candidate's star, finds every
observation that covers a predicted transit of the candidate (ephemeris from results/candidates/summary.csv),
and measures the dip on the target, on the second target aperture (for TOI-218, its wide-binary companion),
and in the ratio of the two, against two baselines (flat; linear in time).

    .venv/bin/python scripts/18_archival_ground_check.py [TIC]        (default 32090583)

Other public light curves (measurement files, joint-fit subsets) are checked for coverage only.
Writes results/hardening/archival/TIC<tic>_ground_check.json and a figure per covering observation.
"""
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tess_search import DATA, RESULTS

OUT = RESULTS / "hardening" / "archival"
EXOFOP = "https://exofop.ipac.caltech.edu/tess"


def fetch(url, dest):
    if not dest.exists():
        subprocess.run(["curl", "-s", "-f", "-L", "-o", str(dest), url], check=True)
    return dest


def dip(t_h, flux, t14_h, baseline):
    """Out-of-transit minus in-transit mean (fractional), with an error from the out-of-transit scatter.
    baseline: 'flat' or 'linear' (linear in time, fitted to out-of-transit points only)."""
    inn = np.abs(t_h) < 0.4 * t14_h
    out = np.abs(t_h) > t14_h / 2 + 0.1
    if inn.sum() < 5 or out.sum() < 10:
        return np.nan, np.nan
    f = flux / np.median(flux[out])
    if baseline == "flat":
        sd = np.std(f[out])
        return float(np.mean(f[out]) - np.mean(f[inn])), float(sd * np.sqrt(1 / inn.sum() + 1 / out.sum()))
    c = np.polyfit(t_h[out], f[out], 1)
    r = f / np.polyval(c, t_h)
    sd = np.std(r[out])
    X = np.vstack([t_h[out], np.ones(out.sum())]).T
    cov = np.linalg.inv(X.T @ X) * sd ** 2
    var_b = np.mean([x @ cov @ x for x in np.vstack([t_h[inn], np.ones(inn.sum())]).T])
    return float(1 - np.mean(r[inn])), float(np.sqrt(sd ** 2 / inn.sum() + var_b))


if __name__ == "__main__":
    tic = int(sys.argv[1]) if len(sys.argv) > 1 else 32090583
    cs = pd.read_csv(RESULTS / "candidates" / "summary.csv").set_index("tic")
    c = cs.loc[tic]
    P, P_err, T0, T0_err, T14_h = c.period, c.period_err, c.t0_bjd, c.t0_err, c.t14_h
    depth = c.depth_ppm * 1e-6
    work = DATA / "exofop" / str(tic)
    work.mkdir(parents=True, exist_ok=True)
    info = json.loads(fetch(f"{EXOFOP}/target.php?id={tic}&json", work / "target.json").read_text())
    tables = [f for f in info.get("files", []) if f["fext"] == "tbl" and f["ftype"] == "Light_Curve"]
    OUT.mkdir(parents=True, exist_ok=True)
    results = []
    for f in tables:
        path = fetch(f"{EXOFOP}/get_file.php?id={f['fid']}", work / f["fname"])
        d = pd.read_csv(path, sep="\t")
        if "BJD_TDB" not in d:
            continue
        t = d["BJD_TDB"].to_numpy()
        n0, n1 = int(np.floor((t.min() - T0) / P)), int(np.ceil((t.max() - T0) / P))
        for n in range(n0, n1 + 1):
            tc = T0 + n * P
            ing, egr = tc - T14_h / 48, tc + T14_h / 48
            covered = (min(egr, t.max()) - max(ing, t.min())) / (T14_h / 24)
            if covered <= 0.5:
                continue
            th = (t - tc) * 24
            rec = {"file": f["fname"], "tag": f["ftag"], "epoch": n, "tc_bjd": tc,
                   "tc_sigma_min": float(np.hypot(T0_err, n * P_err) * 1440),
                   "transit_fraction_covered": round(float(min(covered, 1.0)), 2),
                   "pre_transit_baseline_h": round(float(max(ing - t.min(), 0) * 24), 2),
                   "post_transit_baseline_h": round(float(max(t.max() - egr, 0) * 24), 2),
                   "tess_depth_ppm": round(depth * 1e6)}
            series = {"target (T1)": d["rel_flux_T1"].to_numpy()}
            if "rel_flux_T2" in d:
                series["second aperture (T2)"] = d["rel_flux_T2"].to_numpy()
                series["T1 / T2 ratio"] = (d["Source-Sky_T1"] / d["Source-Sky_T2"]).to_numpy()
            for name, flux in series.items():
                for base in ("flat", "linear"):
                    dp, er = dip(th, flux, T14_h, base)
                    rec[f"{name} | {base}"] = {"depth_ppm": round(dp * 1e6), "err_ppm": round(er * 1e6)}
            results.append(rec)
            if "T1 / T2 ratio" in series:
                r = series["T1 / T2 ratio"]
                out = np.abs(th) > T14_h / 2 + 0.1
                r = r / np.median(r[out])
                fig, ax = plt.subplots(figsize=(7.4, 3.0))
                ax.axvspan(-T14_h / 2, T14_h / 2, color="#e1e0d9", lw=0, label="predicted transit")
                ax.plot(th, (r - 1) * 1e3, "o", ms=3, color="#52514e", label="TOI-218 / companion, per image")
                nb = 6
                k = len(th) // nb * nb
                ax.plot(th[:k].reshape(-1, nb).mean(1), (r[:k].reshape(-1, nb).mean(1) - 1) * 1e3, "o-",
                        color="#2a78d6", ms=6, label="binned")
                ax.axhline(-depth * 1e3, color="#eb6834", lw=1.2, ls=(0, (3, 2)), label="TESS depth")
                ax.set_xlabel("hours from predicted mid-transit")
                ax.set_ylabel("relative flux (ppt)")
                ax.set_title(f"{f['fname'][:40]}...  epoch {n}, timing ±{rec['tc_sigma_min']:.1f} min", fontsize=9)
                ax.legend(fontsize=7.5, loc="lower right")
                fig.tight_layout()
                fig.savefig(OUT / f"TIC{tic}_{Path(f['fname']).stem}_epoch{n}.png", dpi=150)
                plt.close(fig)
    # other public light curves (e.g. TRAPPIST measurement files, joint-fit subsets): coverage only
    others = [f for f in info.get("files", []) if f["ftype"] == "Light_Curve" and f["fext"] in ("txt", "dat")
              and ("measur" in f["fname"].lower() or "subset" in f["fname"].lower())]
    other_cov = []
    for f in others:
        path = fetch(f"{EXOFOP}/get_file.php?id={f['fid']}", work / f["fname"])
        try:
            d = pd.read_csv(path, sep=None, engine="python")
        except Exception:
            continue
        col = next((c for c in d.columns if "BJD" in str(c).upper()), None)
        if col is None:
            continue
        t = pd.to_numeric(d[col], errors="coerce").dropna().to_numpy()
        t = t + 2400000 if t.max() < 100000 else t
        n = np.round((t - T0) / P)
        hours = (t - (T0 + n * P)) * 24
        inside = np.abs(hours) < T14_h / 2
        other_cov.append({"file": f["fname"], "tag": f["ftag"], "points_in_transit": int(inside.sum()),
                          "fraction_of_points_in_transit": round(float(inside.mean()), 3)})
    (OUT / f"TIC{tic}_ground_check.json").write_text(json.dumps(
        {"tic": tic, "tables_checked": len(tables), "covering_observations": results,
         "other_light_curves_coverage": other_cov}, indent=1, default=float))
    print(f"{len(tables)} tables checked; {len(results)} cover more than half of a predicted transit")
    for o in other_cov:
        print("other:", o)
    for r in results:
        print(json.dumps(r, default=float))
