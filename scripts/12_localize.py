"""Step 12: pixel-level localization of the candidates, and the tests that show it works.

For every signal: where on the sky does the light go missing? (tess_search/localize.py)

  candidates   the five candidates
  planets      confirmed planets: the source must come out on the target
  nebs         TOIs that TFOP retired as nearby eclipsing binaries: the source must come
               out off the target (TOI-419.01 on TIC 279251647 specifically)
  injections   for each candidate star, a synthetic eclipse planted in the real pixels on
               its most dangerous neighbours (and on the target as a control), sized to
               reproduce the candidate's observed depth: the method must find each one

    .venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/12_localize.py
    .venv/bin/python scripts/12_localize.py --only 229689348      # one star

Results: results/hardening/localize/<label>.json and .png, plus localization_summary.csv.
"""
import argparse
import json
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for var in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(var, "1")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from tess_search import DATA, RESULTS

OUT = RESULTS / "hardening" / "localize"
CANDIDATE_FILES = ["TIC32090583_5", "TIC229689348_1", "TIC198412174_1", "TIC149390648_1", "TIC294053492_1"]
# (label, tic, signal number in results/all_signals.csv, kind, expected source Gaia DR3 id or None)
VALIDATION = [
    ("TOI-6000b", 259233660, 1, "planet", None),
    ("TOI-5728b", 219875976, 1, "planet", None),
    ("TOI-700c", 150428135, 1, "planet", None),
    ("TOI-700b", 150428135, 2, "planet", None),
    ("TOI-700d", 150428135, 3, "planet", None),
    ("TOI-2094b", 356016119, 1, "planet", None),
    ("TOI-1756b", 364074068, 1, "planet", None),
    ("TOI-2084b", 441738827, 1, "planet", None),
    ("L98-59c", 307210830, 1, "planet", None),
    ("L98-59d", 307210830, 2, "planet", None),
    ("L98-59b", 307210830, 3, "planet", None),
    ("TOI-419.01", 279251651, 1, "neb", 5480636006390943872),   # TIC 279251647 per TFOP
    ("TOI-2084.02", 441738827, 2, "neb", None),
    ("TOI-2283.01", 198211976, 2, "neb", None),
    # context: which star of the TOI-218 wide binary hosts its two known candidates
    ("TOI-218.01", 32090583, 1, "context", None),
    ("TOI-218.02", 32090583, 2, "context", None),
]
# this search's detection of TOI-2283.01 is weak (24 ppm) with a box duration of 4.8 h, half its
# 9.6-h orbit; use the TOI catalogue duration instead (1.578 h)
DURATION_OVERRIDE_H = {"TOI-2283.01": 1.578}
MAX_INJECT_NEIGHBOURS = 3
# Systematic position error (1-D sigma, arcsec) added in quadrature to the statistical error:
# PRF-model and pixel-calibration errors. Measured at the end of each run from every case whose
# true source is known (confirmed planets and injections); the run is repeated if it changes.
SYS_ARCSEC = float(os.environ.get("LOC_SYS_ARCSEC", 1.5))
# The model assumes one source changes brightness in step with the signal. Another variable
# star in the stamp whose own period happens to be commensurate breaks that and shows up as a
# poor fit; such localizations are flagged unreliable (seen in one injection series, TOI-6000).
CHI2_RELIABLE = 2.0


def weak_list():
    """The new weak candidates (passed the hard tests but flagged), from their follow-up fits."""
    c = pd.read_csv(RESULTS / "candidates.csv")
    w = c[(c.verdict == "weak candidate") & c.novelty.isin(["new", "SPOC TCE only (never promoted)"])]
    out = []
    for _, r in w.iterrows():
        d = json.loads((RESULTS / "followup" / f"TIC{r.tic}_{r.signal}.json").read_text())
        fit = d.get("fit", {})
        ok = "error" not in fit and np.isfinite(fit.get("t14_h", np.nan))
        out.append({"label": f"TIC{r.tic}_{r.signal}", "tic": int(r.tic), "kind": "weak", "expected": None,
                    "period": fit["period"] if ok else float(r.period),
                    "t0": fit["t0_btjd"] if ok else float(r.t0_btjd),
                    "t14": (fit["t14_h"] if ok else float(r.duration_h)) / 24,
                    "depth": (fit["depth_ppm"] if ok else float(r.depth_ppm)) * 1e-6})
    return out


def signal_list(only=None, weak=False):
    if weak:
        out = weak_list()
        return [s for s in out if s["tic"] in only] if only else out
    sig = pd.read_csv(RESULTS / "all_signals.csv")
    out = []
    for f in CANDIDATE_FILES:
        d = json.loads((RESULTS / "followup" / f"{f}.json").read_text())
        fit = d["fit"]
        out.append({"label": f"TIC{d['tic']}", "tic": d["tic"], "kind": "candidate", "expected": None,
                    "period": fit["period"], "t0": fit["t0_btjd"], "t14": fit["t14_h"] / 24,
                    "depth": fit["depth_ppm"] * 1e-6})
    for label, tic, k, kind, exp in VALIDATION:
        r = sig[(sig.tic == tic) & (sig.signal == k)].iloc[0]
        t14_h = DURATION_OVERRIDE_H.get(label, float(r.duration_h))
        out.append({"label": label, "tic": tic, "kind": kind, "expected": exp, "period": float(r.period),
                    "t0": float(r.t0_btjd), "t14": t14_h / 24, "depth": float(r.depth_ppm) * 1e-6})
    if only:
        out = [s for s in out if s["tic"] in only]
    return out


def describe(loc, stars, target_idx, sig):
    from tess_search.localize import with_systematics

    best = loc.stars[0]
    t = next(r for r in loc.stars if r["is_target"])
    rec = {k: v for k, v in loc.__dict__.items() if k not in ("stars", "grid")}
    tot = with_systematics({"offset_err_arcsec": loc.offset_err_arcsec, "stars": loc.stars}, SYS_ARCSEC)
    loc.stars = tot["stars"]
    rec.update(label=sig["label"], kind=sig["kind"], period=sig["period"], depth_ppm=sig["depth"] * 1e6,
               best_star=best["source_id"], best_star_sep=best["sep_arcsec"], best_star_is_target=best["is_target"],
               target_dchi2=t["dchi2"],
               amplitude_ratio=loc.amplitude / loc.expected_amplitude if loc.expected_amplitude > 0 else np.nan)
    rec.update(sys_arcsec=SYS_ARCSEC, offset_err_total_arcsec=tot["offset_err_total_arcsec"],
               target_sigma_total=tot["target_sigma_total"], reliable=bool(loc.chi2_red < CHI2_RELIABLE))
    # neighbours the data cannot exclude at 3 sigma (systematic floor included)
    alive = [r for r in loc.stars if not r["is_target"] and r["excluded_sigma_total"] < 3]
    rec["neighbours_not_excluded_3sigma"] = [(r["source_id"], round(r["sep_arcsec"], 1), round(r["tmag"], 2))
                                             for r in alive]
    rec["miss_arcsec"] = float(np.hypot(loc.offset_east - t["east"], loc.offset_north - t["north"]))
    if sig.get("expected"):
        e = next((r for r in loc.stars if r["source_id"] == sig["expected"]), None)
        rec["expected_source_sigma"] = e["excluded_sigma_total"] if e else None
        rec["expected_source_is_best"] = bool(e and e["source_id"] == best["source_id"])
    return rec


def plot(loc, items, stars, target_idx, sig, path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from tess_search.localize import pix, positions_at

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
    from tess_search.localize import R68

    offs = np.array(loc.grid["offsets_arcsec"])
    s_stat = loc.offset_err_arcsec / R68
    d = np.array(loc.grid["dchi2"]) * s_stat ** 2 / (s_stat ** 2 + SYS_ARCSEC ** 2)  # systematic floor included
    a = ax[0]
    im = a.imshow(np.sqrt(np.clip(d, 0, None)), origin="lower", cmap="viridis_r",
                  extent=[offs[0], offs[-1], offs[0], offs[-1]], vmax=12)
    a.contour(offs, offs, d, levels=[2.30, 11.8, 28.7], colors=["w", "w", "w"], linewidths=[1.5, 1, 0.6])
    cosd = np.cos(np.radians(stars.dec.values[target_idx]))
    year = np.median([c.year for c, *_ in items])
    ra, dec = positions_at(stars, year)
    ra0, dec0 = positions_at(stars.iloc[[target_idx]], year)
    e = (ra - ra0[0]) * cosd * 3600
    n = (dec - dec0[0]) * 3600
    m = (np.abs(e) < offs[-1]) & (np.abs(n) < offs[-1])
    a.scatter(e[m], n[m], s=np.clip(200 * 10 ** (-0.4 * (stars.tmag.values[m] - stars.tmag.values[target_idx])), 6, 300),
              facecolors="none", edgecolors="r", lw=1)
    a.plot(0, 0, "r+", ms=14, mew=2)
    a.plot(loc.offset_east, loc.offset_north, "kx", ms=10, mew=2)
    a.set_xlim(offs[-1], offs[0])  # east to the left, as on the sky
    a.set_xlabel('arcsec east of target')
    a.set_ylabel('arcsec north of target')
    a.set_title(f"A. where the light goes missing (incl. {SYS_ARCSEC:.1f}\" systematic)\n"
                "colour: exclusion in sigma; contours 1/3/5 sigma; red circles: Gaia stars", fontsize=9, loc="left")
    fig.colorbar(im, ax=a, fraction=0.046, label="sigma")

    # the sector whose difference image carries the most signal
    snrs = [np.nansum(diff / err ** 2 * np.clip(diff, 0, None)) for _, diff, err, _ in items]
    cal, diff, err, nt = items[int(np.argmax(snrs))]
    x, y = pix(cal, *positions_at(stars, cal.year))
    for k, (img, ttl, cmap) in enumerate(((np.median(cal.stamp.flux, 0), "B. out-of-transit image", "Greys_r"),
                                          (diff / err, "C. difference image / error (out - in)", "RdBu_r"))):
        a = ax[k + 1]
        v = np.nanpercentile(np.abs(img), 99)
        a.imshow(img, origin="lower", cmap=cmap, vmin=-v if k else 0, vmax=v)
        ny, nx = img.shape
        mm = (x > -0.5) & (x < nx - 0.5) & (y > -0.5) & (y < ny - 0.5)
        a.scatter(x[mm], y[mm], s=np.clip(150 * 10 ** (-0.4 * (stars.tmag.values[mm] - stars.tmag.values[target_idx])), 5, 200),
                  facecolors="none", edgecolors="orange", lw=1)
        a.plot(x[target_idx], y[target_idx], "r+", ms=12, mew=2)
        ap = cal.stamp.aperture
        for (r, c) in zip(*np.nonzero(ap)):
            a.add_patch(plt.Rectangle((c - 0.5, r - 0.5), 1, 1, fill=False, ec="c", lw=0.5))
        a.set_title(f"{ttl}\nsector {cal.stamp.sector}, {nt} transits", fontsize=9, loc="left")
    fig.suptitle(title, fontsize=10, x=0.01, ha="left")
    t_tot = next(r for r in loc.stars if r["is_target"])["excluded_sigma_total"]
    err_tot = R68 * np.hypot(s_stat, SYS_ARCSEC)
    fig.text(0.01, 0.01, f"best position {loc.offset_arcsec:.1f} +- {err_tot:.1f}\" from target "
             f"(E {loc.offset_east:+.0f}\", N {loc.offset_north:+.0f}\"); target excluded at {t_tot:.1f} sigma; "
             f"light lost {loc.amplitude:.2f} +- {loc.amplitude_err:.2f} e-/s vs {loc.expected_amplitude:.2f} expected "
             f"if on target; {loc.n_sectors} sectors, {loc.n_transits} transits; reduced chi2 {loc.chi2_red:.2f}",
             fontsize=8)
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(path, dpi=110)
    plt.close(fig)


def run_star(tic, sigs, inject=True):
    from tess_search import localize, pixels

    t_start = time.time()
    tic_row = pd.read_parquet(DATA / "tic_15plus.parquet").set_index("tic").loc[tic]
    stamps = pixels.load_stamps(tic)
    stars = localize.gaia_stars(float(tic_row.ra), float(tic_row.dec))
    match = np.flatnonzero(stars.source_id.values == int(tic_row.GAIA))
    target_idx = int(match[0]) if len(match) else int(np.argmin(stars.sep_arcsec.values))
    cals = localize.calibrate(stamps, stars, target_idx)
    out = []
    for sig in sigs:
        res = localize.localize(tic, cals, stars, target_idx, sig["period"], sig["t0"], sig["t14"], sig["depth"])
        if res is None:
            out.append({"label": sig["label"], "error": "no usable difference images"})
            continue
        loc, items = res
        rec = describe(loc, stars, target_idx, sig)
        rec["calibration"] = [c.fit for c in cals]
        if sig["kind"] == "candidate":
            rec["split_halves"] = split_halves(tic, cals, stars, target_idx, sig)
        (OUT / f"{sig['label']}.json").write_text(json.dumps(
            dict(rec, stars=loc.stars[:40], grid=loc.grid), default=float, indent=1))
        plot(loc, items, stars, target_idx, sig, OUT / f"{sig['label']}.png",
             f"{sig['label']} (TIC {tic}), P = {sig['period']:.5f} d, depth {sig['depth'] * 1e6:.0f} ppm [{sig['kind']}]")
        out.append(rec)

        first_planet = sig["kind"] == "planet" and sig is next(x for x in sigs if x["kind"] == "planet")
        if inject and (sig["kind"] == "candidate" or first_planet):
            out.extend(injections(tic, cals, stars, target_idx, sig))
    print(f"TIC {tic}: {len(sigs)} signals in {time.time() - t_start:.0f} s", flush=True)
    return out


def split_halves(tic, cals, stars, target_idx, sig):
    """Localize again from the odd-numbered and the even-numbered sectors separately: a real
    source position must not depend on which half of the data is used."""
    from tess_search import localize

    out = []
    for name, sel in (("odd sectors", cals[0::2]), ("even sectors", cals[1::2])):
        res = localize.localize(tic, sel, stars, target_idx, sig["period"], sig["t0"], sig["t14"], sig["depth"])
        if res is None:
            continue
        loc, _ = res
        tot = localize.with_systematics({"offset_err_arcsec": loc.offset_err_arcsec, "stars": loc.stars}, SYS_ARCSEC)
        out.append({"half": name, "offset_east": loc.offset_east, "offset_north": loc.offset_north,
                    "offset_err_total_arcsec": tot["offset_err_total_arcsec"],
                    "target_sigma_total": tot["target_sigma_total"], "chi2_red": loc.chi2_red,
                    "n_sectors": loc.n_sectors})
    return out


def injections(tic, cals, stars, target_idx, sig):
    """Plant eclipses that would reproduce this candidate's depth on its riskiest neighbours."""
    from tess_search import localize

    f_t = np.median([c.fit["target_flux_total"] for c in cals])
    frac_t = localize.aperture_fraction(cals, stars, target_idx)
    scale = np.median([c.fit["scale"] for c in cals])
    cand = []
    for k in range(len(stars)):
        if k == target_idx or stars.sep_arcsec.values[k] > 90:
            continue
        frac_k = localize.aperture_fraction(cals, stars, k)
        if frac_k < 1e-3:
            continue
        amp = sig["depth"] * f_t * frac_t / frac_k          # e-/s the neighbour must lose
        needed = amp / (localize.star_flux(stars.tmag.values[k]) * scale)
        if needed < 0.5:                                      # an eclipse this deep is possible
            cand.append((needed, k, amp))
    cand.sort()
    p_inj = sig["period"] * 1.3713
    t0_inj = sig["t0"] + 0.37 * p_inj
    trials = [(None, target_idx, sig["depth"] * f_t)] + cand[:MAX_INJECT_NEIGHBOURS]
    from tess_search.localize import R68, sigma_from_dchi2

    out = []
    for needed, k, amp in trials:
        icals = localize.inject_eclipse(cals, stars, k, p_inj, t0_inj, sig["t14"], amp)
        res = localize.localize(tic, icals, stars, target_idx, p_inj, t0_inj, sig["t14"], sig["depth"])
        if res is None:
            continue
        loc, _ = res
        src = int(stars.source_id.values[k])
        row = next((r for r in loc.stars if r["source_id"] == src), None)
        t_row = next(r for r in loc.stars if r["is_target"])
        out.append({"label": f"{sig['label']}_inject_{'target' if k == target_idx else src}", "kind": "injection",
                    "tic": tic, "injected_on": "target" if k == target_idx else "neighbour",
                    "injected_source": src, "injected_sep_arcsec": float(stars.sep_arcsec.values[k]),
                    "injected_tmag": float(stars.tmag.values[k]), "needed_eclipse_depth": needed,
                    "recovered_best_star": loc.stars[0]["source_id"],
                    "recovered_correct": loc.stars[0]["source_id"] == src,
                    "injected_source_sigma": row["excluded_sigma"] if row else None,
                    "target_sigma": loc.target_sigma, "offset_arcsec": loc.offset_arcsec,
                    "offset_err_arcsec": loc.offset_err_arcsec, "chi2_red": loc.chi2_red,
                    "target_dchi2": t_row["dchi2"], "injected_dchi2": row["dchi2"] if row else None,
                    "miss_arcsec": float(np.hypot(loc.offset_east - row["east"], loc.offset_north - row["north"]))
                    if row else None})
        s_stat = loc.offset_err_arcsec / R68
        shrink = s_stat ** 2 / (s_stat ** 2 + SYS_ARCSEC ** 2)
        out[-1].update(reliable=bool(loc.chi2_red < CHI2_RELIABLE),
                       target_sigma_total=sigma_from_dchi2(t_row["dchi2"] * shrink),
                       injected_source_sigma_total=sigma_from_dchi2(row["dchi2"] * shrink) if row else None)
    return out


def calibrate_systematics(df):
    """Systematic floor s such that, over every case with a known true source, the misses
    are consistent with the total error: mean(miss^2 / (s_stat^2 + s^2)) = 2 (2-D chi-square)."""
    from scipy.optimize import brentq

    from tess_search.localize import R68

    known = df[df.kind.isin(["planet", "injection"]) & df.miss_arcsec.notna() & (df.chi2_red < CHI2_RELIABLE)]
    if len(known) < 5:
        return
    miss = known.miss_arcsec.to_numpy(float)
    s_stat = known.offset_err_arcsec.to_numpy(float) / R68
    f = lambda s: np.mean(miss ** 2 / (s_stat ** 2 + s ** 2)) - 2.0
    sys_fit = 0.0 if f(0.0) <= 0 else brentq(f, 0.0, 60.0)
    pulls = miss / np.sqrt(s_stat ** 2 + sys_fit ** 2)
    flagged = df[df.kind.isin(["planet", "injection"]) & (df.chi2_red >= CHI2_RELIABLE)]
    from scipy import stats as st
    lim = np.sqrt(st.chi2.ppf(0.9545, 2))
    cover = [float(np.mean(miss / np.sqrt(s_stat ** 2 + s ** 2) < lim)) for s in np.arange(0, 5.01, 0.1)]
    s95 = float(np.arange(0, 5.01, 0.1)[next(i for i, c in enumerate(cover) if c >= 0.95)]) if max(cover) >= 0.95 else None
    res = {"n_cases": int(len(known)), "sys_arcsec_fit": float(sys_fit), "sys_arcsec_used": SYS_ARCSEC,
           "sys_arcsec_95pct_coverage": s95, "flagged_unreliable": flagged.label.tolist(),
           "pull_rms_stat_only": float(np.sqrt(np.mean((miss / s_stat) ** 2) / 2)),
           "frac_within_2_sigma_total": float(np.mean(miss / np.sqrt(s_stat ** 2 + SYS_ARCSEC ** 2) < lim)),
           "cases": known[["label", "kind", "miss_arcsec", "offset_err_arcsec"]].to_dict("records")}
    (OUT.parent / "localization_systematics.json").write_text(json.dumps(res, indent=1, default=float))
    print(f"systematic floor: fitted {sys_fit:.2f}\", 95% coverage needs {s95}\" from {len(known)} known cases "
          f"(used {SYS_ARCSEC:.2f}\"); flagged unreliable: {flagged.label.tolist()}", flush=True)


def stats_chi2_ppf(q):
    from scipy import stats
    return stats.chi2.ppf(q, 2)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", type=int, nargs="*")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--no-inject", action="store_true")
    ap.add_argument("--weak", action="store_true", help="the weak candidates (separate summary file)")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    sigs = signal_list(set(args.only) if args.only else None, weak=args.weak)
    by_star = {}
    for s in sigs:
        by_star.setdefault(s["tic"], []).append(s)
    rows = []
    with ProcessPoolExecutor(args.workers) as pool:
        futs = {pool.submit(run_star, tic, ss, not args.no_inject): tic for tic, ss in by_star.items()}
        for f in as_completed(futs):
            try:
                rows.extend(f.result())
            except Exception:
                print(f"TIC {futs[f]} failed:\n{traceback.format_exc()}", flush=True)
    summ = OUT.parent / ("localization_weak.csv" if args.weak else "localization_summary.csv")
    df = pd.DataFrame(rows)
    if summ.exists() and args.only:  # update those rows in the full summary
        old = pd.read_csv(summ)
        df = pd.concat([old[~old.label.isin(df.label)], df], ignore_index=True)
    df.to_csv(summ, index=False)
    if not args.only and not args.weak:
        calibrate_systematics(df)
    cols = [c for c in ["label", "kind", "offset_arcsec", "offset_err_arcsec", "target_sigma", "best_star_is_target",
                        "amplitude_ratio", "chi2_red", "expected_source_is_best", "recovered_correct",
                        "injected_sep_arcsec", "injected_source_sigma"] if c in df.columns]
    print(df[cols].to_string())
