"""Step 8: deeper checks on every signal that passed vetting and is not already known.

    .venv/bin/python scripts/08_followup.py

For each new candidate / weak candidate (and SPOC-TCE-only ones):
  * Gaia DR3 neighbours: could a nearby star produce the dip? (needed eclipse depth)
  * the target's Gaia RUWE and non-single-star flag (hidden companions)
  * the signal searched for again, independently, in the first and second half of
    the data (a real planet must appear in both, at the same period)
  * a fresh diagnostic sheet
Writes results/followup/<tic>_<n>.json and .png, and results/followup/summary.csv.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from tess_search import DATA, RESULTS, crossmatch, followup, lightcurve, report, search, vetting

OUT = RESULTS / "followup"
OUT.mkdir(parents=True, exist_ok=True)


def half_search(lc, period, r_star, m_star):
    """Search each half of the time series near the candidate period; return SNRs."""
    mid = np.median(lc.time)
    out = []
    for name, m in (("first half", lc.time < mid), ("second half", lc.time >= mid)):
        sub = lc.subset(m)
        if len(sub.time) < 5000:
            out.append({"half": name, "snr": None})
            continue
        P, t0, dur, depth, err, _ = search.refine(sub, period, rel_width=0.002)
        out.append({"half": name, "period": P, "snr": float(depth / err) if err > 0 else 0.0,
                    "depth_ppm": float(depth * 1e6), "period_match": bool(abs(P / period - 1) < 0.001)})
    return out


def periodogram(lc, rec, k, r_star, m_star):
    """The stacked periodogram the search saw for signal k (earlier signals masked, as in search_star)."""
    span = max(np.ptp(lc.time[lc.season == s]) for s in np.unique(lc.season))
    periods = search.period_grid(span, r_star, m_star)
    mask = np.zeros(len(lc.time), bool)
    for prev in rec["signals"][:k - 1]:
        p = prev["detection"]
        mask |= search.transit_mask(lc.time, p["period"], p["t0"], p["duration"], factor=2.0)
    return periods, search.sde(search.stacked_periodogram(lc, periods, mask=mask))


def pm(value, err, max_decimals=8):
    """'value ± err' with the error to two significant figures."""
    if err is None or not np.isfinite(err) or err <= 0:
        return f"{value:.6f}"
    dec = int(np.clip(1 - np.floor(np.log10(err)), 0, max_decimals))
    return f"{value:.{dec}f} ± {err:.{dec}f}"


def centroid_sentence(v):
    """What the light-curve centroid test can and cannot say."""
    z, p = v.get("centroid_z", float("nan")), v.get("centroid_p", float("nan"))
    shift_z = v.get("centroid_shift_z", float("nan"))
    if np.isfinite(shift_z) and shift_z > 0:
        return f"centroid motion in transit: chi2 z = {z:.1f} vs fake epochs (p = {p:.2f}), shift z = {shift_z:.1f}"
    return (f"centroid motion in transit: chi2 z = {z:.1f} vs fake epochs (p = {p:.2f}); the shift itself is "
            f"below the noise, which for a dip this shallow does not exclude a neighbour (see pixel localization)")


def dossier(i):
    """One page per candidate with everything a CTOI submission (ExoFOP) asks for."""
    f, v, g = i["fit"], i["vet"], i["gaia"]
    ok = "error" not in f
    lines = [f"# TIC {i['tic']}, signal {i['signal']}", "",
             f"**Verdict:** {i['verdict']} ({i['novelty']})" + (f"; flags: {'; '.join(i['reasons'])}" if i["reasons"] else ""), "",
             "## Star", f"TESS mag {i['tmag']:.2f}, Teff {i['teff']:.0f} K, R* {i['r_star']:.3f} R_sun, "
             f"M* {i['m_star']:.3f} M_sun (TIC v8.2). RA {i['ra']:.5f}, Dec {i['dec']:.5f}. "
             f"{len(i['sectors'])} sectors: {', '.join(map(str, i['sectors']))}.", ""]
    if ok:
        lines += ["## Transit fit (batman, quadratic limb darkening fixed; stellar-density prior on a/R*)",
                  "| parameter | value |", "|---|---|",
                  f"| period (d) | {pm(f['period'], f['period_err'])} |",
                  f"| mid-transit (BJD_TDB) | {pm(f['t0_bjd'], f['t0_err'])} |",
                  f"| depth (ppm) | {f['depth_ppm']:.0f} |",
                  f"| Rp/R* | {f['rp_rs']:.4f} ± {f['rp_rs_err']:.4f} |",
                  f"| planet radius (R_earth) | {f['rp_rearth']:.2f} ± {f['rp_rearth_err']:.2f} (stellar radius error not included) |",
                  f"| impact parameter | {f['b']:.2f} ± {f['b_err']:.2f} |",
                  f"| a/R* | {f['a_rs']:.1f} ± {f['a_rs_err']:.1f} |",
                  f"| duration T14 (h) | {f['t14_h']:.2f} |",
                  f"| reduced chi-square | {f['chi2_reduced']:.3f} ({f['n_points']} points) |", ""]
    lines += ["## Vetting", f"SNR {v.get('snr_red', float('nan')):.1f} (red-noise aware), folded red-noise SNR "
              f"{v.get('red_snr', float('nan')):.1f}, {v.get('n_transits_measured')} transits observed. "
              f"Odd/even difference {v.get('oddeven_sigma', float('nan')):.1f} sigma; phase-0.5 depth "
              f"{v.get('phase05_depth_ppm', float('nan')):.0f} ppm ({v.get('phase05_sigma', float('nan')):.1f} sigma); "
              f"{centroid_sentence(v)}; largest single-transit share "
              f"{v.get('max_single_frac', float('nan')):.2f}; depth in SAP flux {v.get('sap_depth_ppm', float('nan')):.0f} ppm.",
              "Independent searches of each half of the data: " + ", ".join(
                  f"{h['half']} SNR {h['snr']:.1f}" + (" (same period)" if h.get("period_match") else " (different period)")
                  for h in i["halves"] if h.get("snr") is not None) + ".", ""]
    if "error" not in g:
        lines += ["## Neighbouring stars (Gaia DR3, within 63\")",
                  f"{g['n_gaia']} Gaia sources; target RUWE {g.get('target_ruwe')}; neighbours bright enough to mimic the "
                  f"dip if they were eclipsing binaries: {g['possible_sources']}.", ""]
        risky = [n for n in g["neighbours"] if n["could_be_source"]]
        if risky:
            lines += ["| Gaia DR3 | separation (\") | Δmag | eclipse depth it would need |", "|---|---|---|---|"]
            lines += [f"| {n['source_id']} | {n['sep_arcsec']:.1f} | {n['delta_mag']:.2f} | {100 * n['needed_depth']:.1f}% |"
                      for n in sorted(risky, key=lambda n: n["needed_depth"])]
            lines += [""]
    lines += ["## Catalogue matches", ", ".join(f"{h['name']} ({h['relation']} period, {h['disposition']})"
                                                for h in i["matches"]) or "none", "",
              "## What would confirm or refute it",
              "Ground-based photometry of the predicted transits (to see which star dims), high-resolution imaging, "
              "and radial velocities or statistical validation (e.g. TRICERATOPS).", ""]
    return "\n".join(lines)


def main():
    sig = pd.read_csv(RESULTS / "candidates.csv")
    sig = sig[sig.novelty.isin(["new", "SPOC TCE only (never promoted)"])]
    targets = pd.read_csv(DATA / "targets.csv").set_index("tic")
    tic_tab = pd.read_parquet(DATA / "tic_15plus.parquet").set_index("tic")
    known = crossmatch.load_known(targets.index.tolist())
    rows = []
    for tic, group in sig.groupby("tic"):
        row = targets.loc[tic]
        r_star = float(row.rad)
        m_star = float(row["mass"]) if np.isfinite(row.get("mass", np.nan)) else r_star
        lc = lightcurve.prepare(int(tic))
        rec = json.loads((RESULTS / "search" / f"{tic}.json").read_text())
        for _, s in group.iterrows():
            d = rec["signals"][int(s.signal) - 1]["detection"]
            v = vetting.vet(lc, d, r_star, m_star)
            verdict, reasons = vetting.classify(v)
            hits = crossmatch.match(tic, d["period"], known)
            t = tic_tab.loc[tic]
            try:
                gaia = followup.neighbour_check(float(t.ra), float(t.dec), v["depth_ppm"] * 1e-6, t.get("GAIA"))
            except Exception as exc:  # archive hiccups should not stop the rest
                gaia = {"error": repr(exc)}
            halves = half_search(lc, d["period"], r_star, m_star)
            try:
                fit = followup.fit_transit(lc, d["period"], d["t0"], d["duration"], d["depth"], r_star, m_star)
            except Exception as exc:
                fit = {"error": repr(exc)}
            png = OUT / f"TIC{tic}_{int(s.signal)}.png"
            periods, sde_curve = periodogram(lc, rec, int(s.signal), r_star, m_star)
            report.plot_detection(lc, d, v, verdict, reasons, hits, periods=periods, sde_curve=sde_curve,
                                  path=png, title_extra=f"signal {int(s.signal)}: {verdict} ({s.novelty})")
            info = {"tic": int(tic), "signal": int(s.signal), "verdict": verdict, "reasons": reasons,
                    "novelty": s.novelty, "matches": hits, "vet": v, "gaia": gaia, "halves": halves, "fit": fit,
                    "tmag": float(row.Tmag), "teff": float(row.Teff), "r_star": r_star, "m_star": m_star,
                    "ra": float(t.ra), "dec": float(t.dec), "sectors": rec["sectors"]}
            (OUT / f"TIC{tic}_{int(s.signal)}.json").write_text(json.dumps(info, default=float, indent=1))
            (OUT / f"TIC{tic}_{int(s.signal)}.md").write_text(dossier(info))
            both = all(h.get("snr") and h["snr"] > 3 and h.get("period_match") for h in halves)
            rows.append({"tic": int(tic), "signal": int(s.signal), "period": d["period"], "depth_ppm": v["depth_ppm"],
                         "rp_rearth": v["rp_rearth"], "snr_red": v.get("snr_red"), "verdict": verdict,
                         "novelty": s.novelty, "gaia_possible_sources": gaia.get("possible_sources"),
                         "ruwe": gaia.get("target_ruwe"), "in_both_halves": both,
                         "half_snrs": [h.get("snr") for h in halves], "reasons": "; ".join(reasons)})
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False)


if __name__ == "__main__":
    main()
