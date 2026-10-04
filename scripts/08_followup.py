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
            png = OUT / f"TIC{tic}_{int(s.signal)}.png"
            report.plot_detection(lc, d, v, verdict, reasons, hits, path=png,
                                  title_extra=f"signal {int(s.signal)}: {verdict} ({s.novelty})")
            info = {"tic": int(tic), "signal": int(s.signal), "verdict": verdict, "reasons": reasons,
                    "novelty": s.novelty, "matches": hits, "vet": v, "gaia": gaia, "halves": halves,
                    "tmag": float(row.Tmag), "teff": float(row.Teff), "r_star": r_star}
            (OUT / f"TIC{tic}_{int(s.signal)}.json").write_text(json.dumps(info, default=float, indent=1))
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
