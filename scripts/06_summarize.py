"""Step 6: collect per-star results into tables and a summary.

    .venv/bin/python scripts/06_summarize.py [results_dir]

Writes:
  results/all_signals.csv     every vetted signal with its metrics and verdict
  results/candidates.csv      signals that passed vetting (candidate / weak candidate)
  results/summary.json        counts used in the write-up
"""
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

import re

from tess_search import RESULTS, vetting

out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else RESULTS
records = [json.loads(p.read_text()) for p in sorted((out_dir / "search").glob("*.json"))]
print(f"{len(records)} stars processed")

# Re-apply the current classification rules to the stored metrics, so every star is
# judged by the same (final) rules even if the search ran while rules were refined.
for r in records:
    for s in r["signals"]:
        if s.get("vet"):
            s["verdict"], s["reasons"] = vetting.classify(s["vet"])


def reason_category(text):
    """'one transit carries 63% of the signal' -> 'one transit carries N% of the signal'."""
    t = text.split(" (")[0].split(":")[0]
    return re.sub(r"[-+]?\d+(\.\d+)?", "N", t)

rows = []
for r in records:
    for k, s in enumerate(r["signals"], 1):
        d = s["detection"]
        v = s.get("vet") or {}
        rows.append({
            "tic": r["tic"], "signal": k, "tmag": r["tmag"], "teff": r["teff"], "r_star": r["r_star"],
            "n_sectors": r["n_sectors"], "noise_ppm": r["diagnostics"].get("noise_ppm"),
            "period": d["period"], "t0_btjd": d["t0"], "depth_ppm": v.get("depth_ppm", d["depth"] * 1e6),
            "duration_h": (v.get("trap") or {}).get("t14", d["duration"]) * 24,
            "rp_rearth": v.get("rp_rearth"), "snr_bls": d["snr"], "snr_red": v.get("snr_red"),
            "red_snr_fold": v.get("red_snr"), "sde": d["sde"], "n_transits": v.get("n_transits_measured", d["n_transits"]),
            "oddeven_sigma": v.get("oddeven_sigma"), "phase05_sigma": v.get("phase05_sigma"),
            "centroid_shift_z": v.get("centroid_shift_z"), "centroid_offset_arcsec": v.get("centroid_offset_arcsec"),
            "duration_ratio": v.get("duration_ratio"), "max_single_frac": v.get("max_single_frac"),
            "verdict": s.get("verdict", "not vetted"), "reasons": "; ".join(s.get("reasons", [])),
            "novelty": s.get("novelty", ""),
            "matches": "; ".join(f"{h['name']} ({h['relation']})" for h in s.get("matches", [])),
            "plot": s.get("plot", ""),
        })
sig = pd.DataFrame(rows)
sig.to_csv(out_dir / "all_signals.csv", index=False)

passed = sig[sig.verdict.isin(["candidate", "weak candidate"])].sort_values(["verdict", "snr_red"], ascending=[True, False])
passed.to_csv(out_dir / "candidates.csv", index=False)

vetted = sig[sig.verdict != "not vetted"]
summary = {
    "stars": len(records),
    "signals_total": int(len(sig)),
    "signals_vetted": int(len(vetted)),
    "verdicts": dict(Counter(vetted.verdict)),
    "novelty_of_passed": dict(Counter(passed.novelty)),
    "new_candidates": passed[(passed.verdict == "candidate") & (passed.novelty == "new")][
        ["tic", "period", "depth_ppm", "rp_rearth", "snr_red"]].to_dict("records"),
    "new_weak_candidates": int(((passed.verdict == "weak candidate") & (passed.novelty == "new")).sum()),
    "tce_only_candidates": passed[passed.novelty.str.startswith("SPOC TCE")][
        ["tic", "period", "depth_ppm", "rp_rearth", "snr_red", "verdict"]].to_dict("records"),
    "fp_reason_counts": dict(Counter(
        reason_category(r) for rs in vetted[vetted.verdict == "false positive"].reasons for r in rs.split("; ") if r
    ).most_common()),
    "median_runtime_s": float(np.median([r["runtime_s"] for r in records])) if records else None,
}
(out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=float))

print(json.dumps({k: v for k, v in summary.items() if k not in ("new_candidates", "tce_only_candidates")}, indent=2, default=float))
print("\nNEW candidates (passed every test, not in any catalogue):")
cols = ["tic", "signal", "period", "depth_ppm", "rp_rearth", "duration_h", "snr_red", "red_snr_fold", "n_transits", "tmag"]
print(passed[(passed.verdict == "candidate") & (passed.novelty == "new")][cols].round(3).to_string(index=False))
print("\nSPOC TCE only (found by NASA's pipeline, never promoted):")
print(passed[passed.novelty.str.startswith("SPOC TCE")][cols + ["verdict"]].round(3).to_string(index=False))
print("\nNEW weak candidates:")
print(passed[(passed.verdict == "weak candidate") & (passed.novelty == "new")][cols + ["reasons"]].round(3).to_string(index=False, max_colwidth=70))
