"""Step 19: which neighbouring stars could still produce each candidate's signal?

Seeing-limited follow-up (TFOP SG1) clears a candidate by showing that no nearby star has an eclipse deep
enough to produce the TESS signal. For every Gaia DR3 star in the localization field of each candidate
(results/hardening/localize/TIC<tic>.json, about 100 arcsec around the target), this lists its separation
and position angle from the target, its magnitudes, the eclipse depth it would need to produce the observed
signal (from the joint PRF fit with the source placed on that star), and how strongly the TESS pixels
already exclude it (including the 1.5 arcsec systematic floor).

    .venv/bin/python scripts/19_neighbour_checklist.py

Writes results/followup_planning/neighbours.csv. A star "could cause it" if the eclipse depth it needs is
positive and below 100%; it is "open" if, in addition, the pixels exclude it at less than 3 sigma.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from tess_search import RESULTS

OUT = RESULTS / "followup_planning"

if __name__ == "__main__":
    cs = pd.read_csv(RESULTS / "candidates" / "summary.csv")
    cs = cs[cs.verdict.str.startswith("Candidate")]
    rows = []
    for _, c in cs.iterrows():
        loc = json.loads((RESULTS / "hardening" / "localize" / f"TIC{int(c.tic)}.json").read_text())
        for s in loc["stars"]:
            if s["is_target"]:
                continue
            depth = s["implied_eclipse_depth"]
            could = bool(0 < depth < 1)
            rows.append({"candidate_tic": int(c.tic), "gaia_dr3": str(s["source_id"]),
                         "sep_arcsec": round(float(np.hypot(s["east"], s["north"])), 1),
                         "pa_deg": round(float(np.degrees(np.arctan2(s["east"], s["north"])) % 360), 0),
                         "gmag": round(float(s["gmag"]), 2),
                         "tmag": round(float(s["tmag"]), 2) if s.get("tmag") is not None else None,
                         "eclipse_depth_needed_pct": float(f"{100 * depth:.3g}") if could else None,
                         "pixel_exclusion_sigma": round(float(s["excluded_sigma_total"]), 1),
                         "could_cause_signal": could,
                         "open": could and s["excluded_sigma_total"] < 3})
    df = pd.DataFrame(rows).sort_values(["candidate_tic", "could_cause_signal", "pixel_exclusion_sigma"],
                                        ascending=[True, False, True])
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "neighbours.csv", index=False)
    for tic, g in df.groupby("candidate_tic", sort=False):
        could = g[g.could_cause_signal]
        print(f"TIC {tic}: {len(g)} neighbours, {len(could)} could cause the signal, {int(g.open.sum())} open")
        print(could.head(4).to_string(index=False))
