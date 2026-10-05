"""Step 11: why did SPOC's own SNR for two candidates fall as data were added?

For TIC 229689348 and TIC 198412174 the SPOC multi-sector TCE SNR dropped run after
run. Hypothesis: SPOC's period (and its 0.01-d-rounded epoch) drifted slightly
between runs; over ~1,000-4,000 orbits a tiny period error smears a 0.5-h transit
across many hours and the folded SNR collapses even though the signal is unchanged.
Test: coherent box SNR of *our* cleaned light curve as a function of trial period,
best epoch at each period, with SPOC's periods marked.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares

from tess_search import DATA, RESULTS, lightcurve
from tess_search.lightcurve import robust_std

OUT = RESULTS / "hardening"
CASES = {229689348: 0.465377, 198412174: 1.390160, 149390648: 2.836887}

tce = pd.read_parquet(DATA / "catalogs" / "tce_targets.parquet")
summary = {}
fig, axes = plt.subplots(len(CASES), 1, figsize=(8, 2.6 * len(CASES)))
for ax, (tic, p_ours) in zip(axes, CASES.items()):
    lc = lightcurve.prepare(tic)
    span = np.ptp(lc.time)
    runs = tce[(tce.ticid == tic) & (abs(tce.tce_period / p_ours - 1) < 1e-3)]
    lo = min(runs.tce_period.min(), p_ours) - 3e-5
    hi = max(runs.tce_period.max(), p_ours) + 3e-5
    dur = 0.55 / 24 if p_ours < 1 else 0.6 / 24 if p_ours < 2 else 1.1 / 24
    grid = np.arange(lo, hi, p_ours * dur / span / 8)
    bls = BoxLeastSquares(lc.time, lc.flux, dy=np.full(len(lc.time), robust_std(lc.flux - 1)))
    res = bls.power(grid, dur, objective="snr", oversample=20)
    snr = np.asarray(res.depth_snr)
    ax.plot(grid, snr, "k-", lw=0.7)
    ax.axvline(p_ours, color="C2", lw=1.5, label=f"this work {p_ours:.6f} d")
    rows = []
    for _, r in runs.sort_values("run").iterrows():
        j = int(np.argmin(abs(grid - r.tce_period)))
        ax.axvline(r.tce_period, color="C3", ls="--", lw=1)
        ax.text(r.tce_period, snr.max() * 0.95, f" {r.run}\n SPOC SNR {r.tce_model_snr:.1f}", fontsize=7,
                color="C3", va="top")
        drift_h = abs(r.tce_period - p_ours) * span / p_ours * 24
        rows.append({"run": r.run, "spoc_period": r.tce_period, "spoc_snr": r.tce_model_snr,
                     "our_snr_at_spoc_period": float(snr[j]), "drift_over_baseline_h": drift_h})
    k = int(np.argmin(abs(grid - p_ours)))
    summary[tic] = {"our_period": p_ours, "our_snr_box": float(snr[k]), "baseline_d": span, "runs": rows}
    ax.set_title(f"TIC {tic}: box SNR of this work's light curve vs trial period", fontsize=9, loc="left")
    ax.set_ylabel("SNR")
    ax.legend(fontsize=7, loc="upper left")
axes[-1].set_xlabel("trial period (d)")
fig.tight_layout()
fig.savefig(OUT / "spoc_period_check.png", dpi=130)
(OUT / "spoc_period_check.json").write_text(json.dumps(summary, indent=1, default=float))
print(json.dumps(summary, indent=1, default=float))
