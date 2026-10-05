"""Step 13a: inputs for TRICERATOPS (statistical validation), written from the main environment.

TRICERATOPS (Giacalone et al. 2021) runs in its own environment (.venv-tri, it pins
numpy 1.x), so this step exports what it needs for each candidate:
  * the folded light curve within 3 transit durations of mid-transit, in 2-minute
    phase bins (flux normalised, PDCSAP so already corrected for known contamination)
  * the period and depth
  * the SPOC photometric aperture (CCD column, row of each pixel) for 4 sectors spread
    over the mission, which TRICERATOPS uses to work out how much light each
    neighbouring star puts in the aperture

    .venv/bin/python scripts/13_triceratops_inputs.py
then
    .venv-tri/bin/python scripts/13_triceratops_run.py
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from tess_search import RESULTS, lightcurve, pixels

OUT = RESULTS / "hardening" / "triceratops"
CANDIDATE_FILES = ["TIC32090583_5", "TIC229689348_1", "TIC198412174_1", "TIC149390648_1", "TIC294053492_1"]
N_SECTORS = 4


def folded(lc, period, t0, t14, bin_min=2.0):
    dt = ((lc.time - t0 + 0.5 * period) % period) - 0.5 * period
    m = np.abs(dt) < 3 * t14
    edges = np.arange(-3 * t14, 3 * t14 + 1e-9, bin_min / 1440)
    idx = np.clip(np.digitize(dt[m], edges) - 1, 0, len(edges) - 2)
    n = np.bincount(idx, minlength=len(edges) - 1)
    f = np.bincount(idx, weights=lc.flux[m], minlength=len(edges) - 1)
    ok = n > 0
    centres = 0.5 * (edges[1:] + edges[:-1])
    flux = f[ok] / n[ok]
    # one error for all bins, as TRICERATOPS takes: out-of-transit scatter of the binned points
    out = np.abs(centres[ok]) > 0.75 * t14
    return centres[ok], flux, float(np.std(flux[out]))


def apertures(tic):
    stamps = pixels.load_stamps(tic)
    pick = np.unique(np.linspace(0, len(stamps) - 1, N_SECTORS).round().astype(int))
    out = {}
    for i in pick:
        st = stamps[i]
        rows, cols = np.nonzero(st.aperture)
        out[int(st.sector)] = np.column_stack([st.col0 + cols, st.row0 + rows]).tolist()
    return out


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for f in CANDIDATE_FILES:
        d = json.loads((RESULTS / "followup" / f"{f}.json").read_text())
        fit = d["fit"]
        lc = lightcurve.prepare(d["tic"])
        t14 = fit["t14_h"] / 24
        time, flux, err = folded(lc, fit["period"], fit["t0_btjd"], t14)
        rec = {"tic": d["tic"], "label": f"TIC{d['tic']}", "ra": d["ra"], "dec": d["dec"],
               "period": fit["period"], "depth": fit["depth_ppm"] * 1e-6,
               "t14": t14, "time": time.tolist(), "flux": flux.tolist(), "flux_err": err,
               "apertures": apertures(d["tic"])}
        (OUT / f"TIC{d['tic']}_input.json").write_text(json.dumps(rec))
        print(f"TIC {d['tic']}: {len(time)} bins, err {err * 1e6:.0f} ppm, sectors {list(rec['apertures'])}", flush=True)
