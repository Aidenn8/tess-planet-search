"""Step 5: how complete and how reliable is the pipeline?

Two experiments, both running the exact pipeline used on the real data:

  inject   plant one synthetic planet (batman model) into a real star's raw
           light curve, run everything, and record whether it came back out
           (detected at the right period and epoch) and whether it survived vetting.
           -> completeness as a function of period and planet size.

  invert   flip each light curve upside down (dips become bumps, bumps become
           dips) and run everything. Real planets cannot appear; whatever passes
           vetting is a false alarm made by noise and systematics.
           -> false-alarm rate per star.

    .venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/05_reliability.py inject --n 300
    .venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/05_reliability.py invert --n 200

Results are appended to results/reliability/<mode>.jsonl (resumable by trial id).
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

OUT = RESULTS / "reliability"
_known = None


def init_worker(tics):
    global _known
    from tess_search import crossmatch
    _known = crossmatch.load_known(tics)


def make_trials(mode, n, seed=42):
    """Deterministic list of trials. Stars with known planets/TOIs are skipped so an
    injected planet is not confused with a real one."""
    from tess_search import crossmatch

    targets = pd.read_csv(DATA / "targets.csv")
    have = {int(p.stem) for p in (DATA / "lc").glob("*.npz") if ".part" not in p.name}
    known = crossmatch.load_known(targets.tic.tolist())
    hosts = set(known[known.source.isin(["TOI", "CTOI", "CONFIRMED", "EB"])].tic)
    pool = targets[targets.tic.isin(have) & ~targets.tic.isin(hosts)].reset_index(drop=True)
    rng = np.random.default_rng(seed)
    trials = []
    for i in range(n):
        row = pool.iloc[int(rng.integers(len(pool)))].to_dict()
        trial = {"id": i, "row": row}
        if mode == "inject":
            trial.update(period=float(np.exp(rng.uniform(np.log(0.5), np.log(40)))),
                         rp=float(rng.uniform(0.6, 4.0)), b=float(rng.uniform(0, 0.9)),
                         phase=float(rng.uniform(0, 1)))
        trials.append(trial)
    return trials


def run_trial(mode, trial):
    from tess_search import download, inject, pipeline

    row = trial["row"]
    raw = download.load_star(int(row["tic"]))
    rec = {"id": trial["id"], "tic": int(row["tic"]), "tmag": row["Tmag"], "r_star": row["rad"]}
    r_star = float(row["rad"])
    m_star = float(row["mass"]) if np.isfinite(row.get("mass", np.nan)) else r_star
    if mode == "inject":
        good = np.isfinite(raw["time"])
        t0 = float(raw["time"][good][0] + trial["phase"] * trial["period"])
        raw = inject.inject(raw, trial["period"], t0, trial["rp"], r_star, m_star, trial["b"])
        depth = float(1 - inject.transit_model(np.array([t0]), trial["period"], t0, trial["rp"],
                                               r_star, m_star, trial["b"])[0])
        rec.update(period=trial["period"], rp=trial["rp"], b=trial["b"], t0=t0, depth_true=depth)
    else:
        raw = dict(raw)
        for key in ("pdcsap", "sap"):
            f = raw[key].astype(np.float64).copy()
            for s in np.unique(raw["sector"]):
                m = raw["sector"] == s
                med = np.nanmedian(f[m])
                f[m] = 2 * med - f[m]
            raw[key] = f.astype(np.float32)
    result = pipeline.process_star(row, known=_known, raw=raw, plots=False)
    rec["noise_ppm"] = result["diagnostics"]["noise_ppm"]
    rec["signals"] = [
        {"period": s["detection"]["period"], "t0": s["detection"]["t0"], "snr": s["detection"]["snr"],
         "snr_red": (s.get("vet") or {}).get("snr_red"), "depth": s["detection"]["depth"],
         "verdict": s.get("verdict", "not vetted"), "reasons": s.get("reasons", [])}
        for s in result["signals"]
    ]
    if mode == "inject":
        rec["recovered"], rec["verdict"] = False, None
        for s in rec["signals"]:
            if abs(s["period"] / rec["period"] - 1) < 0.005:
                n = np.round((s["t0"] - t0) / rec["period"])
                if abs(s["t0"] - t0 - n * rec["period"]) < 0.1:
                    rec["recovered"], rec["verdict"] = True, s["verdict"]
                    rec["snr_found"] = s["snr_red"]
                    break
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["inject", "invert"])
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--workers", type=int, default=3)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{args.mode}.jsonl"
    done = set()
    if path.exists():
        done = {json.loads(l)["id"] for l in path.read_text().splitlines() if l.strip()}
    trials = [t for t in make_trials(args.mode, args.n) if t["id"] not in done]
    targets = pd.read_csv(DATA / "targets.csv")
    print(f"{args.mode}: {len(trials)} trials to run ({len(done)} already done)", flush=True)
    t_start = time.time()
    with ProcessPoolExecutor(args.workers, initializer=init_worker, initargs=(targets.tic.tolist(),)) as pool:
        futs = {pool.submit(run_trial, args.mode, t): t["id"] for t in trials}
        for k, fut in enumerate(as_completed(futs), 1):
            try:
                rec = fut.result()
            except Exception:
                with open(OUT / "errors.log", "a") as fh:
                    fh.write(f"trial {futs[fut]}\n{traceback.format_exc()}\n")
                continue
            with open(path, "a") as fh:
                fh.write(json.dumps(rec, default=float) + "\n")
            rate = k / (time.time() - t_start) * 3600
            extra = f"recovered={rec['recovered']} verdict={rec['verdict']}" if args.mode == "inject" else \
                f"candidates={sum(s['verdict'] == 'candidate' for s in rec['signals'])}"
            print(f"[{k}/{len(trials)}] trial {rec['id']} TIC {rec['tic']} {extra} ({rate:.0f}/h)", flush=True)


if __name__ == "__main__":
    main()
