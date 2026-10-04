"""Step 4: search, vet, crossmatch and plot every downloaded target (resumable, parallel).

    .venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/04_search.py [--workers 4] [--tics 1,2,3] [--out results]

Writes <out>/search/<tic>.json per star and <out>/plots/TIC<tic>_<n>.png per signal.
"""
import argparse
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for var in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(var, "1")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from tess_search import DATA, RESULTS

_known = None


def init_worker(tics):
    global _known
    from tess_search import crossmatch
    _known = crossmatch.load_known(tics)


def run_one(row, out_dir):
    from tess_search import pipeline

    tic = int(row["tic"])
    if (out_dir / "search" / f"{tic}.json").exists():
        return tic, "skip", 0.0
    t0 = time.time()
    try:
        rec = pipeline.process_star(row, known=_known, out_dir=out_dir)
        pipeline.save(rec, out_dir)
        return tic, "ok", time.time() - t0
    except Exception:
        with open(out_dir / "search_errors.log", "a") as fh:
            fh.write(f"TIC {tic}\n{traceback.format_exc()}\n")
        return tic, "error", time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--tics", type=str, default=None)
    ap.add_argument("--out", type=str, default=str(RESULTS))
    args = ap.parse_args()
    out_dir = Path(args.out)
    (out_dir / "search").mkdir(parents=True, exist_ok=True)
    targets = pd.read_csv(DATA / "targets.csv")
    all_tics = targets.tic.tolist()
    if args.tics:
        wanted = {int(x) for x in args.tics.split(",")}
        targets = targets[targets.tic.isin(wanted)]
    have = {int(p.stem) for p in (DATA / "lc").glob("*.npz") if ".part" not in p.name}
    targets = targets[targets.tic.isin(have)]
    rows = targets.to_dict("records")
    print(f"{len(rows)} downloaded targets to process with {args.workers} workers -> {out_dir}", flush=True)
    t_start = time.time()
    done = 0
    with ProcessPoolExecutor(args.workers, initializer=init_worker, initargs=(all_tics,)) as pool:
        futures = [pool.submit(run_one, r, out_dir) for r in rows]
        for fut in as_completed(futures):
            tic, status, dt = fut.result()
            done += 1
            if status != "skip":
                rate = done / (time.time() - t_start) * 3600
                print(f"[{done}/{len(rows)}] TIC {tic} {status} {dt:.0f}s ({rate:.0f}/h)", flush=True)


if __name__ == "__main__":
    main()
