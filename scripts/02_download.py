"""Step 2: download every target's light curves (resumable).

Run under the thermal guard:
    .venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/02_download.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from tess_search import DATA, download

index = pd.read_parquet(DATA / "lc_index.parquet")
targets = pd.read_csv(DATA / "targets.csv")
limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
tics = targets.tic.tolist()[:limit]
download.download_many(index, tics, star_workers=3, threads=4)
