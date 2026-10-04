"""Step 1: build the light-curve index and choose the M-dwarf sample.

Writes data/lc_index.parquet, data/tic_15plus.parquet and data/targets.csv.
Expects the per-sector MAST scripts in data/scripts/ (see README for the curl line).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from tess_search import DATA, targets

MIN_SECTORS = 20
TMAG_MAX = 13.5

index_path = DATA / "lc_index.parquet"
if index_path.exists():
    index = pd.read_parquet(index_path)
else:
    index = targets.build_lc_index()
    index.to_parquet(index_path)
counts = targets.sector_counts(index)
print(f"{len(index)} light-curve files, {len(counts)} stars, latest sector {index.sector.max()}")

tic_path = DATA / "tic_15plus.parquet"
if tic_path.exists():
    tic = pd.read_parquet(tic_path)
else:
    tic = targets.query_tic(counts[counts >= 15].index)
    tic.to_parquet(tic_path)

sample = targets.select_m_dwarfs(tic, counts, min_sectors=MIN_SECTORS, tmag_max=TMAG_MAX)
sample.to_csv(DATA / "targets.csv", index=False)
print(f"selected {len(sample)} M dwarfs with >= {MIN_SECTORS} sectors and Tmag <= {TMAG_MAX}")
print(sample[["Teff", "rad", "Tmag", "n_sectors"]].describe().round(2))
