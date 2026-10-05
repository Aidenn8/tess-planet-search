"""Step 10: download target pixel files for the candidates and the localization test stars.

    .venv/bin/python scripts/10_download_tpfs.py            # candidates + validation set
    .venv/bin/python scripts/10_download_tpfs.py 123 456    # specific TIC IDs
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from tess_search import DATA, pixels

CANDIDATES = [32090583, 229689348, 198412174, 149390648, 294053492]
# Localization tests. On target: confirmed planets spanning the candidates' signal strengths
# (TOI-6000 b is an ultra-short period like TIC 229689348). Off target: TOIs that TFOP retired
# as nearby eclipsing binaries (TOI-419.01 "may be on TIC 279251647"; TOI-2084.02 "NEB"; the
# same star also hosts the confirmed on-target planet TOI-2084 b) and TOI-2283.01 ("possible
# offset towards another star to the north").
VALIDATION = [259233660, 219875976, 150428135, 356016119, 364074068, 441738827, 307210830,
              279251651, 198211976]

if __name__ == "__main__":
    tics = [int(a) for a in sys.argv[1:]] or CANDIDATES + VALIDATION
    index = pd.read_parquet(DATA / "lc_index.parquet")
    for tic in tics:
        print(pixels.download_tpfs(tic, index), flush=True)
