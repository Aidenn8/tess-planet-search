"""Print the TIC IDs of targets that host a known planet, TOI or CTOI (and are downloaded).

Used to check that the search rediscovers known signals before trusting it on the rest.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from tess_search import DATA, crossmatch

targets = pd.read_csv(DATA / "targets.csv")
known = crossmatch.load_known(targets.tic.tolist())
hosts = sorted(known[known.source.isin(["TOI", "CTOI", "CONFIRMED"])].tic.unique())
have = {int(p.stem) for p in (DATA / "lc").glob("*.npz") if ".part" not in p.name}
ready = [t for t in hosts if t in have]
print(f"{len(hosts)} known hosts, {len(ready)} downloaded", file=sys.stderr)
print(",".join(str(t) for t in ready))
