"""Step 3: download every catalogue used to decide whether a signal is already known.

    .venv/bin/python scripts/03_catalogs.py

Into data/catalogs/:
  toi.csv          TESS Objects of Interest (ExoFOP)
  ctoi.csv         Community TOIs (ExoFOP)
  confirmed.csv    NASA Exoplanet Archive confirmed planets (pscomppars)
  tce/*.csv        every SPOC Threshold Crossing Event file (single- and multi-sector)
  tess_ebs.csv     TESS eclipsing binaries (Prsa et al. 2022), fetched on first use
"""
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tess_search import DATA

CAT = DATA / "catalogs"
(CAT / "tce").mkdir(parents=True, exist_ok=True)

SOURCES = {
    "toi.csv": "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv",
    "ctoi.csv": "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv",
    "confirmed.csv": ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name,hostname,tic_id,"
                      "pl_orbper,pl_rade,pl_tranmid,disc_facility,disc_year,tran_flag+from+pscomppars&format=csv"),
}


def curl(url, dest):
    subprocess.run(["curl", "-s", "-f", "-L", "-o", str(dest), url], check=True)


for name, url in SOURCES.items():
    curl(url, CAT / name)
    print(f"{name}: {sum(1 for _ in open(CAT / name)) - 1} rows")

page = subprocess.run(["curl", "-s", "https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_tce.html"],
                      capture_output=True, text=True, check=True).stdout
links = sorted(set(re.findall(r'href="([^"]*tce[^"]*\.csv)"', page)))
for link in links:
    dest = CAT / "tce" / Path(link).name
    if not dest.exists():
        curl("https://archive.stsci.edu" + link, dest)
print(f"TCE files: {len(links)}")

from tess_search import crossmatch  # noqa: E402

crossmatch.load_ebs()
cache = CAT / "tce_targets.parquet"
if cache.exists():
    cache.unlink()  # rebuilt for the current target list on next use
print("done")
