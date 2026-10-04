"""Target selection.

The MAST bulk-download scripts (one per sector) list every SPOC 2-minute light
curve file. Counting sectors per TIC ID tells us which stars TESS has watched
longest; the TESS Input Catalog (TIC v8.2) gives their temperature, radius and
brightness so we can keep the small, cool, bright ones.
"""
import glob
import re
import time

import numpy as np
import pandas as pd

from . import DATA

FILE_RE = re.compile(r"-o (tess\d+-s(\d{4})-(\d{16})-\d{4}-s_lc\.fits)")
URL_PREFIX = "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/"

TIC_COLUMNS = [
    "ID", "ra", "dec", "Tmag", "Teff", "e_Teff", "logg", "rad", "e_rad", "mass",
    "rho", "lumclass", "contratio", "disposition", "duplicate_id", "d", "plx",
    "GAIA", "gaiaqflag", "ebv",
]


def build_lc_index(script_glob=None):
    """Parse every tesscurl_sector_N_lc.sh into a (sector, tic, filename) table."""
    script_glob = script_glob or str(DATA / "scripts" / "tesscurl_sector_*_lc.sh")
    rows = []
    for path in glob.glob(script_glob):
        with open(path) as fh:
            for line in fh:
                m = FILE_RE.search(line)
                if m:
                    rows.append((int(m.group(2)), int(m.group(3)), m.group(1)))
    df = pd.DataFrame(rows, columns=["sector", "tic", "filename"])
    return df.drop_duplicates().sort_values(["tic", "sector"]).reset_index(drop=True)


def sector_counts(index):
    return index.groupby("tic").sector.nunique().rename("n_sectors")


def query_tic(tic_ids, chunk=400, pause=0.5, retries=4):
    """Fetch TIC v8.2 rows for many TIC IDs from MAST, in chunks."""
    from astroquery.mast import Catalogs

    tic_ids = [int(t) for t in tic_ids]
    out = []
    for i in range(0, len(tic_ids), chunk):
        ids = tic_ids[i : i + chunk]
        for attempt in range(retries):
            try:
                tab = Catalogs.query_criteria(catalog="Tic", ID=ids)
                break
            except Exception as exc:  # network hiccups: back off and retry
                if attempt == retries - 1:
                    raise
                time.sleep(5 * (attempt + 1))
                print(f"TIC query retry {attempt + 1} after {exc!r}")
        df = tab.to_pandas()
        keep = [c for c in TIC_COLUMNS if c in df.columns]
        out.append(df[keep])
        time.sleep(pause)
        print(f"TIC: {min(i + chunk, len(tic_ids))}/{len(tic_ids)}", flush=True)
    tic = pd.concat(out, ignore_index=True)
    tic["ID"] = tic["ID"].astype(np.int64)
    return tic.rename(columns={"ID": "tic"}).drop_duplicates("tic")


def select_m_dwarfs(tic, counts, min_sectors=20, teff_max=3900.0, rad_max=0.65,
                    tmag_max=13.5, contratio_max=0.2):
    """Keep bright, uncrowded, main-sequence M dwarfs with many sectors.

    teff_max 3900 K and rad_max 0.65 R_sun roughly bound the M spectral class;
    contratio is the TIC's estimate of flux from neighbours divided by the
    target's flux in the TESS aperture (0.2 = neighbours add 20%).
    """
    df = tic.merge(counts.reset_index(), on="tic")
    keep = (
        (df.n_sectors >= min_sectors)
        & (df.Teff <= teff_max)
        & (df.rad <= rad_max)
        & (df.Tmag <= tmag_max)
        & (df.lumclass == "DWARF")
        & ((df.contratio.isna()) | (df.contratio <= contratio_max))
        & (df.disposition.isna() | (df.disposition == ""))  # drop TIC duplicates/artifacts
    )
    return df[keep].sort_values("Tmag").reset_index(drop=True)
