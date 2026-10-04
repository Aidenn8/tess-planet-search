"""Download SPOC 2-minute light curves and keep only what the search needs.

Each star's sectors are fetched from MAST, the useful columns are copied into one
compressed file per star (data/lc/<tic>.npz), and the original FITS files are
deleted (they can always be downloaded again). A star whose .npz exists is
skipped, so the job can be stopped and restarted at any time.
"""
import json
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits

from . import DATA
from .targets import URL_PREFIX

LC_DIR = DATA / "lc"
TMP_DIR = DATA / "tmp_fits"
FAIL_LOG = DATA / "download_failures.jsonl"

COLUMNS = {  # FITS column -> (key in npz, dtype)
    "TIME": ("time", np.float64),
    "PDCSAP_FLUX": ("pdcsap", np.float32),
    "PDCSAP_FLUX_ERR": ("pdcsap_err", np.float32),
    "SAP_FLUX": ("sap", np.float32),
    "SAP_BKG": ("bkg", np.float32),
    "QUALITY": ("quality", np.int32),
    "MOM_CENTR1": ("centr1", np.float32),
    "MOM_CENTR2": ("centr2", np.float32),
}
SECTOR_KEYS = ["CROWDSAP", "FLFRCSAP", "CDPP0_5", "CDPP1_0", "CDPP2_0", "PDCMETHD"]


def fetch(filename, dest, retries=5):
    url = URL_PREFIX + filename
    for attempt in range(retries):
        r = subprocess.run(
            ["curl", "-s", "-f", "-L", "--max-time", "300", "-o", str(dest), url],
            capture_output=True,
        )
        if r.returncode == 0 and dest.exists() and dest.stat().st_size > 10_000:
            return
        time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"curl failed for {filename} (code {r.returncode})")


def read_sector(path):
    with fits.open(path, memmap=False) as hdul:
        h0, h1 = hdul[0].header, hdul[1].header
        d = hdul[1].data
        ok = np.isfinite(d["TIME"])
        arrays = {key: np.asarray(d[col][ok], dtype=dt) for col, (key, dt) in COLUMNS.items()}
        meta = {"sector": int(h0["SECTOR"]), "camera": int(h0["CAMERA"]), "ccd": int(h0["CCD"])}
        for k in SECTOR_KEYS:
            meta[k.lower()] = h1.get(k)
        meta["n_aperture_pix"] = int(np.sum((hdul[2].data & 2) > 0))
        star = {k.lower(): h0.get(k) for k in ["TICID", "TESSMAG", "TEFF", "RADIUS", "RA_OBJ", "DEC_OBJ"]}
    return arrays, meta, star


def download_star(tic, filenames, threads=4):
    """Fetch all sectors of one star and write data/lc/<tic>.npz. Returns status string."""
    out = LC_DIR / f"{tic}.npz"
    if out.exists():
        return "skip"
    tmp = TMP_DIR / str(tic)
    tmp.mkdir(parents=True, exist_ok=True)
    paths = [tmp / f for f in filenames]
    with ThreadPoolExecutor(threads) as pool:
        list(pool.map(lambda fp: fetch(fp[0], fp[1]) if not fp[1].exists() else None,
                      zip(filenames, paths)))
    parts, metas, star = [], [], None
    for p in paths:
        arrays, meta, star = read_sector(p)
        arrays["sector"] = np.full(len(arrays["time"]), meta["sector"], dtype=np.int16)
        parts.append(arrays)
        metas.append(meta)
    merged = {k: np.concatenate([a[k] for a in parts]) for k in parts[0]}
    order = np.argsort(merged["time"], kind="stable")
    merged = {k: v[order] for k, v in merged.items()}
    merged["sector_meta"] = np.array(json.dumps(sorted(metas, key=lambda m: m["sector"])))
    merged["star_meta"] = np.array(json.dumps(star))
    part = out.with_suffix(".part.npz")
    np.savez_compressed(part, **merged)
    part.rename(out)
    shutil.rmtree(tmp, ignore_errors=True)
    return "ok"


def load_star(tic):
    """Load a downloaded star: dict of arrays plus parsed sector/star metadata."""
    with np.load(LC_DIR / f"{tic}.npz") as z:
        d = {k: z[k] for k in z.files}
    d["sector_meta"] = json.loads(str(d["sector_meta"]))
    d["star_meta"] = json.loads(str(d["star_meta"]))
    return d


def download_many(index, tics, star_workers=2, threads=4):
    LC_DIR.mkdir(parents=True, exist_ok=True)
    groups = index[index.tic.isin(tics)].groupby("tic").filename.apply(list)
    todo = [t for t in tics if t in groups.index]

    def work(tic):
        try:
            return tic, download_star(tic, groups[tic], threads=threads)
        except Exception as exc:
            with open(FAIL_LOG, "a") as fh:
                fh.write(json.dumps({"tic": int(tic), "error": repr(exc), "t": time.time()}) + "\n")
            return tic, "fail"

    done = 0
    t0 = time.time()
    with ThreadPoolExecutor(star_workers) as pool:
        for tic, status in pool.map(work, todo):
            done += 1
            if status != "skip" or done % 100 == 0:
                rate = done / (time.time() - t0)
                print(f"[{done}/{len(todo)}] TIC {tic}: {status}  ({rate * 3600:.0f} stars/h)", flush=True)
