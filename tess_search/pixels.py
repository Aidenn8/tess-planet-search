"""Target pixel files (TPFs): the small postage-stamp images behind each light curve.

A light curve says *that* something dimmed; the pixels say *where*. TESS pixels
are 21" wide, so a neighbouring eclipsing binary a few pixels away can leak a
planet-sized dip into the target's aperture. Comparing images taken during
transit with images taken just before and after shows which part of the stamp
lost light (`localize.py` turns that into a position on the sky).

This module only downloads and reads the files. SPOC TPFs share their names with
the light curves (`..._lc.fits` -> `..._tp.fits`); they are kept as FITS under
data/tpf/<tic>/ because the sky-to-pixel mapping (WCS) lives in their headers.
"""
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

from . import DATA
from .targets import URL_PREFIX

TPF_DIR = DATA / "tpf"


def tpf_filenames(tic, index=None):
    index = index if index is not None else pd.read_parquet(DATA / "lc_index.parquet")
    rows = index[index.tic == int(tic)].sort_values("sector")
    return [(int(s), f.replace("_lc.fits", "_tp.fits")) for s, f in zip(rows.sector, rows.filename)]


def _fetch(filename, dest, retries=5):
    tmp = dest.with_suffix(".part")
    for attempt in range(retries):
        r = subprocess.run(["curl", "-s", "-f", "-L", "--max-time", "600", "-o", str(tmp),
                            URL_PREFIX + filename], capture_output=True)
        if r.returncode == 0 and tmp.exists() and tmp.stat().st_size > 100_000:
            tmp.rename(dest)
            return "ok"
        time.sleep(3 * (attempt + 1))
    return f"fail ({r.returncode})"


def download_tpfs(tic, index=None, threads=4):
    """Fetch every sector's TPF for one star (skips files already present)."""
    out = TPF_DIR / str(int(tic))
    out.mkdir(parents=True, exist_ok=True)
    todo = [(f, out / f) for _, f in tpf_filenames(tic, index) if not (out / f).exists()]
    with ThreadPoolExecutor(threads) as pool:
        status = list(pool.map(lambda fd: _fetch(*fd), todo))
    return {"tic": int(tic), "fetched": status.count("ok"), "failed": len(status) - status.count("ok"),
            "have": len(list(out.glob("*_tp.fits")))}


@dataclass
class Stamp:
    """One sector of pixel data, cleaned to good cadences."""
    tic: int
    sector: int
    camera: int
    ccd: int
    time: np.ndarray        # BTJD
    flux: np.ndarray        # (n_cadence, ny, nx), background-subtracted e-/s
    aperture: np.ndarray    # (ny, nx) bool, SPOC photometric aperture
    wcs: WCS
    col0: int               # CCD column of pixel [.., 0]
    row0: int               # CCD row of pixel [0, ..]

    def sky_to_pix(self, ra, dec):
        """Sky position(s) -> (x, y) in stamp pixel coordinates (0-based, x = column)."""
        x, y = self.wcs.all_world2pix(np.atleast_1d(ra), np.atleast_1d(dec), 0)
        return x, y

    def pix_to_sky(self, x, y):
        ra, dec = self.wcs.all_pix2world(np.atleast_1d(x), np.atleast_1d(y), 0)
        return ra, dec


def load_stamps(tic):
    """All sectors of a star as Stamp objects, time-sorted."""
    stamps = []
    for path in sorted((TPF_DIR / str(int(tic))).glob("*_tp.fits")):
        with fits.open(path, memmap=False) as h:
            h0, h1 = h[0].header, h[1].header
            d = h[1].data
            good = (d["QUALITY"] == 0) & np.isfinite(d["TIME"])
            flux = np.asarray(d["FLUX"][good], dtype=np.float32)
            good2 = np.isfinite(flux).all(axis=(1, 2))
            # WCS keywords are stored per column in the binary table; astropy reads them via keysel
            wcs = WCS(h[2].header)
            stamps.append(Stamp(
                tic=int(tic), sector=int(h0["SECTOR"]), camera=int(h0["CAMERA"]), ccd=int(h0["CCD"]),
                time=np.asarray(d["TIME"][good][good2], dtype=np.float64), flux=flux[good2],
                aperture=(np.asarray(h[2].data) & 2) > 0, wcs=wcs,
                col0=int(h1["1CRV4P"]), row0=int(h1["2CRV4P"]),
            ))
    return stamps
