"""A transit search of heavily observed TESS M dwarfs.

Modules, in pipeline order:
    targets     which stars to search (MAST download lists + TESS Input Catalog)
    download    fetch SPOC 2-minute light curves from MAST
    lightcurve  load, clean, stitch and detrend all sectors of one star
    search      iterative box least squares (BLS) search for periodic dips
    vetting     tests that separate planet-like signals from false alarms
    crossmatch  compare against known planets, TOIs, CTOIs and eclipsing binaries
    inject      injection-recovery tests (how sensitive is the search?)
    report      plots and per-candidate summary pages
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
