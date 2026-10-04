# TESS M-dwarf deep search

A search for transiting planets around the **1,279 bright M dwarfs that NASA's TESS
telescope has watched the longest** (20 to 44 sectors of 2-minute data each, through
sector 107), with automated vetting, comparison against every public catalogue,
and measured completeness and false-alarm rates.

> **Status:** pipeline complete and validated; results are in `results/` and summarised
> in [`REPORT.md`](REPORT.md). Plain-language overview: [`EXPLAINER.md`](EXPLAINER.md).

## Why these stars

* **M dwarfs are small**, so an Earth-sized planet blocks ~0.1% of their light, which is
  10x more than the same planet in front of the Sun. That is the easiest place to find
  small planets, including temperate ones.
* **Many sectors of data** means many transits even for long orbits, so smaller and
  longer-period planets rise above the noise.
* **New data**: NASA's latest combined multi-sector search (SPOC, sectors 1-96, June
  2026) predates sectors 97-107. This search uses all of them.

## Pipeline

| step | script | what it does |
|---|---|---|
| 1 | `scripts/01_select_targets.py` | parse MAST's per-sector download lists (1.7 M light curves), count sectors per star, query the TESS Input Catalog, keep M dwarfs (Teff <= 3900 K, R <= 0.65 R_sun, Tmag <= 13.5, >= 20 sectors, low contamination) |
| 2 | `scripts/02_download.py` | download SPOC 2-min light curves, keep the needed columns in one compressed file per star |
| 4 | `scripts/04_search.py` | clean, detrend, search, vet, crossmatch and plot every star |
| 5 | `scripts/05_reliability.py` | injection-recovery (completeness) and inverted light curves (false alarms) |
| 6 | `scripts/06_summarize.py` | candidate tables and counts |
| 7 | `scripts/07_figures.py` | figures for the report |

### Cleaning and detrending (`tess_search/lightcurve.py`)
Good-quality cadences only (QUALITY = 0), each sector normalised, flares removed (two or
more consecutive cadences above 3 sigma, or one above 5 sigma, plus a short decay tail),
then a robust biweight filter (wotan) with a 0.5-day window, shortened for strongly
spotted fast rotators so rotation does not leak into the search.

### Search (`tess_search/search.py`)
Eight years of data with big gaps would need millions of trial periods for a single
coherent search. Instead each **season** (a run of sectors without long gaps) is searched
separately with box least squares on a shared, physically motivated period grid
(Ofir 2014; 0.4 to 40 days), the season log-likelihoods are **added**, and each peak is
then refit coherently on all data. The search repeats after masking each signal (up to
5 per star) to find multi-planet systems.

### Vetting (`tess_search/vetting.py`)
Odd/even depths, secondary eclipses, transit shape, implied radius, duration versus the
star's density, per-transit consistency, single-event dominance, data-gap pile-up,
spacecraft-orbit and stellar-rotation periods, centroid motion and background changes.
TESS noise is correlated, so flux tests use noise measured at the transit timescale, and
centroid/background tests are **calibrated against fake transit epochs** rather than
assumed error bars.

### Crossmatch (`tess_search/crossmatch.py`)
TOIs, community TOIs, confirmed planets (NASA Exoplanet Archive), the TESS eclipsing-binary
catalogue, and **every SPOC Threshold Crossing Event** (132 single- and multi-sector
runs), including signals NASA's pipeline found but never promoted.

## Validation on known planets
On the stars hosting known planets, the pipeline independently recovers every transiting
planet in its period range, including all four TOI-700 planets (TOI-700 d and e are
Earth-sized, in or near the habitable zone) and L 98-59 b, c, d, and its vetting passes
them while rejecting signals the TESS team had already classified as false positives.
See `REPORT.md` for the full table.

## Running it

```bash
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python \
  numpy scipy astropy astroquery pandas matplotlib wotan transitleastsquares numba \
  pytest requests tqdm pyarrow lightkurve setuptools
# MAST per-sector lists (light; ~350 MB):
mkdir -p data/scripts && seq 1 107 | xargs -P 6 -I{} curl -s -f -o data/scripts/tesscurl_sector_{}_lc.sh \
  https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{}_lc.sh
.venv/bin/python scripts/01_select_targets.py
.venv/bin/python scripts/03_catalogs.py   # TOI/CTOI/confirmed/TCE/EB catalogues
# heavy steps under the thermal guard (see tools/thermal/):
caffeinate -i -s tools/thermal/daemon.sh &
.venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/02_download.py
.venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/04_search.py --workers 3
.venv/bin/python -m pytest tests/
```

The full run downloads ~13 GB (compressed) and takes several hours on a MacBook Air.

### Thermal guard
This was developed on a fanless MacBook Air. `tools/thermal/daemon.sh` reads macOS's
own thermal-pressure level and the battery temperature every 20 s, freezes heavy jobs
(SIGSTOP) at "serious" pressure or 40 C, resumes them when cool, and terminates them
at "critical" or 45 C. Jobs started through `guarded_run.py` run at low priority in
their own process group so the guard can control them.

## Honest limitations
* A "candidate" here passed automated tests on TESS data alone. Confirming a planet
  needs follow-up (ground-based photometry to rule out nearby eclipsing binaries,
  high-resolution imaging, radial velocities or statistical validation).
* Only 2-minute-cadence stars with >= 20 sectors were searched, and only periods
  0.4-40 days with at least 3 transits.
* Vetting thresholds were tuned on a small validation set (11 stars with known planets)
  and checked with injection-recovery; they are not the TESS team's thresholds.
