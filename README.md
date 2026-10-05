# TESS M-dwarf deep search

A search for transiting planets around the **1,279 bright M dwarfs that NASA's TESS
telescope has watched the longest** (20 to 44 sectors of 2-minute data each, through
sector 107), with automated vetting, comparison against every public catalogue,
and measured completeness and false-alarm rates.

> **Results** (full report: [`REPORT.md`](REPORT.md); plain-language version:
> [`EXPLAINER.md`](EXPLAINER.md)):
>
> * **5 Earth-sized candidates in no planet catalogue** passed every test and a deeper
>   follow-up (physical transit fit, Gaia neighbours, independent half-data searches),
>   including a **third signal in the TOI-218 system** (2.147 d, 1.0 R_earth) and an
>   **11-hour orbit** around TIC 229689348 (1.2 R_earth). Each has a dossier in
>   `results/followup/` with what a community-TOI submission needs.
> * Recovers **25 of 26 confirmed transiting planets** in range, including all four of TOI-700.
> * **Completeness:** 74% of 300 injected planets found and kept (89% for 2-4 R_earth
>   inside 15 days). **Reliability:** no false candidates in 200 flipped light curves.
>
> Candidates are signals worth follow-up observations, not confirmed planets.

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
| 8 | `scripts/08_followup.py` | for each surviving candidate: fresh vetting, physical transit fit, Gaia DR3 neighbours, half-data searches, dossier |
| 9 | `scripts/09_report.py` | writes `REPORT.md` from the result files |
| 10 | `scripts/10_download_tpfs.py` | target pixel files for the candidates, the test stars and the weak signals |
| 11 | `scripts/11_spoc_period_check.py` | SNR of this light curve at each SPOC TCE period (why SPOC's SNR fell) |
| 12 | `scripts/12_localize.py` | pixel-level source localization, with its validation (confirmed planets, TFOP nearby EBs, injected eclipses); `--weak` for the weak signals |
| 13 | `scripts/13_triceratops_inputs.py`, `13_triceratops_run.py` | TRICERATOPS false-positive probabilities (separate environment); `--cleared` treats neighbours excluded by step 12 as cleared |
| 14 | `scripts/14_mcmc.py` | MCMC transit fits (transits masked in the detrending, stellar-density prior, stellar-radius error propagated), checked on TOI-700 d and L 98-59 c |
| 15 | `scripts/15_dossiers.py` | final dossiers (`results/candidates/`) and the draft ExoFOP upload |
| 16 | `scripts/16_weak_recheck.py` | do the weak signals' dips appear in the pixels? |
| paper | `paper/make_figure.py`, `paper/make_note.py` | figure and draft Research Note of the AAS |

Step 3 (`scripts/03_catalogs.py`) downloads the TOI, CTOI, confirmed-planet, SPOC TCE and
eclipsing-binary catalogues used for crossmatching.

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

### Hardening: pixels, statistics and NASA's own pipeline
* **Pixel-level localization** (`tess_search/localize.py`): for every sector, the images taken during
  transit are subtracted from those just before and after (difference images; per-pixel errors from ~60
  fake transits per sector, flares removed). All sectors are fitted together with the SPOC pixel response
  function, calibrated per sector on the Gaia DR3 stars in the stamp (proper motions applied), to find
  where on the sky the light went missing. Errors include a 1.5" systematic floor measured on 46 sources of
  known position. A second test asks whether the target loses, in the pixels, the light the light-curve
  depth predicts (confirmed planets: 0.85-1.29).
* **TRICERATOPS** (Giacalone et al. 2021) false-positive probabilities with the Gaia DR3 field population
  (queried through VizieR), with and without the neighbours the localization excludes.
* **NASA SPOC Data Validation reports** for the candidates SPOC flagged (`tess_search/spoc_dv.py`).
* **MCMC transit fits** (`tess_search/mcmc.py`).

### Crossmatch (`tess_search/crossmatch.py`)
TOIs, community TOIs, confirmed planets (NASA Exoplanet Archive), the TESS eclipsing-binary
catalogue, and **every SPOC Threshold Crossing Event** (132 single- and multi-sector
runs), including signals NASA's pipeline found but never promoted.

## Validation on known planets
Across the searched stars the pipeline independently recovers 25 of the 26 confirmed
transiting planets with periods of 0.4-40 days (the miss is TOI-1752 c at 32.7 d), including
all four TOI-700 planets (TOI-700 d and e are Earth-sized, in or near the habitable zone)
and L 98-59 b, c and d. Its vetting passes them while rejecting signals the TESS team had
already classified as false positives. The MCMC transit fit gives TOI-700 d 1.18 ± 0.05 R_earth
(published 1.07 ± 0.06, Gilbert et al. 2023) and L 98-59 c 1.34 ± 0.04 R_earth (published 1.39 ± 0.09,
Demangeon et al. 2021). See `REPORT.md`.

## Running it

```bash
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python -r requirements.txt
uv venv --python 3.12 .venv-tri && uv pip install --python .venv-tri/bin/python -r requirements-triceratops.txt
# MAST per-sector lists (light; ~350 MB):
mkdir -p data/scripts && seq 1 107 | xargs -P 6 -I{} curl -s -f -o data/scripts/tesscurl_sector_{}_lc.sh \
  https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{}_lc.sh
.venv/bin/python scripts/01_select_targets.py
.venv/bin/python scripts/03_catalogs.py   # TOI/CTOI/confirmed/TCE/EB catalogues
# heavy steps under the thermal guard (see tools/thermal/):
caffeinate -i -s tools/thermal/daemon.sh &
.venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/02_download.py
.venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/04_search.py --workers 3
.venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/05_reliability.py inject --n 300
.venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/05_reliability.py invert --n 200
.venv/bin/python scripts/06_summarize.py && .venv/bin/python scripts/check_known.py
.venv/bin/python scripts/08_followup.py && .venv/bin/python scripts/07_figures.py && .venv/bin/python scripts/09_report.py
# hardening (pixel files ~20 GB for the candidates and test stars, ~28 GB more for the weak signals)
.venv/bin/python scripts/10_download_tpfs.py && .venv/bin/python scripts/11_spoc_period_check.py
.venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/12_localize.py --workers 3
.venv/bin/python scripts/13_triceratops_inputs.py
.venv/bin/python tools/thermal/guarded_run.py -- .venv-tri/bin/python scripts/13_triceratops_run.py
.venv/bin/python tools/thermal/guarded_run.py -- .venv-tri/bin/python scripts/13_triceratops_run.py --cleared
.venv/bin/python tools/thermal/guarded_run.py -- .venv/bin/python scripts/14_mcmc.py
.venv/bin/python scripts/15_dossiers.py && .venv/bin/python scripts/09_report.py
.venv/bin/python paper/make_figure.py && .venv/bin/python paper/make_note.py
.venv/bin/python -m pytest tests/
```

The full run downloads ~12 GB (compressed light curves) and took about 6 hours of
computing on a fanless MacBook Air M2 (search ~4.5 h with 3-4 workers, reliability tests
alongside). The bulk diagnostic sheets (`results/plots/`, ~200 MB) are not in git;
`04_search.py` regenerates them.

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
  and checked with injection-recovery and inverted light curves; they are not the TESS
  team's thresholds. Some rules were refined while the search ran; every stored result
  is re-classified with the final rules (`06_summarize.py`), so all stars are judged alike.
* The light-curve centroid test is weak for dips this shallow; the pixel-level localization
  replaces it for the candidates, but it cannot separate sources closer than ~5" (TIC 198412174's
  4.7" neighbour) and it assumes one variable source per image (a variable star whose period is
  commensurate with the signal breaks this; such fits are flagged by their poor reduced chi2).
* Planet radii use TIC v8.2 stellar parameters (radius error propagated). The detrending choice
  moves radii by 2-4% for the candidates' short transits and by up to ~10% for long ones
  (TOI-700 d: 1.18 ± 0.05 R_earth with transits masked vs 1.07 published).
* TRICERATOPS was run without high-resolution imaging, so its FPPs are upper-end values; no
  candidate is statistically validated.
