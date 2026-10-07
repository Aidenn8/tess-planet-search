# 8. Limitations and reproducibility

## 8.1 Limitations

**Scope of the search.**

* Only stars with SPOC 2-minute photometry in at least 20 sectors, T ≤ 13.5 and TIC contamination ratio
  ≤ 0.2 were searched. Fainter, more crowded and less-observed M dwarfs, and full-frame-image photometry,
  are outside the sample.
* Periods of 0.4–40 days with at least three observed transits. Sensitivity falls steeply below 1 R⊕
  (20% of injected planets found and kept) and towards long periods (Chapter 3).
* Catalogue crossmatches reflect the TOI, CTOI, confirmed-planet, eclipsing-binary and SPOC TCE lists as
  downloaded on 2026 October 4. A signal labelled new may have been reported since.

**Vetting.**

* Thresholds were tuned on a small validation set (11 stars with known planets) and checked with
  injection–recovery and inverted light curves; they are not the TESS team's thresholds.
* Some rules were refined while the search ran; all stored results were re-classified with the final
  rules so that every star is judged alike, but the rules were not frozen before the search began.
* The light-curve centroid test is weak for dips of ~0.1%; the pixel-level localization replaces it for
  the signals that matter.

**Pixel-level localization.**

* Sources about 5″ apart cannot be separated (TIC 198412174's neighbour 5.3″ away).
* The model assumes that only one source in the stamp varies in step with the transit ephemeris. A
  second variable star with a commensurate period breaks this; such fits are flagged by a reduced
  χ² ≥ 2 (Chapter 4), but a weaker version of the same effect would not be flagged.
* Position errors include a 1.5″ systematic floor calibrated on 46 known cases; the calibration sample
  is dominated by bright planet hosts and by synthetic eclipses in the candidates' own fields.
* Gaia DR3 is complete only to G ≈ 21; fainter blended sources are not in the star list, although the
  position-grid fit does not need them to localize the light loss.

**Transit parameters.**

* Stellar parameters come from TIC v8.2 (Mann et al. relations for M dwarfs). The 5% mass uncertainty
  is assumed rather than catalogued. Radii carry the TIC radius error; systematic errors in the stellar
  radius scale would shift all planet radii together.
* Detrending choices move radii: by 1–4% for the candidates' short transits, and by up to ~18% for the
  3.2-hour transit of TOI-700 d (Chapter 5).
* Orbits are assumed circular; limb darkening uses broad generic M-dwarf priors; transit-timing
  variations are not modelled.
* The shape-only stellar-density check is weak, as the L 98-59 c test shows.

**Statistical validation.**

* Only TOI-218 has public high-resolution imaging; for the other three, unresolved companions are
  constrained only by Gaia and the reported FPPs are upper-end values. TOI-218 falls below the numerical
  validation thresholds once its 2020 speckle contrast curve is included, but is not called validated
  (Chapter 6).
* The field-star population was queried from VizieR rather than the ESA Gaia archive, using TRICERATOPS's
  own conversion from Gaia photometry to stellar properties.
* Treating pixel-excluded neighbours as cleared is analogous to clearing by ground-based photometry, but
  it is this work's own method; the TESS-only FPPs are reported alongside.

**Status.** The four candidates are signals that survived every test available from TESS data. They are
not confirmed planets. Confirmation needs ground-based photometry, high-resolution imaging, and
radial-velocity measurements or a full statistical validation.

## 8.2 Reproducing the results

### Environment

Python 3.12 with two virtual environments managed by [`uv`](https://docs.astral.sh/uv/): the main one
for everything except TRICERATOPS, which pins NumPy 1.x and runs separately.

```bash
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python -r requirements.txt
uv venv --python 3.12 .venv-tri && uv pip install --python .venv-tri/bin/python -r requirements-triceratops.txt
```

Principal packages: NumPy 2.5, SciPy 1.18, Astropy 8.0, Astroquery 0.4.11, pandas 3.0, wotan 1.10,
batman 2.5.3, emcee 3.1.6, TESS_PRF 0.1.3, lightkurve 2.5.1, Matplotlib 3.11; TRICERATOPS 1.1.0 (with
NumPy 1.26).

### Pipeline

```bash
# 1. sample and catalogues (MAST download lists, ~350 MB)
mkdir -p data/scripts && seq 1 107 | xargs -P 6 -I{} curl -s -f -o data/scripts/tesscurl_sector_{}_lc.sh \
  https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{}_lc.sh
.venv/bin/python scripts/01_select_targets.py
.venv/bin/python scripts/03_catalogs.py

# heavy steps run under the thermal guard on a fanless laptop (see tools/thermal/)
caffeinate -i -s tools/thermal/daemon.sh &
G=".venv/bin/python tools/thermal/guarded_run.py --"

# 2. light curves (~12 GB), search, reliability
$G .venv/bin/python scripts/02_download.py
$G .venv/bin/python scripts/04_search.py --workers 3
$G .venv/bin/python scripts/05_reliability.py inject --n 300
$G .venv/bin/python scripts/05_reliability.py invert --n 200
.venv/bin/python scripts/06_summarize.py && .venv/bin/python scripts/check_known.py
.venv/bin/python scripts/08_followup.py

# 3. pixel files (~44 GB), localization, NASA cross-check
.venv/bin/python scripts/10_download_tpfs.py                       # candidates and test stars
.venv/bin/python scripts/10_download_tpfs.py <TIC IDs of the weak signals>
.venv/bin/python scripts/11_spoc_period_check.py
$G .venv/bin/python scripts/12_localize.py --workers 3
$G .venv/bin/python scripts/12_localize.py --weak --workers 3
.venv/bin/python scripts/16_weak_recheck.py

# 4. statistical validation and transit fits
.venv/bin/python scripts/13_triceratops_inputs.py
$G .venv-tri/bin/python scripts/13_triceratops_run.py
TRI_RUNS=3 $G .venv-tri/bin/python scripts/13_triceratops_run.py --cleared
$G .venv/bin/python scripts/14_mcmc.py
$G .venv/bin/python scripts/14_mcmc.py --only L98-59c TOI-700d --scale 3

# 5. outputs
.venv/bin/python scripts/15_dossiers.py && .venv/bin/python scripts/09_report.py
.venv/bin/python scripts/07_figures.py && .venv/bin/python docs/make_figures.py
.venv/bin/python docs/make_figures.py --dark && .venv/bin/python site/build.py   # website figures and pages
.venv/bin/python paper/make_figure.py && .venv/bin/python paper/make_note.py
.venv/bin/python -m pytest tests/
```

### Compute

On an Apple M2 MacBook Air (8 cores, 16 GB, no fan, thermal guard active): search about 4.5 hours;
injection and inversion tests a few hours alongside; localization of all candidates, test cases and
injections about 10 minutes; TRICERATOPS about 3 hours (dominated by 10⁶-draw runs); MCMC under an hour.
Raw data: 12 GB of compressed light curves and 44 GB of target pixel files, neither stored in this
repository; every file can be re-downloaded from MAST by the scripts above.

### Determinism

Random elements use fixed seeds: injection draws (seed 42), null-event placement in difference imaging
(per-sector seeds), MCMC initialisation (seeds 1 and 2). TRICERATOPS's Monte Carlo is not seeded; its
run-to-run scatter is reported with every FPP. Results can differ slightly if MAST reprocesses data or if
the catalogues are re-downloaded at a later date.

### Tests

`tests/` contains 13 tests: the search recovers an injected planet and both planets of a two-planet
system in synthetic data and finds nothing strong in pure noise; flares are removed while transits
survive; the period grid, transit masking and catalogue period relations behave as specified; vetting
passes a clean planet, rejects an eclipsing binary with unequal eclipses and flags a secondary eclipse;
the timescale-noise estimate matches white noise when it should; and 19-digit Gaia source IDs survive
the neighbour bookkeeping exactly.

## 8.3 Repository map

| path | contents |
|---|---|
| `tess_search/` | library: target selection, download, cleaning, search, vetting, crossmatch, injection, pixels, localization, SPOC DV parsing, MCMC, assessments |
| `scripts/` | numbered pipeline steps 01–16 |
| `results/` | every result: per-star search records, signal tables, reliability runs, follow-up and hardening outputs, candidate dossiers |
| `results/candidates/` | one dossier per signal and a summary table |
| `results/hardening/` | localization, MCMC, TRICERATOPS, SPOC cross-check and weak-signal re-check outputs; draft ExoFOP upload |
| `docs/` | this technical write-up and its figures (`docs/make_figures.py`) |
| `paper/` | draft Research Note of the AAS and its figure |
| `REPORT.md` | generated summary of every number (`scripts/09_report.py`) |
| `EXPLAINER.md` | plain-language summary |
| `tools/thermal/` | thermal guard for long jobs on a fanless laptop |
| `tests/` | unit tests |

## 8.4 Data and acknowledgements

This work uses data from the TESS mission, funded by NASA's Science Mission Directorate, obtained from
the Mikulski Archive for Space Telescopes (MAST) at the Space Telescope Science Institute; light curves,
target pixel files and Data Validation products from the Science Processing Operations Center (SPOC) at
NASA Ames Research Center; the TESS Input Catalog; data from the European Space Agency mission Gaia,
processed by the Gaia Data Processing and Analysis Consortium, accessed through the VizieR service at CDS
Strasbourg; the Exoplanet Follow-up Observing Program (ExoFOP) and the NASA Exoplanet Archive, both
operated by the California Institute of Technology under contract with NASA; and the open-source
packages Astropy, Astroquery, NumPy, SciPy, pandas, Matplotlib, wotan, batman, emcee, lightkurve,
TESS_PRF and TRICERATOPS.
