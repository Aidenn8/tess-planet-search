# Follow-up guide

This page is for observers and collaborators. For each candidate it states what is already established,
what is still open, which observation would settle it, and when that observation can be made. Every
number comes from the result files in this repository; the status is as of October 2026.

## At a glance

| candidate | P (d) | depth (ppm) | T | open question → next observation |
|---|---|---|---|---|
| **TOI-218**, new signal | 2.1468 | 1273 | 13.3 | no ground-based detection yet → multi-band photometry of a predicted transit resolving the 13.5″ binary |
| **TIC 229689348** | 0.4654 | 693 | 13.0 | unresolved bound companion → high-resolution imaging |
| **TIC 149390648** | 2.8369 | 798 | 13.4 | crowded field, unresolved companion → seeing-limited photometry of the field; imaging |
| **TIC 198412174** | 1.3902 | 382 | 11.8 | a T = 18.2 star 5.3″ away is not excluded → photometry resolving it (it would need an ~11% eclipse); imaging |

The pixel-level localization puts all four signals on their target stars; for TIC 198412174 one close
neighbour remains possible.

TRICERATOPS false-positive probabilities, with the neighbours that the pixel localization excludes treated
as cleared: TOI-218 0.088 (**1.4 × 10⁻⁵** with the existing speckle imaging), TIC 229689348 0.197,
TIC 149390648 0.088, TIC 198412174 0.478. All four meet the TRICERATOPS *likely planet* criteria
(FPP < 0.5, NFPP < 0.001). The evidence for each is in [Chapter 6](docs/06-candidates.md) and the
per-candidate dossiers in [`results/candidates/`](results/candidates/).

## Ephemerides

| TIC | RA, Dec (°) | T₀ (BJD_TDB) | P (d) | T₁₄ (h) |
|---|---|---|---|---|
| 32090583 | 58.42035, −68.73945 | 2459318.6971 ± 0.0010 | 2.1467982 ± 0.0000019 | 0.93 |
| 229689348 | 263.23532, +57.09234 | 2459027.6289 ± 0.0005 | 0.4653769 ± 0.0000005 | 0.53 |
| 149390648 | 84.07329, −62.20426 | 2459180.6342 ± 0.0015 | 2.8368867 ± 0.0000038 | 1.14 |
| 198412174 | 259.21165, +59.82059 | 2459644.3198 ± 0.0008 | 1.3901602 ± 0.0000014 | 0.61 |

Coordinates are TIC positions (epoch J2000). Periods and epochs come from least-squares transit fits to
all TESS sectors of each star (period errors inflated for correlated noise; T₀ near the middle of each
baseline); durations and depths (1273 ± 106, 693 ± 65, 798 ± 88 and 382 ± 43 ppm) are MCMC medians.
Mid-transit uncertainties propagated to 2027 are 3.2, 3.7, 5.0 and 3.0 minutes, so no ephemeris
refinement is needed before observing. Machine-readable:
[`results/candidates/summary.csv`](results/candidates/summary.csv).

## Upcoming transits

[`scripts/17_transit_predictions.py`](scripts/17_transit_predictions.py) lists every transit between
2026-10-06 and 2027-03-31 for which the whole transit plus 30 minutes of baseline on each side occurs with
the target above 30° and the Sun more than 12° below the horizon, at the six Las Cumbres Observatory (LCO)
sites:

| candidate | transits | CTIO | SAAO | SSO | Teide | McDonald | Haleakala |
|---|---|---|---|---|---|---|---|
| TOI-218 | 82 | 16 | 17 | 14 | – | – | – |
| TIC 229689348 | 378 | – | – | – | 18 | 12 | 6 |
| TIC 149390648 | 62 | 11 | 12 | 11 | – | – | – |
| TIC 198412174 | 126 | – | – | – | 6 | 6 | 3 |

The full list (132 event–site pairs, UTC windows and minimum altitudes) is in
[`results/followup_planning/transits.csv`](results/followup_planning/transits.csv), and the
[project page](https://aidenn8.github.io/tess-planet-search/follow-up.html#next-transits) lists the next
eight per candidate. These are planning estimates: confirm with a scheduling tool such as the Swarthmore
Transit Finder (TAPIR), which accepts a user-supplied ephemeris.

<!-- transit-tables -->

## What each candidate needs

### TOI-218: a third signal

**Established.** The signal is on TIC 32090583 in the pixels (2.2 ± 3.0″) and the near-twin wide-binary
companion, Gaia DR3 4667466549703138304 at 13.5″, is excluded at 4.4σ; the two known TOIs localize to the
same star. Gemini-South 'Zorro' speckle imaging from 2020 (on ExoFOP, PI S. Howell) reaches Δ = 4.4 mag
at 562 nm and 5.5 mag at 832 nm at 0.5″. With that contrast curve and the pixel-excluded neighbours cleared,
TRICERATOPS gives FPP = 1.4 × 10⁻⁵ and NFPP < 10⁻⁵, below the thresholds Giacalone et al. (2021) use
for statistical validation (FPP < 0.015, NFPP < 0.001).

**Open.** No ground-based light curve has seen the transit yet, and TFOP notes describe the host as a
flaring (eruptive) variable. Of the 22 TFOP time-series observations on ExoFOP, the public photometry
(27 AstroImageJ tables and one TRAPPIST light curve, all from observations up to September 2025) contains
one light curve that covers more than half of a predicted transit: LCO-CTIO 1 m on 2018-12-03 (87% coverage, timing ±1.8 min). It starts at ingress with no
pre-transit baseline and is inconclusive: the target-to-companion flux ratio dips by 1874 ± 608 ppm with
a flat baseline and 586 ± 1221 ppm with a linear one, against 1273 ppm expected
([`scripts/18_archival_ground_check.py`](scripts/18_archival_ground_check.py),
[`results/hardening/archival/`](results/hardening/archival/)). TFOP's four observations of TOI-218 since
October 2025 have no public photometry yet; checking them against this ephemeris costs nothing.

**Observation that settles it.** Photometry of a predicted transit with at least an hour of baseline on
both sides, with apertures that separate the 13.5″ pair. At 1.3 ppt in a 56-minute transit this is
within reach of a 2-m-class telescope; simultaneous multi-band photometry (e.g. LCO's MuSCAT4 at Siding
Spring, which has observed TOI-218 for TFOP since 2024) also tests for colour-dependent depths. The
companion would need a 0.13% eclipse to produce the signal, so the same data clear it directly.

### TIC 229689348: an 11.2-hour orbit

**Established.** On target in the pixels at the correct period with all 25 sectors (3.6 ± 3.3″); the
55″ offset in NASA's DV report came from difference images that failed SPOC's own quality test. The
nearest neighbour (9.7″, T = 15.9) is excluded at 4.6σ and the bright star SPOC's offset pointed towards
(48.3″) at 8.0σ. FPP = 0.197, almost all of it a planet on an unresolved bound companion.

**Observation that settles it.** High-resolution imaging (AO or speckle) to detect or exclude a close
companion. A seeing-limited check of the neighbours is straightforward: the 9.7″ star would need a
0.73% eclipse. The transit itself (0.69 ppt, 32 min) is below the ~1 ppt general threshold of TFOP's seeing-limited
photometry sub-group (SG1) for detecting a transit on the target
([Collins 2019](https://tsc.mit.edu/2019/talks/Karen%20Collins.pdf)). Northern target (Dec +57°).

### TIC 149390648: a crowded field

**Established.** On target (4.2 ± 3.6″), in agreement with SPOC's own difference-image offset; all 89
Gaia neighbours in the ~100″ field are excluded at ≥ 3.5σ. FPP = 0.088, NFPP = 0.00093 with the
pixel-excluded neighbours cleared. The only neighbour scenario left is TIC 149390646, a TIC entry with
T = 17.0 at a position where Gaia DR3 lists only a G = 21.6 source.

**Open.** The pixels lose 1.47 ± 0.23 times the light the light-curve depth predicts, the highest ratio
in the set, which in a crowded field can mean the SPOC crowding correction is imperfect. The first half of the data alone does not lock onto the period.

**Observation that settles it.** Seeing-limited photometry of a predicted transit covering the field:
the least-excluded neighbour (14.7″, T = 17.4) would need a 2.8% eclipse. High-resolution imaging for the
bound-companion scenarios. Southern target (Dec −62°).

### TIC 198412174: which star?

**Established.** On or near the target (2.2 ± 3.5″); every neighbour beyond 18″ is excluded at ≥ 5.4σ.

**Open.** A T = 18.2 star 5.3″ away at position angle 115° is excluded at only 1.3σ: two sources this close
are below the localization's resolution. TRICERATOPS (FPP 0.478) weighs a planet on an unresolved bound
companion (48%) about equally with the target (43%), because the transit is short for this star.

**Observation that settles it.** Two independent observations. (1) Seeing-limited photometry that
resolves the 5.3″ pair during a predicted transit: that star would need an ~11% eclipse, easily seen
with a 1-m telescope in good seeing; a flat light curve on it clears it. (2) AO or speckle imaging to test
for a bound companion. The 382-ppm transit on the target is too shallow for most ground-based
photometry. Northern target (Dec +60°).

## Neighbours that could produce each signal

For each candidate, the Gaia DR3 stars in the localization field that are bright enough to produce the
signal (required eclipse depth below 100%), least-excluded first. Separations and position angles (east
of north) are from the target at the median TESS epoch; "pixels" is how strongly the TESS difference
images exclude the star as the source, including a 1.5″ systematic floor. Full list:
[`results/followup_planning/neighbours.csv`](results/followup_planning/neighbours.csv)
([`scripts/19_neighbour_checklist.py`](scripts/19_neighbour_checklist.py)). For signals shallower than
about 1 ppt, SG1's aim is to check every Gaia star within 2.5′ of the target for a nearby eclipsing binary
(Collins 2019); the localization field covers about 100″, so stars between 100″ and 150″ still need that
check.

**TOI-218** (11 neighbours in the field; 5 could produce the signal; the first is the wide-binary companion)

| Gaia DR3 | sep | PA | G | T | eclipse needed | pixels |
|---|---|---|---|---|---|---|
| 4667466549703138304 | 13.5″ | 226° | 14.8 | 13.6 | 0.13% | 4.4σ |
| 4667466652782354176 | 36.1″ | 65° | 17.4 | 17.0 | 0.53% | 8.2σ |
| 4667466476687967488 | 37.0″ | 158° | 20.3 | 19.9 | 5.6% | 8.4σ |

**TIC 229689348** (20 neighbours; 10 could produce the signal)

| Gaia DR3 | sep | PA | G | T | eclipse needed | pixels |
|---|---|---|---|---|---|---|
| 1422141810046979968 | 9.7″ | 49° | 16.4 | 15.9 | 0.73% | 4.6σ |
| 1422141805751675264 | 22.2″ | 148° | 20.7 | 20.2 | 22% | 6.7σ |
| 1422142462882048000 | 25.3″ | 267° | 18.7 | 17.7 | 1.7% | 7.2σ |

**TIC 149390648** (89 neighbours; 54 could produce the signal)

| Gaia DR3 | sep | PA | G | T | eclipse needed | pixels |
|---|---|---|---|---|---|---|
| 4757656533594500352 | 14.7″ | 311° | 18.0 | 17.4 | 2.8% | 3.5σ |
| 4757656529296176896 | 16.4″ | 357° | 20.3 | 20.1 | 29% | 4.5σ |
| 4757656533592367872 | 16.6″ | 57° | 20.7 | 20.5 | 29% | 6.0σ |

**TIC 198412174** (16 neighbours; 9 could produce the signal)

| Gaia DR3 | sep | PA | G | T | eclipse needed | pixels |
|---|---|---|---|---|---|---|
| 1437353828295869824 | 5.3″ | 115° | 18.6 | 18.1 | 11% | **1.3σ** |
| 1437353858359596928 | 18.3″ | 12° | 20.3 | 19.6 | 17% | 5.4σ |
| 1437353789640128000 | 35.2″ | 146° | 20.5 | 19.6 | 3.2% | 6.3σ |

## Reporting the candidates

**ExoFOP community TOIs (CTOIs).** ExoFOP paused community-candidate uploads on 2026-03-31 and reopened
them on 2026-08-19 with two new conditions: candidates must first be published in a peer-reviewed
journal with online access, and uploaders must be approved through ExoFOP's *Published Candidate Upload
Request* form. A Research Note of the AAS qualifies only if the methodology used to detect and vet the
candidate was previously published in a peer-reviewed journal and is cited by the note
([ExoFOP news](https://exofop.ipac.caltech.edu/tess/news.php),
[candidate upload help](https://exofop.ipac.caltech.edu/tess/candidate_help.php)). An upload file in
ExoFOP's bulk format is ready in
[`results/hardening/exofop/params_planet_DRAFT.txt`](results/hardening/exofop/params_planet_DRAFT.txt);
it needs the paper URL and the uploader's tag. Candidate numbers follow ExoFOP's convention: TIC 32090583.01
and .02 are already used by TOI-218.01 and .02, so the new signal is TIC 32090583.03. Its period is not a
multiple of either known period (ratios 4.90 and 3.89), so it is a new candidate rather than an
additional parameter set for an existing one.

**TESS Follow-up Observing Program (TFOP).** Seeing-limited photometry (SG1), reconnaissance
spectroscopy (SG2), high-resolution imaging (SG3), precise radial velocities (SG4) and space photometry
(SG5) for TESS candidates are coordinated within TFOP's sub-groups
([TFOP](https://tess.mit.edu/followup/), [how to join](https://tess.mit.edu/followup/apply-join-tfop/)).
TOI-218 is already observed within TFOP for TOI-218.01 and .02 (22 time-series and two imaging entries
on ExoFOP).

## Contact

Questions, corrections and observations are welcome through
[GitHub issues](https://github.com/Aidenn8/tess-planet-search/issues) on this repository.
