# 6. Candidates

Five signals that are in no planet catalogue passed every light-curve test. After the pixel-level,
statistical and pipeline checks of Chapters 4 and 5, four remain planet candidates and one is a nearby
eclipsing binary. This chapter assembles the evidence for each. Full per-signal dossiers, with every
number and figure, are in [`results/candidates/`](../results/candidates/); a draft ExoFOP community-TOI
upload for the four candidates is in [`results/hardening/exofop/`](../results/hardening/exofop/).

A candidate is a signal that has survived every test available from TESS data and is worth follow-up
observations. None of these is a confirmed planet.

| | TOI-218 (new signal) | TIC 229689348 | TIC 149390648 | TIC 198412174 | TIC 294053492 |
|---|---|---|---|---|---|
| verdict | candidate (strongest) | candidate | candidate | candidate (ambiguous host) | nearby eclipsing binary |
| period (d) | 2.1467982 | 0.4653769 | 2.8368867 | 1.3901602 | 1.0852647 |
| radius (R⊕) | 1.05 ± 0.06 | 1.25 ± 0.08 | 1.00 ± 0.07 | 1.35 ± 0.11 | – |
| T_eq (K) | 578 | 1131 | 581 | 954 | – |
| SNR | 11.1 | 10.4 | 9.1 | 9.7 | 7.9 |
| source offset | 2.2 ± 3.0″ | 3.6 ± 3.3″ | 4.2 ± 3.6″ | 2.2 ± 3.5″ | 21.6 ± 3.1″ |
| FPP (cleared) | 0.088 | 0.197 | 0.088 | 0.478 | – |
| prior status | in no list | SPOC TCE, never a TOI | SPOC TCE, never a TOI | SPOC TCE, never a TOI | in no list |

---

## 6.1 TOI-218: a third transiting signal

**Star.** TIC 32090583 = Gaia DR3 4667466549703138176; RA 58.42035°, Dec −68.73945°; T = 13.33;
T_eff = 3249 ± 157 K, R⋆ = 0.285 ± 0.009 R☉, M⋆ = 0.259 M☉ (TIC v8.2); 52 pc; Gaia RUWE 1.12; 41 sectors
(Sectors 1–98). TOI-218 already has two TOIs: TOI-218.01 (0.438 d) and TOI-218.02 (8.352 d).

**Signal.** 2.1467982 ± 0.0000019 d, depth 1273 ± 106 ppm, T₁₄ = 0.93 h, b = 0.51; R_p = 1.05 ± 0.06 R⊕,
T_eq ≈ 578 K. Found by the iterative search after both known TOIs were masked. Its period is not a
simple ratio of either known period (4.90 and 3.89). It matches no TOI, CTOI, confirmed planet,
eclipsing binary or SPOC TCE.

**Evidence.**

* 437 transits with data; red-noise-aware SNR 11.1; recovered independently in both halves of the data
  (SNR 11.0 and 8.4, same period).
* Odd and even transits agree (0.9σ); no secondary eclipse; no single transit carries more than 3% of the
  signal; flat-bottomed shape consistent with a planet.
* **A wide binary.** TOI-218 has a near-twin companion 13.5″ away: Gaia DR3 4667466549703138304, with the
  same parallax (19.03 mas) and proper motion (139.9, 191.2 mas yr⁻¹) and only 0.17 mag fainter in the
  TESS band. The two are blended in every TESS aperture. The pixel localization places the new signal on
  TIC 32090583 (2.2 ± 3.0″) and excludes the companion at 4.4σ; the two known TOIs also come from
  TIC 32090583 (companion excluded at 8.7σ and 6.8σ). Before the Gaia-ID fix described in the
  repository history, this companion was missing from the candidate's neighbour table.
* The target loses 1.25 ± 0.18 times the light the light-curve depth predicts, within the range of
  confirmed planets.
* TRICERATOPS: FPP = 0.41 from photometry alone, almost entirely the scenario that the planet orbits the
  companion; with the companion cleared by the pixel localization, **FPP = 0.088, NFPP < 10⁻⁵**; adding the
  existing speckle imaging, **FPP = 1.4 × 10⁻⁵, NFPP < 10⁻⁵** (below).
* Additional candidates in systems that already host candidates are statistically more likely to be
  real (Lissauer et al. 2012).

**Open questions.** TFOP notes describe the host as an eruptive (flaring) variable and suspect that
TOI-218.01 may be stellar variability. Flares are removed before the search, and the new signal is a
periodic, flat-bottomed, 0.9-hour dip present in hundreds of transits and in both halves of the data,
which does not resemble spot modulation or flaring. No ground-based light curve has yet seen the transit.
ExoFOP lists 22 TFOP time-series observations of TOI-218, obtained for the two known TOIs, with public
photometry for those taken before October 2025: 27 AstroImageJ tables, a TRAPPIST light curve and
joint-fit subsets, all checked against this ephemeris (`scripts/18_archival_ground_check.py`). Only one,
LCO-CTIO 1 m on 2018-12-03 (tag 1415), covers more than half of a predicted transit (87%, timing
uncertainty ±1.8 min). It starts at ingress, with no pre-transit baseline, and the answer depends on the
baseline model: the target-to-companion flux ratio dips by 1874 ± 608 ppm with a flat baseline and by
586 ± 1221 ppm with a linear one, against 1273 ppm expected. It neither confirms nor rules out the signal.
The four TFOP observations since October 2025 have no public photometry yet.

**High-resolution imaging.** ExoFOP also lists Gemini-South 'Zorro' speckle imaging of TOI-218
(2020-11-27; PI S. Howell), with contrast limits of Δ = 4.4 mag (562 nm) and 5.5 mag (832 nm) at 0.5″;
ExoFOP lists no detected companion. Adding the public 562-nm contrast curve to the cleared TRICERATOPS run gives
**FPP = 1.4 × 10⁻⁵ and NFPP < 10⁻⁵**, far below the thresholds Giacalone et al. (2021) use for
validation (FPP < 0.015, NFPP < 0.001). This work does not call the planet validated: the host is a
flaring star, the 2020 imaging was taken for the two known TOIs, and a ground-based detection of the
transit is the standard independent check.

**What would settle it.** Ground-based photometry of a predicted transit that resolves the 13.5″ pair,
with baseline on both sides: at 1.3 ppt the transit is within reach of a 2-m telescope such as LCO's
MuSCAT4, which has observed TOI-218 for TFOP since 2024 (`FOLLOW_UP.md` lists observable windows). The
TFOP team's own unpublished light curves of TOI-218 may already cover transits of the new signal.

---

## 6.2 TIC 229689348: an 11.2-hour orbit

**Star.** Gaia DR3 1422141810046979584; RA 263.23532°, Dec +57.09234°; T = 12.99; T_eff = 3365 ± 157 K,
R⋆ = 0.430 ± 0.013 R☉, M⋆ = 0.425 M☉; 77 pc; RUWE 1.02; 25 sectors (Sectors 15–86).

**Signal.** 0.4653769 ± 0.0000005 d (11.17 hours), depth 693 ± 65 ppm, T₁₄ = 0.53 h, b = 0.80;
R_p = 1.25 ± 0.08 R⊕; a = 0.0089 AU; T_eq ≈ 1130 K; 272 times Earth's insolation. An ultra-short-period
(USP) candidate.

**Evidence.**

* 1,088 transits with data; SNR 10.4; both halves recover it (8.5 and 7.8); odd/even 1.0σ; no secondary;
  depths in SAP and PDCSAP flux agree (670 and 658 ppm); the folded transit is flat-bottomed.
* **NASA's pipeline saw it.** SPOC flagged it in three multi-sector runs (Sectors 14–50, 14–55, 14–86) but
  it was never promoted to a TOI. SPOC's reported SNR (5.8, 18.4, 3.8) tracks the accuracy of SPOC's period
  estimate in each run, not the signal (Chapter 5). In the best run SPOC's difference-image centroid put the
  source 55 ± 10″ away, which may be why it was not promoted; none of that run's 22 per-sector difference
  images passed SPOC's quality metric.
* **The pixels put it on the target.** The joint localization at the correct period with all 25 sectors
  (1,138 transit events) places the source 3.6 ± 3.3″ from the target and excludes the neighbours at 9.7″
  (4.6σ) and 48.3″ (the bright star SPOC's offset points towards, 8.0σ). Synthetic eclipses planted on
  those neighbours are correctly traced to them, so the method could have seen an off-target source here.
  The pixel/light-curve depth ratio is 1.02 ± 0.15.
* TRICERATOPS: FPP = 0.199 with NFPP below 10⁻⁵ once the six pixel-excluded neighbours are cleared. The
  planet-on-target scenario carries 68%; the remaining probability is almost all unresolved bound
  companions (20% for a planet transiting a companion).

**Open questions.** The transit is short for this star, which the fit explains with b ≈ 0.80; the
shape-only density leans towards a denser host (3.1 +1.4/−2.0 × the TIC value), which is what gives the
companion scenarios their weight.

**What would settle it.** High-resolution imaging to exclude a close bound companion; ground-based
photometry of a transit (a 0.07% dip on a T = 13 star is at the edge of what 1-m class telescopes
achieve, but ruling out deep eclipses on the neighbours is straightforward). At about 1.3 R⊕ and 11 hours,
it would join the small population of Earth-sized USP planets around M dwarfs.

---

## 6.3 TIC 149390648

**Star.** Gaia DR3 4757656533594499840; RA 84.07329°, Dec −62.20426°; T = 13.36; T_eff = 3382 ± 157 K,
R⋆ = 0.345 ± 0.010 R☉, M⋆ = 0.328 M☉; 74 pc; RUWE 1.23; 32 sectors (Sectors 2–96).

**Signal.** 2.8368867 ± 0.0000038 d, depth 798 ± 88 ppm, T₁₄ = 1.14 h, b = 0.51; R_p = 1.00 ± 0.07 R⊕;
T_eq ≈ 581 K.

**Evidence.**

* 242 transits; SNR 9.1; odd/even 0.1σ; no secondary; SAP and PDCSAP depths agree (818 and 829 ppm).
* SPOC flagged it once (Sectors 1–96, SNR 8.5) and never promoted it. SPOC's difference-image offset was
  2.2 ± 3.3″, consistent with the target, and this work's localization agrees (4.2 ± 3.6″), with every
  neighbour excluded at more than 3σ even though the field is crowded (28 Gaia neighbours within 63″
  bright enough to mimic the dip if they were eclipsing binaries).
* TRICERATOPS: FPP = 0.088 and NFPP = 0.00093 once 15 pixel-excluded neighbours are cleared. The last
  neighbour scenario left is a TIC entry with no Gaia counterpart whose TIC brightness is inconsistent with
  the G = 21.6 source Gaia sees there (Chapter 5).

**Open questions.** The first half of the data alone does not lock onto the period (SNR 5.4, at a
different period), while the second half does (7.6). The target loses 1.47 ± 0.23 times the expected
light, the highest ratio among the candidates; in a crowded field this can indicate an imperfect crowding
correction in the PDCSAP light curve.

**What would settle it.** Ground-based photometry that resolves the crowded field during a transit.

---

## 6.4 TIC 198412174

**Star.** Gaia DR3 1437353828295347712; RA 259.21165°, Dec +59.82059°; T = 11.82; T_eff = 3713 ± 157 K,
R⋆ = 0.573 ± 0.017 R☉, M⋆ = 0.565 M☉; 77 pc; RUWE 1.21; 37 sectors (Sectors 14–86). The brightest and
largest host in the set.

**Signal.** 1.3901602 ± 0.0000014 d, depth 382 ± 43 ppm, T₁₄ = 0.61 h, b = 0.93 (near-grazing);
R_p = 1.35 ± 0.11 R⊕; T_eq ≈ 954 K.

**Evidence.**

* 557 transits; SNR 9.7; both halves recover it (7.6 and 7.4); odd/even 0.0σ; no secondary.
* SPOC flagged it in three runs (Sectors 14–50, 14–78, 14–86) with falling SNR (8.5, 6.0, 4.4); as for
  TIC 229689348, the fall follows SPOC's period error, and at this work's period the light curve gives
  SNR 10.2.
* The pixels place the source on or very near the target (2.2 ± 3.5″); every other neighbour (the next
  nearest is 18″ away) is excluded.

**Open questions.** Two separate ambiguities.

* A T = 18.2 star 5.3″ from the target is excluded at only 1.3σ; it would need an ~11% eclipse. Two
  sources this close are below the localization's resolution (Chapter 4).
* TRICERATOPS gives FPP = 0.478, dominated by a planet transiting an unresolved bound companion (48%)
  rather than the target (43%). The transit is short for this star's size, which a smaller, denser host
  would explain; the shape-only density (6.8 +9.4/−6.6 × TIC) leans the same way but is unconstrained.

**What would settle it.** High-resolution imaging (adaptive optics or speckle) to detect or exclude both
a close bound companion and the 5.3″ star's contribution. Seeing-limited photometry that resolves the
5.3″ pair can also test that star directly: it would need an ~11% eclipse, easy to see in a 1-m telescope. If a companion hosts the transit, the planet
would be larger than 1.35 R⊕.

---

## 6.5 TIC 294053492: a nearby eclipsing binary

**Star.** Gaia DR3 5486535779627128448; RA 107.26759°, Dec −58.43336°; T = 13.43; 97 pc; 21 sectors.

**Signal.** 1.0852647 d, depth ~860 ppm in the light curve.

This signal passed every light-curve test and was listed as a candidate in the first version of this
work. Its light-curve cautions were mild: a somewhat V-shaped dip, a lower depth in SAP than in PDCSAP
flux (833 vs. 1,107 ppm in the vetting measurement), a weaker second half of the data (SNR 4.8 vs. 8.3),
and a crowded field with 17 Gaia neighbours bright enough to mimic the dip.

**The pixels.** The light loss lies **21.6 ± 3.1″ north-east of the target**, which is excluded as the
source at 7.7σ. Localizing from the odd and even sectors separately gives the same off-target position
(target excluded at 6.1σ and 6.8σ). More light goes missing than the target could lose for the observed
depth (1.76 ×). The best-fit position lies 2.1″ from Gaia DR3 5486535775331400832, a G = 19.8 background
star (parallax 1.2 mas) that would need a ~16% eclipse to produce the dip: an ordinary eclipsing binary.
That star is fainter than the G < 19.5 limit of the first localization run, which is why the limit was
extended to G = 21.

**The statistics.** TRICERATOPS, without pixel information, also gives substantial weight to a neighbour
(NFPP = 0.30) but attributes it to the bright star 18.4″ to the north-west, which the pixels exclude at
8.3σ. The faint true source carries little prior probability in TRICERATOPS. The two approaches agree that
the target is not the host and disagree about which neighbour is, and the pixel data are the direct
measurement.

This signal is reported as a false positive and is not included in the community-TOI draft.

![Localization maps](figures/fig06_localization_maps.png)

*Figure 6 (repeated from Chapter 4). Localization of the four candidates and of TIC 294053492.*

## 6.6 Follow-up priorities

| priority | candidate | most useful observation | why |
|---|---|---|---|
| 1 | TOI-218 | photometry of a predicted transit resolving the 13.5″ pair; TFOP's unpublished light curves checked against the new ephemeris | strongest signal; meets TRICERATOPS's validation thresholds with the existing speckle imaging; third member of a known system |
| 2 | TIC 229689348 | high-resolution imaging; ground-based check of the 10″ and 48″ neighbours | Earth-sized USP; SPOC's offset not reproduced |
| 3 | TIC 149390648 | seeing-limited photometry of the crowded field | likely planet; crowding is the main remaining question |
| 4 | TIC 198412174 | adaptive-optics or speckle imaging; photometry resolving the 5.3″ neighbour | host ambiguous between the target and a possible companion |

## References

* Lissauer, J. J., et al. 2012, ApJ, 750, 112 (multiplicity boost for candidates in multi-candidate systems)
