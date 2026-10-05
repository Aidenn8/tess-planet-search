# Technical write-up

A complete account of the search, its calibration and its results, in eight chapters. Each chapter
stands on its own; together they document every step from raw TESS data to the four candidates.

| | chapter | contents |
|---|---|---|
| 1 | [Data and search](01-data-and-search.md) | sample selection, light-curve cleaning and detrending, the stacked seasonal BLS search, catalogue crossmatch |
| 2 | [Vetting](02-vetting.md) | red-noise-aware and null-calibrated tests, thresholds, verdicts, what happened to all 1,139 vetted signals |
| 3 | [Sensitivity and reliability](03-sensitivity-and-reliability.md) | recovery of known planets, 300 injected planets, 200 inverted light curves |
| 4 | [Pixel-level localization](04-pixel-level-localization.md) | difference imaging with the SPOC PRF, systematic-error calibration, validation on planets, nearby eclipsing binaries and injected eclipses |
| 5 | [Transit fits and false-positive probabilities](05-transit-fits-and-false-positive-probabilities.md) | MCMC transit model, validation on TOI-700 d and L 98-59 c, NASA's pipeline history, TRICERATOPS |
| 6 | [Candidates](06-candidates.md) | the four candidates and one false positive in detail, with follow-up priorities |
| 7 | [Weak signals](07-weak-signals.md) | 21 weak signals, their expected false-alarm rate, and a pixel re-check |
| 8 | [Limitations and reproducibility](08-limitations-and-reproducibility.md) | caveats, environment, commands, compute, repository map, acknowledgements |

Figures are generated from the result files by [`make_figures.py`](make_figures.py).
Per-candidate dossiers are in [`../results/candidates/`](../results/candidates/); a generated summary of
every number is in [`../REPORT.md`](../REPORT.md).
