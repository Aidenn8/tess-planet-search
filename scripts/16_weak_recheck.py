"""Step 16: the weak candidates, re-checked in the pixels.

For each weak signal (localized by `12_localize.py --weak`): does the target lose, in the
difference images, the light the light-curve depth predicts? Classes:
  consistent     ratio within 2 sigma of 1 and the light loss detected at > 3 sigma in the pixels
  not in pixels  ratio more than 3 sigma below 1: the light-curve dip is not reproduced in the
                 pixels (most likely a light-curve artefact)
  marginal       in between
Errors include the 12% scatter measured on confirmed planets. A search picks noise peaks with
inflated depths, so even real weak signals would tend to come out below 1.
Writes results/hardening/weak_recheck.csv.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from tess_search import RESULTS
from tess_search.localize import pixel_depth_ratio

HARD = RESULTS / "hardening"

if __name__ == "__main__":
    loc = pd.read_csv(HARD / "localization_weak.csv")
    rows = []
    for lab in loc.label:
        rec = json.loads((HARD / "localize" / f"{lab}.json").read_text())
        ratio, err, z1, z0 = pixel_depth_ratio(rec)
        cls = "consistent" if (z1 > -2 and z0 > 3) else ("not in pixels" if z1 < -3 else "marginal")
        tic, sig = lab[3:].split("_")
        rows.append({"tic": int(tic), "signal": int(sig), "period": rec["period"], "depth_ppm": rec["depth_ppm"],
                     "pixel_depth_ratio": ratio, "ratio_err": err, "sigma_from_1": z1, "pixel_detection_sigma": z0,
                     "class": cls, "offset_arcsec": rec["offset_arcsec"],
                     "target_sigma_total": rec["target_sigma_total"], "reduced_chi2": rec["chi2_red"]})
    df = pd.DataFrame(rows).sort_values("sigma_from_1", ascending=False)
    df.to_csv(HARD / "weak_recheck.csv", index=False)
    print(df.round(2).to_string(index=False))
    print(df["class"].value_counts().to_dict())
