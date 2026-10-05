"""Figure for the research note: folded TESS transits of the four candidates with the MCMC model.

    .venv/bin/python paper/make_figure.py      ->  paper/figure1.png
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tess_search import DATA, RESULTS, lightcurve, mcmc

CANDIDATES = [(32090583, "TOI-218 (TIC 32090583)"), (229689348, "TIC 229689348"), (149390648, "TIC 149390648"),
              (198412174, "TIC 198412174")]


def main():
    tic_tab = pd.read_parquet(DATA / "tic_15plus.parquet").set_index("tic")
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4), sharey=False)
    for ax, (tic, name) in zip(axes.ravel(), CANDIDATES):
        m = json.loads((RESULTS / "hardening" / "mcmc" / f"TIC{tic}.json").read_text())
        s = m["prior"]
        period = m["period"]
        t0 = s["t0_btjd"]["median"]
        t14 = s["t14_h"]["median"] / 24
        lc = mcmc.remask_detrend(lightcurve.prepare(tic), period, t0, t14)
        t, f, e, _ = mcmc.fold_bin(lc, period, t0, t14, bin_min=1.0)
        nb = 8                               # 8-minute bins for display
        n = len(t) // nb * nb
        tb = t[:n].reshape(-1, nb).mean(1)
        fb = f[:n].reshape(-1, nb).mean(1)
        eb = np.sqrt((e[:n].reshape(-1, nb) ** 2).sum(1)) / nb
        post = mcmc.TransitPosterior(t, f, e, period, t14, None)
        theta = [0.0, s["rp_rs"]["median"], s["b"]["median"], np.log10(s["rho_cgs"]["median"]),
                 s["u1"]["median"], s["u2"]["median"]]
        ax.errorbar(tb * 24, (fb - 1) * 1e6, eb * 1e6, fmt="o", ms=2.5, color="0.2", lw=0.6, capsize=0)
        ax.plot(t * 24, (post.curve(theta) - 1) * 1e6, color="C3", lw=1.5)
        r = s["rp_rearth"]
        ax.set_title(f"{name}\nP = {period:.4f} d, Rp = {r['median']:.2f} ± {0.5 * (r['plus'] + r['minus']):.2f} R$_\\oplus$",
                     fontsize=8.5)
        ax.set_xlim(-3 * t14 * 24, 3 * t14 * 24)
        ax.tick_params(labelsize=8)
    for ax in axes[1]:
        ax.set_xlabel("hours from mid-transit", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("relative flux (ppm)", fontsize=9)
    fig.tight_layout()
    out = Path(__file__).resolve().parent / "figure1.png"
    fig.savefig(out, dpi=200)
    print(out)


if __name__ == "__main__":
    main()
