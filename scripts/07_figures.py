"""Step 7: figures for the write-up (results/figures/*.png).

    .venv/bin/python scripts/07_figures.py

  sample.png         the searched stars: temperature, brightness, sectors, noise
  completeness.png   injection-recovery: fraction of fake planets found vs period and size
  snr_recovery.png   recovery fraction vs expected signal-to-noise
  verdicts.png       what happened to every detected signal
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from tess_search import DATA, RESULTS  # noqa: E402

FIG = RESULTS / "figures"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


def sample_figure():
    targets = pd.read_csv(DATA / "targets.csv")
    recs = [json.loads(p.read_text()) for p in (RESULTS / "search").glob("*.json")]
    noise = pd.DataFrame([{"tic": r["tic"], "noise": r["diagnostics"]["noise_ppm"]} for r in recs])
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.8))
    ax[0].scatter(targets.Teff, targets.rad, s=4, c=targets.Tmag, cmap="viridis_r")
    ax[0].set_xlabel("effective temperature (K)")
    ax[0].set_ylabel("radius (R_sun)")
    ax[0].set_title(f"{len(targets):,} M dwarfs (darker = fainter)", loc="left")
    ax[1].hist(targets.n_sectors, bins=np.arange(19.5, 46.5, 1), color="0.3")
    ax[1].set_xlabel("sectors of 2-minute data")
    ax[1].set_ylabel("stars")
    ax[1].set_title(f"median {int(targets.n_sectors.median())} sectors (~2 years) per star", loc="left")
    if len(noise):
        m = targets.merge(noise, on="tic")
        ax[2].scatter(m.Tmag, m.noise, s=4, color="0.3")
        ax[2].set_yscale("log")
        ax[2].set_xlabel("TESS magnitude")
        ax[2].set_ylabel("noise per 2-min point (ppm)")
        ax[2].set_title("photometric noise after detrending", loc="left")
    fig.tight_layout()
    fig.savefig(FIG / "sample.png", dpi=130)
    plt.close(fig)


def load_jsonl(name):
    p = RESULTS / "reliability" / f"{name}.jsonl"
    if not p.exists():
        return pd.DataFrame()
    return pd.DataFrame([json.loads(l) for l in p.read_text().splitlines() if l.strip()])


def completeness_figure():
    inj = load_jsonl("inject")
    if inj.empty:
        return None
    from tess_search import vetting

    def passed(r):
        # final vetting rules re-applied to the stored metrics of the matching detection
        if not r["recovered"]:
            return False
        s = next((s for s in r["signals"] if abs(s["period"] / r["period"] - 1) < 0.005), None)
        v = vetting.classify(s["vet"])[0] if s and s.get("vet") else (s or {}).get("verdict")
        return v in ("candidate", "weak candidate")

    inj["passed"] = [passed(r) for r in inj.to_dict("records")]
    pb = np.exp(np.linspace(np.log(0.5), np.log(40), 7))
    rb = np.linspace(0.6, 4.0, 7)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for a, col, title in ((ax[0], "recovered", "detected at the right period"),
                          (ax[1], "passed", "detected AND passed vetting")):
        grid = np.full((len(rb) - 1, len(pb) - 1), np.nan)
        cnt = np.zeros_like(grid)
        for i in range(len(rb) - 1):
            for j in range(len(pb) - 1):
                m = (inj.rp >= rb[i]) & (inj.rp < rb[i + 1]) & (inj.period >= pb[j]) & (inj.period < pb[j + 1])
                cnt[i, j] = m.sum()
                if m.sum() >= 3:
                    grid[i, j] = inj.loc[m, col].mean()
        im = a.pcolormesh(pb, rb, grid, vmin=0, vmax=1, cmap="magma")
        for i in range(len(rb) - 1):
            for j in range(len(pb) - 1):
                if np.isfinite(grid[i, j]):
                    a.text(np.sqrt(pb[j] * pb[j + 1]), 0.5 * (rb[i] + rb[i + 1]),
                           f"{100 * grid[i, j]:.0f}%\n(n={int(cnt[i, j])})", ha="center", va="center", fontsize=7,
                           color="w" if grid[i, j] < 0.6 else "k")
        a.set_xscale("log")
        a.set_xlabel("orbital period (days)")
        a.set_ylabel("planet radius (Earth radii)")
        a.set_title(title, loc="left")
    fig.colorbar(im, ax=ax, label="fraction of injected planets")
    fig.savefig(FIG / "completeness.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # expected white-noise SNR = depth / noise * sqrt(points in transit), with the
    # number of in-transit points estimated from the star's sector count
    from tess_search.inject import expected_duration
    targets = pd.read_csv(DATA / "targets.csv")[["tic", "n_sectors"]]
    inj = inj.merge(targets, on="tic", how="left")
    inj["dur"] = [expected_duration(r.period, r.r_star, r.r_star, r.b) for r in inj.itertuples()]
    observed_days = inj.n_sectors * 27.4 * 0.85  # ~85% of each sector is usable 2-min data
    n_in = observed_days / inj.period * inj.dur / (2 / 1440)
    inj["snr_expected"] = inj.depth_true / (inj.noise_ppm * 1e-6) * np.sqrt(np.maximum(n_in, 1))
    bins = np.array([0, 4, 6, 8, 10, 12, 15, 20, 30, 50, 1e4])
    mid, f_rec, f_pass, n = [], [], [], []
    for a_, b_ in zip(bins[:-1], bins[1:]):
        m = (inj.snr_expected >= a_) & (inj.snr_expected < b_)
        if m.sum() >= 3:
            mid.append(min(np.sqrt(max(a_, 1) * b_), 80))
            f_rec.append(inj.loc[m, "recovered"].mean())
            f_pass.append(inj.loc[m, "passed"].mean())
            n.append(m.sum())
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(mid, f_rec, "o-", label="detected")
    ax.plot(mid, f_pass, "s-", label="detected and passed vetting")
    ax.set_xscale("log")
    ax.set_xlabel("expected signal-to-noise (white-noise estimate)")
    ax.set_ylabel("fraction recovered")
    ax.axvline(7.3, color="0.6", ls=":")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "snr_recovery.png", dpi=130)
    plt.close(fig)
    return inj


def verdict_figure():
    p = RESULTS / "all_signals.csv"
    if not p.exists():
        return
    sig = pd.read_csv(p)
    sig = sig[sig.verdict != "not vetted"]
    totals = sig.verdict.value_counts()
    sig = sig[sig.depth_ppm > 1]  # zero/negative fitted depths cannot be drawn on a log axis
    fig, ax = plt.subplots(figsize=(8, 4.5))
    colors = {"false positive": "0.7", "below threshold": "0.85", "weak candidate": "tab:orange",
              "candidate": "tab:blue"}
    for v, c in colors.items():
        m = sig.verdict == v
        ax.scatter(sig.period[m], sig.depth_ppm[m].clip(1, 1e5), s=8 if v in ("false positive", "below threshold") else 22,
                   color=c, label=f"{v} ({totals.get(v, 0)})", zorder=3 if "candidate" in v else 1)
    known = sig.novelty.fillna("").str.startswith("known")
    ax.scatter(sig.period[known], sig.depth_ppm[known].clip(1, 1e5), s=60, facecolor="none", edgecolor="k",
               label=f"matches known planet/TOI ({known.sum()})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("period (days)")
    ax.set_ylabel("depth (ppm)")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "verdicts.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    sample_figure()
    completeness_figure()
    verdict_figure()
    print("figures written to", FIG)
