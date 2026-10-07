"""Figures for the technical documents in docs/ and for the project site (every value is read from results/).

    .venv/bin/python docs/make_figures.py [names]            ->  docs/figures/*.png   (light, for the Markdown)
    .venv/bin/python docs/make_figures.py --dark [names]     ->  site/figures/*.png   (dark, for the website)

Style: one chart surface, recessive hairline grid, 2 px lines, >= 8 px markers, a fixed categorical
order (the Okabe-Ito blue, orange, green and pink, distinguishable under colour-vision deficiency),
a single ramp for the localization maps, and no second y-axis. Text is set in Inter when
data/fonts/Inter.ttf is present (the site's typeface), otherwise Helvetica Neue.
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
from matplotlib.colors import LinearSegmentedColormap

from tess_search import DATA, RESULTS

HARD = RESULTS / "hardening"
DOCS = Path(__file__).resolve().parent

THEMES = {
    # light: the Markdown documents (GitHub renders them on white)
    "light": dict(out=DOCS / "figures", surface="#ffffff", ink="#0b0b0b", ink2="#52514e", muted="#898781",
                  grid="#e6e5df", axis="#c3c2b7", band="#eeede7", err="#b4b3aa",
                  c=("#2a78d6", "#eb6834", "#1baf7a", "#eda100"),
                  # localization maps: allowed = dark blue, excluded = surface
                  # (position on the 0-8 sigma scale, colour): reaches the surface colour by ~6 sigma
                  ramp=[(0, "#0d366b"), (0.12, "#184f95"), (0.3, "#3987e5"), (0.5, "#9ec5f4"), (0.7, "#dceaf9"),
                        (0.85, "#f7fafd"), (1, "#ffffff")],
                  map_note="Shading: where the source can be\n(dark = allowed, light = excluded)."),
    # dark: the project site (black page, white Inter text)
    "dark": dict(out=DOCS.parent / "site" / "figures", surface="#000000", ink="#ffffff", ink2="#d4d4d4", muted="#8e8e8e",
                 grid="#1f1f1f", axis="#3a3a3a", band="#1a1a1a", err="#5e5e5e",
                 c=("#56b4e9", "#e69f00", "#2fc495", "#da8fc2"),
                 # localization maps: allowed = bright, excluded = black
                 ramp=[(0, "#f2f8fd"), (0.12, "#b9d9f5"), (0.3, "#56b4e9"), (0.5, "#1d5f8a"), (0.7, "#0a2436"),
                       (0.85, "#02090e"), (1, "#000000")],
                 map_note="Shading: where the source can be\n(bright = allowed, dark = excluded)."),
}
THEME = THEMES["light"]
OUT = THEME["out"]
SURFACE = INK = INK2 = MUTED = GRID = AXIS = BAND = ERR = C1 = C2 = C3 = C4 = ""
MAP_RAMP, MAP_NOTE = [], ""
SYS_ARCSEC = 1.5
R68 = 1.5151947600179898

CANDS = [(32090583, "TOI-218 (TIC 32090583)"), (229689348, "TIC 229689348"), (149390648, "TIC 149390648"),
         (198412174, "TIC 198412174")]


def style(theme="light"):
    global THEME, OUT, SURFACE, INK, INK2, MUTED, GRID, AXIS, BAND, ERR, C1, C2, C3, C4, MAP_RAMP, MAP_NOTE
    THEME = THEMES[theme]
    OUT = THEME["out"]
    SURFACE, INK, INK2, MUTED = THEME["surface"], THEME["ink"], THEME["ink2"], THEME["muted"]
    GRID, AXIS, BAND, ERR = THEME["grid"], THEME["axis"], THEME["band"], THEME["err"]
    C1, C2, C3, C4 = THEME["c"]
    MAP_RAMP, MAP_NOTE = THEME["ramp"], THEME["map_note"]
    inter = DATA / "fonts" / "Inter.ttf"
    if inter.exists():
        from matplotlib import font_manager
        font_manager.fontManager.addfont(str(inter))
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "text.color": INK, "font.size": 9, "axes.titlesize": 9.5, "axes.titleweight": "regular", "axes.titlecolor": INK,
        "axes.titlelocation": "left", "axes.labelcolor": INK2, "axes.labelsize": 9,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
        "axes.axisbelow": True, "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2,
        "ytick.labelcolor": INK2, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.frameon": False,
        "legend.fontsize": 8, "lines.linewidth": 2.0, "lines.solid_capstyle": "round",
        "lines.solid_joinstyle": "round", "errorbar.capsize": 0,
    })


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name, dpi=200, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(OUT / name)


def load(p):
    return json.loads(Path(p).read_text())


# ------------------------------------------------------------------ 1. sample


def fig_sample():
    t = pd.read_csv(DATA / "targets.csv")
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.6))
    bins = np.arange(19.5, t.n_sectors.max() + 1.5)
    ax[0].hist(t.n_sectors, bins=bins, color=C1, rwidth=0.8)
    ax[0].set_xlabel("sectors of 2-minute data")
    ax[0].set_ylabel("stars")
    ax[0].set_title(f"Observing baseline: {len(t):,} M dwarfs, median {t.n_sectors.median():.0f} sectors")
    ax[1].scatter(t.Teff, t.rad, s=6, color=C1, alpha=0.55, linewidths=0)
    ax[1].set_xlabel("effective temperature (K)")
    ax[1].set_ylabel("stellar radius (R$_\\odot$)")
    ax[1].set_title("Stellar parameters (TIC v8.2)")
    ax[1].invert_xaxis()
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig01_sample.png")


# ------------------------------------------------------------------ 2. funnel


def fig_funnel():
    s = load(RESULTS / "summary.json")
    c = pd.read_csv(RESULTS / "candidates.csv")
    passed = len(c)
    novel = int(c.novelty.isin(["new", "SPOC TCE only (never promoted)"]).sum())
    novel_strong = int(((c.verdict == "candidate") & c.novelty.isin(["new", "SPOC TCE only (never promoted)"])).sum())
    cs = pd.read_csv(RESULTS / "candidates" / "summary.csv")
    final = int(cs.verdict.str.startswith("Candidate").sum())
    stages = [("periodic signals detected", s["signals_total"]),
              ("strong enough to vet (SNR \u2265 6)", s["signals_vetted"]),
              ("passed light-curve vetting", passed),
              ("not a known planet, TOI or CTOI", novel),
              ("... at full candidate strength", novel_strong),
              ("on target after pixel and statistical checks", final)]
    fig, ax = plt.subplots(figsize=(7.4, 2.5))
    y = np.arange(len(stages))[::-1]
    vals = [v for _, v in stages]
    ax.barh(y, vals, height=0.5, color=C1)
    ax.set_xscale("log")
    ax.set_xlim(1, 4000)
    ax.set_xticks([1, 10, 100, 1000])
    ax.set_xticklabels(["1", "10", "100", "1,000"])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for yi, (lab, v) in zip(y, stages):
        ax.text(v * 1.12, yi, f"{v:,}", va="center", ha="left", color=INK, fontsize=8.5)
    ax.set_yticks(y)
    ax.set_yticklabels([lab for lab, _ in stages])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("number of signals (log scale)")
    ax.set_title("From 1,560 periodic signals to four candidates")
    fig.tight_layout()
    save(fig, "fig02_funnel.png")


# ------------------------------------------------------------------ 3. completeness


def wilson(k, n, z=1.0):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, c - h, c + h


def fig_completeness():
    from tess_search import vetting

    rows = [json.loads(line) for line in (RESULTS / "reliability" / "inject.jsonl").read_text().splitlines() if line.strip()]
    for r in rows:
        r["passed"] = False
        if r.get("recovered"):
            sg = next((x for x in r["signals"] if abs(x["period"] / r["period"] - 1) < 0.005), None)
            if sg is not None and sg.get("vet"):
                r["passed"] = vetting.classify(sg["vet"])[0] in ("candidate", "weak candidate")
    df = pd.DataFrame(rows)
    edges = np.array([0.6, 1.0, 1.5, 2.0, 3.0, 4.0])
    mids = 0.5 * (edges[1:] + edges[:-1])
    groups = [("P < 5 d", df.period < 5, C1), ("5-15 d", (df.period >= 5) & (df.period < 15), C2),
              ("15-40 d", df.period >= 15, C3)]
    fig, ax = plt.subplots(figsize=(7.4, 3.0))
    for k, (lab, m, col) in enumerate(groups):
        p, lo, hi = [], [], []
        for a, b in zip(edges[:-1], edges[1:]):
            sel = m & (df.rp >= a) & (df.rp < b)
            pk = wilson(int(df[sel].passed.sum()), int(sel.sum()))
            p.append(pk[0]); lo.append(pk[1]); hi.append(pk[2])
        x = mids + (k - 1) * 0.04
        ax.plot(x, p, color=col, marker="o", ms=5, mec=SURFACE, mew=1.5, label=lab, zorder=3)
        ax.vlines(x, lo, hi, color=col, lw=1.2, alpha=0.6, zorder=2)
        ax.text(x[-1] + 0.08, p[-1] + (0, 0.035, -0.035)[k], lab, color=INK2, va="center", fontsize=8)
    tot = [wilson(int(df[(df.rp >= a) & (df.rp < b)].passed.sum()), int(((df.rp >= a) & (df.rp < b)).sum()))[0]
           for a, b in zip(edges[:-1], edges[1:])]
    ax.plot(mids, tot, color=MUTED, lw=1.2, ls=(0, (1, 2)), zorder=1)
    ax.text(mids[3], tot[3] + 0.045, "all periods", color=MUTED, va="center", ha="center", fontsize=8)
    ax.set_xticks(edges)
    ax.set_xlim(0.55, 4.45)
    ax.set_ylim(0, 1.02)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax.set_xlabel("injected planet radius (R$_\\oplus$)")
    ax.set_ylabel("found and kept")
    ax.set_title(f"Completeness from {len(df)} planets injected into real light curves (1-sigma binomial intervals)")
    fig.tight_layout()
    save(fig, "fig03_completeness.png")


# ------------------------------------------------------------------ 4. localization validation


def sigma_from_dchi2(d):
    from scipy import stats

    p = stats.chi2.sf(np.clip(d, 0, None), 2)
    return np.where(p > 1e-300, stats.norm.isf(np.clip(p, 1e-300, 1) / 2), 37.0)


def fig_localization_validation():
    loc = pd.read_csv(HARD / "localization_summary.csv")
    known = loc[loc.kind.isin(["planet", "injection"]) & loc.miss_arcsec.notna() & (loc.chi2_red < 2)]
    s_stat = known.offset_err_arcsec / R68
    pull = np.sort((known.miss_arcsec / np.sqrt(s_stat ** 2 + SYS_ARCSEC ** 2)).to_numpy())
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.0), gridspec_kw={"width_ratios": [1, 1.25]})
    x = np.linspace(0, 4, 200)
    ax[0].plot(x, 1 - np.exp(-x ** 2 / 2), color=MUTED, lw=1.5, ls=(0, (3, 2)), label="expected (2-D Gaussian)")
    ax[0].step(pull, np.arange(1, len(pull) + 1) / len(pull), where="post", color=C1, label=f"{len(pull)} known sources")
    ax[0].set_xlabel("position error / predicted 1-sigma")
    ax[0].set_ylabel("cumulative fraction")
    ax[0].set_title("Position errors are conservative (1.5\" floor)")
    ax[0].set_xlim(0, 4)
    ax[0].set_ylim(0, 1.02)
    ax[0].legend(loc="lower right")

    groups = [("confirmed planets", loc[loc.kind == "planet"], C1),
              ("TFOP nearby EBs", loc[loc.kind == "neb"], C2),
              ("this work's signals", loc[loc.kind == "candidate"], C3)]
    for k, (lab, d, col) in enumerate(groups):
        vals = np.clip(d.target_sigma_total.to_numpy(), 0.06, 40)
        yy = np.full(len(vals), 2 - k) + np.linspace(-0.12, 0.12, len(vals)) if len(vals) > 1 else [2 - k]
        ax[1].scatter(vals, yy, s=34, color=col, edgecolor=SURFACE, linewidths=1.5, zorder=3)
        if lab == "this work's signals":
            for v, yv, name in zip(vals, yy, d.label):
                if v > 3:
                    ax[1].annotate(name.replace("TIC", "TIC "), (v, yv), xytext=(0, -11), textcoords="offset points",
                                   ha="center", fontsize=7.5, color=INK2)
    ax[1].axvline(3, color=MUTED, lw=1, ls=(0, (3, 2)))
    ax[1].text(3.2, 2.45, "3 sigma", color=MUTED, fontsize=7.5)
    ax[1].set_xscale("log")
    ax[1].set_xlim(0.035, 60)
    ax[1].set_xticks([0.1, 1, 3, 10, 37])
    ax[1].set_xticklabels(["0.1", "1", "3", "10", "$\\geq$37"])
    ax[1].set_yticks([2, 1, 0])
    ax[1].set_yticklabels([g[0] for g in groups])
    ax[1].set_ylim(-0.6, 2.6)
    ax[1].grid(axis="y", visible=False)
    ax[1].set_xlabel("target excluded as the source at (sigma)")
    ax[1].set_title("Planets stay on target; nearby EBs do not")
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig04_localization_validation.png")


# ------------------------------------------------------------------ 5. folded transits


def fig_transits():
    from tess_search import lightcurve, mcmc

    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.0))
    for ax, (tic, name) in zip(axes.ravel(), CANDS):
        m = load(HARD / "mcmc" / f"TIC{tic}.json")
        s = m["prior"]
        period = m["period"]
        t0 = s["t0_btjd"]["median"]
        t14 = s["t14_h"]["median"] / 24
        lc = mcmc.remask_detrend(lightcurve.prepare(tic), period, t0, t14)
        t, f, e, _ = mcmc.fold_bin(lc, period, t0, t14, bin_min=1.0)
        nb = 8
        n = len(t) // nb * nb
        tb, fb = t[:n].reshape(-1, nb).mean(1), f[:n].reshape(-1, nb).mean(1)
        eb = np.sqrt((e[:n].reshape(-1, nb) ** 2).sum(1)) / nb
        post = mcmc.TransitPosterior(t, f, e, period, t14, None)
        th = [0.0, s["rp_rs"]["median"], s["b"]["median"], np.log10(s["rho_cgs"]["median"]),
              s["u1"]["median"], s["u2"]["median"]]
        ax.vlines(tb * 24, (fb - eb - 1) * 1e6, (fb + eb - 1) * 1e6, color=ERR, lw=1, zorder=1)
        ax.scatter(tb * 24, (fb - 1) * 1e6, s=10, color=INK2, zorder=2, linewidths=0)
        ax.plot(t * 24, (post.curve(th) - 1) * 1e6, color=C1, zorder=3)
        r = s["rp_rearth"]
        ax.set_title(f"{name}\nP = {period:.4f} d    R$_p$ = {r['median']:.2f} $\\pm$ "
                     f"{0.5 * (r['plus'] + r['minus']):.2f} R$_\\oplus$")
        ax.set_xlim(-3 * t14 * 24, 3 * t14 * 24)
    for ax in axes[1]:
        ax.set_xlabel("hours from mid-transit")
    for ax in axes[:, 0]:
        ax.set_ylabel("relative flux (ppm)")
    fig.tight_layout(h_pad=1.6, w_pad=2.0)
    save(fig, "fig07_transits.png")


# ------------------------------------------------------------------ 6. localization maps


def fig_localization_maps():
    panels = CANDS + [(294053492, "TIC 294053492 (false positive)")]
    cmap = LinearSegmentedColormap.from_list("likelihood", MAP_RAMP)
    fig, axes = plt.subplots(2, 3, figsize=(7.6, 5.4))
    for ax, (tic, name) in zip(axes.ravel(), panels):
        d = load(HARD / "localize" / f"TIC{tic}.json")
        offs = np.array(d["grid"]["offsets_arcsec"])
        s_stat = d["offset_err_arcsec"] / R68
        dchi2 = np.array(d["grid"]["dchi2"]) * s_stat ** 2 / (s_stat ** 2 + SYS_ARCSEC ** 2)
        sig = sigma_from_dchi2(dchi2)
        ax.imshow(np.clip(sig, 0, 8), origin="lower", extent=[offs[0], offs[-1], offs[0], offs[-1]], cmap=cmap,
                  vmin=0, vmax=8, interpolation="bilinear")
        ax.contour(offs, offs, sig, levels=[3], colors=[INK2], linewidths=[0.8], linestyles=[(0, (3, 2))])
        tmag_t = next(s["tmag"] for s in d["stars"] if s["is_target"])
        for s in d["stars"]:
            if abs(s["east"]) < 60 and abs(s["north"]) < 60 and not s["is_target"]:
                size = float(np.clip(90 * 10 ** (-0.4 * (s["tmag"] - tmag_t)), 6, 160))
                ax.scatter(s["east"], s["north"], s=size, facecolors="none", edgecolors=INK2, linewidths=0.8)
        ax.plot(0, 0, marker="+", color=INK, ms=11, mew=1.6)
        ax.plot(d["offset_east"], d["offset_north"], marker="x", color=C2, ms=8, mew=2)
        ax.set_xlim(60, -60)
        ax.set_ylim(-60, 60)
        ax.set_xticks([-50, 0, 50])
        ax.set_yticks([-50, 0, 50])
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_title(f"{name}\noffset {d['offset_arcsec']:.1f} $\\pm$ {d['offset_err_total_arcsec']:.1f}\"; "
                     f"target at {d['target_sigma_total']:.1f}$\\sigma$", fontsize=8.5)
    for ax in axes[1]:
        ax.set_xlabel("arcsec east")
    for ax in axes[:, 0]:
        ax.set_ylabel("arcsec north")
    leg = axes.ravel()[-1]
    leg.axis("off")
    leg.plot([], [], marker="+", color=INK, ms=11, mew=1.6, ls="none", label="target")
    leg.plot([], [], marker="x", color=C2, ms=8, mew=2, ls="none", label="best-fit source position")
    leg.scatter([], [], s=60, facecolors="none", edgecolors=INK2, label="Gaia DR3 star (size ~ brightness)")
    leg.plot([], [], color=INK2, lw=0.8, ls=(0, (3, 2)), label="3-sigma region")
    leg.legend(loc="center left", fontsize=8)
    leg.text(0.0, 0.08, MAP_NOTE + "\n1.5\" systematic floor included.", transform=leg.transAxes, fontsize=7.5,
             color=INK2)
    fig.tight_layout(h_pad=1.2, w_pad=1.0)
    save(fig, "fig06_localization_maps.png")


# ------------------------------------------------------------------ 7. pixel depth ratio


def fig_depth_ratio():
    from tess_search.localize import pixel_depth_ratio

    loc = pd.read_csv(HARD / "localization_summary.csv")
    weak = pd.read_csv(HARD / "weak_recheck.csv")
    rows = []
    for _, r in loc[loc.kind.isin(["planet", "candidate"])].iterrows():
        rec = load(HARD / "localize" / f"{r.label}.json")
        ratio, err, _, _ = pixel_depth_ratio(rec)
        rows.append((r.kind, r.label, ratio, err, ""))
    for _, r in weak.iterrows():
        rows.append(("weak", f"TIC{r.tic}_{r.signal}", r.pixel_depth_ratio, r.ratio_err, r["class"]))
    df = pd.DataFrame(rows, columns=["kind", "label", "ratio", "err", "cls"])
    df = df[~df.label.eq("TIC294053492")]
    fig, ax = plt.subplots(figsize=(7.4, 3.0))
    ax.axhspan(0.85, 1.29, color=BAND, lw=0, zorder=0)
    ax.axhline(1, color=MUTED, lw=1, zorder=1)
    x0 = 0
    ticks, ticklabels = [], []
    for kind, lab, col in (("planet", "confirmed planets", C1), ("candidate", "candidates", C3),
                           ("weak", "weak signals", C2)):
        d = df[df.kind == kind].sort_values("ratio")
        x = x0 + np.arange(len(d))
        ax.vlines(x, d.ratio - d.err, d.ratio + d.err, color=col, lw=1.2, alpha=0.5)
        if kind == "weak":
            for cls, marker in (("consistent", "o"), ("marginal", "s"), ("not in pixels", "v")):
                m = (d.cls == cls).to_numpy()
                ax.scatter(x[m], d.ratio.to_numpy()[m], s=30, color=col, marker=marker, edgecolor=SURFACE,
                           linewidths=1.2, zorder=3, label=f"weak: {cls} ({m.sum()})")
        else:
            ax.scatter(x, d.ratio, s=30, color=col, edgecolor=SURFACE, linewidths=1.2, zorder=3,
                       label=f"{lab} ({len(d)})")
        ticks.append(x.mean())
        ticklabels.append(lab)
        x0 += len(d) + 3
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels)
    ax.grid(axis="x", visible=False)
    ax.set_ylabel("pixel light loss / light-curve depth")
    ax.set_ylim(-0.3, 2.0)
    ax.set_title("Does the target lose, in the pixels, the light the light curve says it does?  "
                 "(band: range for confirmed planets)")
    ax.legend(loc="upper right", ncol=2, fontsize=7.5)
    fig.tight_layout()
    save(fig, "fig05_depth_ratio.png")


# ------------------------------------------------------------------ 8. SPOC period drift


def fig_spoc_period():
    from astropy.timeseries import BoxLeastSquares

    from tess_search import lightcurve
    from tess_search.lightcurve import robust_std

    chk = load(HARD / "spoc_period_check.json")
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.8))
    for ax, tic in zip(axes, ("229689348", "198412174")):
        c = chk[tic]
        lc = lightcurve.prepare(int(tic))
        span = np.ptp(lc.time)
        p0 = c["our_period"]
        runs = c["runs"]
        lo = min([r["spoc_period"] for r in runs] + [p0]) - 3e-5
        hi = max([r["spoc_period"] for r in runs] + [p0]) + 3e-5
        dur = 0.55 / 24 if p0 < 1 else 0.6 / 24
        grid = np.arange(lo, hi, p0 * dur / span / 8)
        bls = BoxLeastSquares(lc.time, lc.flux, dy=np.full(len(lc.time), robust_std(lc.flux - 1)))
        snr = np.asarray(bls.power(grid, dur, objective="snr", oversample=20).depth_snr)
        ax.plot((grid - p0) * 1e5, snr, color=C1, lw=1.5)
        ax.axvline(0, color=INK, lw=1)
        ax.text(0, snr.max() * 1.19, "this work", fontsize=7.5, color=INK, ha="center", va="center",
                bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.5))
        lines = []
        for j, r in enumerate(sorted(runs, key=lambda r: r["spoc_period"])):
            dx = (r["spoc_period"] - p0) * 1e5
            tag = r["run"].replace("s00", "s")
            ax.axvline(dx, color=C2, lw=1, ls=(0, (3, 2)))
            ax.text(dx, snr.max() * (1.08 - 0.09 * j), f" {tag}", fontsize=7, color=INK2, va="center",
                    ha="left" if dx < 0.5 else "left")
            lines.append(f"{tag}: SPOC SNR {r['spoc_snr']:.1f}, here {r['our_snr_at_spoc_period']:.1f}")
        ax.text(0.98, 0.04, "\n".join(lines), transform=ax.transAxes, ha="right", va="bottom", fontsize=7,
                color=INK2, linespacing=1.4, bbox=dict(facecolor=SURFACE, edgecolor="none", pad=2.5))
        ax.set_xlabel("trial period minus this work's period (10$^{-5}$ d)")
        ax.set_title(f"TIC {tic}")
        ax.set_ylim(0, snr.max() * 1.25)
    axes[0].set_ylabel("box SNR of this light curve")
    fig.tight_layout(w_pad=2.0)
    save(fig, "fig08_spoc_period_drift.png")


# ------------------------------------------------------------------ 9. TRICERATOPS


def scen_group(code):
    if code == "TP":
        return 0
    if code in ("PTP", "PEB", "PEBx2P", "STP", "SEB", "SEBx2P"):
        return 1
    if code in ("NTP", "NEB", "NEBx2P"):
        return 2
    return 3


def fig_triceratops():
    labels = ["planet on target", "unresolved bound companion", "resolved neighbour", "other (EB, background)"]
    colors = [C1, C2, C3, C4]
    rows = []
    for tic, name in CANDS:
        for suffix, tag in (("", "TESS only"), ("_cleared", "cleared"), ("_cleared_cc", "cleared + imaging")):
            if suffix == "_cleared_cc" and not (HARD / "triceratops" / f"TIC{tic}_result{suffix}.json").exists():
                continue
            r = load(HARD / "triceratops" / f"TIC{tic}_result{suffix}.json")
            share = np.zeros(4)
            for key, v in r["top_scenarios"]:
                share[scen_group(key.split(":")[0])] += v
            share[3] += max(0.0, 1 - share.sum())
            rows.append((f"{name.split(' (')[0]}", tag, share, r["FPP_mean"]))
    fig, ax = plt.subplots(figsize=(7.4, 3.9))
    y, yy = [], 0.0
    for i, (name, tag, share, fpp) in enumerate(rows):
        if i and name != rows[i - 1][0]:
            yy -= 0.6                      # gap between candidates
        y.append(yy)
        left = 0.0
        for k in range(4):
            w = share[k]
            if w > 0:
                ax.barh(yy, w - (0.004 if w > 0.01 else 0), left=left, height=0.62, color=colors[k],
                        label=labels[k] if i == 0 else None)
            left += w
        ax.text(1.01, yy, f"FPP {fpp:.2f}" if fpp >= 0.01 else f"FPP {fpp:.0e}".replace("e-0", "e-"), va="center", fontsize=7.5,
                color=INK2)
        yy -= 1.0
    ax.set_yticks(y)
    ax.set_yticklabels([f"{n}  ·  {t}" for n, t, _, _ in rows], fontsize=7.5)
    ax.set_xlim(0, 1.1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "25%", "50%", "75%", "100%"])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("TRICERATOPS scenario probability")
    ax.set_title("What else could the dip be?  TESS alone, pixel-excluded neighbours cleared, and imaging")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=4, fontsize=7.5)
    fig.tight_layout()
    save(fig, "fig09_triceratops.png")


# ------------------------------------------------------------------ 10. weak-signal periods


def fig_weak_periods():
    from tess_search import vetting

    weak = pd.read_csv(HARD / "weak_recheck.csv")
    inv = [json.loads(line) for line in (RESULTS / "reliability" / "invert.jsonl").read_text().splitlines() if line.strip()]
    inv_p = [s["period"] for r in inv for s in r["signals"] if s.get("vet") and vetting.classify(s["vet"])[0] == "weak candidate"]
    fig, ax = plt.subplots(figsize=(7.4, 2.3))
    ax.axvspan(36, 41, color=BAND, lw=0)
    ax.text(38.3, 1.42, "36-41 d", fontsize=7.5, color=INK2, ha="center")
    rng = np.random.default_rng(3)
    ax.scatter(weak.period, 1 + rng.uniform(-0.22, 0.22, len(weak)), s=30, color=C1, edgecolor=SURFACE,
               linewidths=1.2, zorder=3)
    ax.scatter(inv_p, rng.uniform(-0.22, 0.22, len(inv_p)), s=30, color=C2, edgecolor=SURFACE, linewidths=1.2,
               zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(0.4, 45)
    ax.set_xticks([0.5, 1, 2, 5, 10, 20, 40])
    ax.set_xticklabels(["0.5", "1", "2", "5", "10", "20", "40"])
    ax.set_yticks([1, 0])
    ax.set_yticklabels([f"weak candidates,\nreal data ({len(weak)})", f"false weak candidates,\ninverted data ({len(inv_p)})"])
    ax.set_ylim(-0.45, 1.55)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("period (d)")
    ax.set_title("Weak signals pile up where few transits are observed, in real and in inverted data alike")
    fig.tight_layout()
    save(fig, "fig10_weak_periods.png")


if __name__ == "__main__":
    args = sys.argv[1:]
    style("dark" if "--dark" in args else "light")
    only = [a for a in args if a != "--dark"]
    for name, fn in [("sample", fig_sample), ("funnel", fig_funnel), ("completeness", fig_completeness),
                     ("locval", fig_localization_validation), ("transits", fig_transits),
                     ("locmaps", fig_localization_maps), ("ratio", fig_depth_ratio), ("spoc", fig_spoc_period),
                     ("tri", fig_triceratops), ("weak", fig_weak_periods)]:
        if not only or name in only:
            fn()
