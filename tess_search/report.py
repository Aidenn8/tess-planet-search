"""Diagnostic plot sheet for one detection (one PNG per signal).

Panels:
  A  full light curve (raw, with the detrending curve) and transit times
  B  stacked periodogram (SDE vs period) with the detected period marked
  C  folded transit with binned points and the trapezoid fit
  D  odd vs even transits side by side
  E  whole orbit folded (looking for a secondary eclipse)
  F  depth of each individual transit over time
  G  text: metrics, flags, catalogue matches
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .lightcurve import bin_lightcurve  # noqa: E402
from .vetting import phase_offset, trapezoid  # noqa: E402


def binned(x, y, width):
    order = np.argsort(x)
    x, y = x[order], y[order]
    edges = np.arange(x.min(), x.max() + width, width)
    idx = np.clip(np.digitize(x, edges) - 1, 0, len(edges) - 2)
    cnt = np.bincount(idx, minlength=len(edges) - 1)
    ok = cnt > 2
    yb = np.bincount(idx, weights=y, minlength=len(edges) - 1)[ok] / cnt[ok]
    y2 = np.bincount(idx, weights=y ** 2, minlength=len(edges) - 1)[ok] / cnt[ok]
    err = np.sqrt(np.maximum(y2 - yb ** 2, 0) / cnt[ok])
    return 0.5 * (edges[1:] + edges[:-1])[ok], yb, err


def plot_detection(lc, det, vet, verdict, reasons, matches, periods=None, sde_curve=None,
                   path=None, title_extra=""):
    P, t0 = det["period"], det["t0"]
    trap = vet.get("trap") or {}
    dur = trap.get("t14", det["duration"]) if trap else det["duration"]
    dt = phase_offset(lc.time, P, t0)

    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(4, 3, height_ratios=[1, 1, 1, 0.9], hspace=0.45, wspace=0.25)

    # A: full light curve
    ax = fig.add_subplot(gs[0, :])
    tb, fb, _ = bin_lightcurve(lc.time, lc.raw, 30 / 1440)
    ax.plot(tb, fb, ",", color="0.4", rasterized=True)
    tt, tr, _ = bin_lightcurve(lc.time, lc.trend, 30 / 1440)
    ax.plot(tt, tr, ",", color="tab:orange", rasterized=True)
    n0 = np.floor((lc.time.min() - t0) / P)
    n1 = np.ceil((lc.time.max() - t0) / P)
    tc = t0 + np.arange(n0, n1 + 1) * P
    tc = tc[np.min(np.abs(tc[:, None] - lc.time[None, ::50]), axis=1) < 0.1] if len(tc) < 5000 else tc
    lo, hi = np.percentile(fb, [0.5, 99.5])
    ax.vlines(tc, lo, lo + 0.15 * (hi - lo), color="tab:red", lw=0.6)
    ax.set_ylim(lo - 0.1 * (hi - lo), hi + 0.1 * (hi - lo))
    ax.set_xlabel("time (BTJD = BJD - 2457000)")
    ax.set_ylabel("PDCSAP flux (30-min bins)")
    ax.set_title(f"A. Full light curve, {len(np.unique(lc.sector))} sectors; red ticks = predicted transits; "
                 "orange = detrending curve", fontsize=10, loc="left")

    # B: periodogram
    ax = fig.add_subplot(gs[1, 0])
    if periods is not None and sde_curve is not None:
        ax.plot(periods, sde_curve, lw=0.4, color="k")
        ax.axvline(P, color="tab:red", alpha=0.4, lw=3)
        ax.set_xscale("log")
        ax.set_xlabel("period (days)")
        ax.set_ylabel("SDE")
    ax.set_title("B. Stacked periodogram", fontsize=10, loc="left")

    # C: folded transit
    ax = fig.add_subplot(gs[1, 1])
    w = 3 * dur
    m = np.abs(dt) < w
    ax.plot(dt[m] * 24, (lc.flux[m] - 1) * 1e6, ".", ms=1, color="0.75", rasterized=True)
    xb, yb, eb = binned(dt[m] * 24, (lc.flux[m] - 1) * 1e6, dur * 24 / 8)
    ax.errorbar(xb, yb, eb, fmt="o", ms=3, color="k")
    if trap:
        xx = np.linspace(-w, w, 500)
        ax.plot(xx * 24, (trapezoid(xx, trap["depth"], trap["t14"], trap["ingress_frac"], trap["tc_offset"]) - 1) * 1e6,
                color="tab:red")
    depth_ppm = vet["depth_ppm"]
    ax.set_ylim(-2.5 * max(depth_ppm, 50), 1.5 * max(depth_ppm, 50))
    ax.set_xlabel("hours from transit centre")
    ax.set_ylabel("ppm")
    ax.set_title(f"C. Folded at P = {P:.5f} d", fontsize=10, loc="left")

    # D: odd vs even
    ax = fig.add_subplot(gs[1, 2])
    ep = np.round((lc.time - t0) / P).astype(int)
    for parity, off, col in ((1, -w, "tab:blue"), (0, w, "tab:green")):
        mm = m & (ep % 2 == parity)
        if mm.sum() > 10:
            xb, yb, eb = binned(dt[mm] * 24, (lc.flux[mm] - 1) * 1e6, dur * 24 / 6)
            ax.errorbar(xb + off * 24 * 1.1, yb, eb, fmt="o", ms=3, color=col)
    ax.axhline(-vet.get("depth_odd_ppm", np.nan), xmax=0.5, color="tab:blue", ls="--")
    ax.axhline(-vet.get("depth_even_ppm", np.nan), xmin=0.5, color="tab:green", ls="--")
    ax.set_ylim(-2.5 * max(depth_ppm, 50), 1.5 * max(depth_ppm, 50))
    ax.set_xticks([])
    ax.set_title(f"D. Odd (blue) vs even (green): {vet.get('oddeven_sigma', np.nan):.1f} sigma apart",
                 fontsize=10, loc="left")

    # E: full phase
    ax = fig.add_subplot(gs[2, 0:2])
    ph = (dt / P + 0.25) % 1 - 0.25
    xb, yb, eb = binned(ph, (lc.flux - 1) * 1e6, max(dur / P / 3, 1e-4))
    ax.errorbar(xb, yb, eb, fmt=".", ms=2, color="k", elinewidth=0.5)
    ax.axvline(0.5, color="tab:purple", alpha=0.3)
    if np.isfinite(vet.get("sec_phase", np.nan)):
        ax.axvline((vet["sec_phase"] + 0.25) % 1 - 0.25, color="tab:orange", ls=":", label="strongest other dip")
        ax.legend(fontsize=8, loc="lower right")
    ax.set_ylim(-2.5 * max(depth_ppm, 50), 1.5 * max(depth_ppm, 50))
    ax.set_xlabel("orbital phase (0 = transit)")
    ax.set_ylabel("ppm")
    ax.set_title("E. Whole orbit (secondary eclipse would appear near 0.5 for a circular binary)",
                 fontsize=10, loc="left")

    # F: per-transit depths
    ax = fig.add_subplot(gs[2, 2])
    pt = vet.get("per_transit")
    if pt:
        ax.errorbar(pt["time"], np.array(pt["depth"]) * 1e6, np.array(pt["err"]) * 1e6, fmt="o", ms=2,
                    color="k", elinewidth=0.5)
        ax.axhline(depth_ppm, color="tab:red")
        ax.axhline(0, color="0.5", lw=0.5)
    ax.set_xlabel("time (BTJD)")
    ax.set_ylabel("depth (ppm)")
    ax.set_title(f"F. Individual transits ({vet.get('n_transits_measured', 0)})", fontsize=10, loc="left")

    # G: text
    ax = fig.add_subplot(gs[3, :])
    ax.axis("off")
    meta = lc.star_meta or {}
    lines = [
        f"Verdict: {verdict.upper()}" + (f"  ({'; '.join(reasons)})" if reasons else ""),
        f"P = {P:.6f} d   T0 = {t0:.5f} BTJD   depth = {depth_ppm:.0f} ppm   duration = {dur * 24:.2f} h "
        f"(expected {vet.get('expected_duration_h', np.nan):.2f} h for b=0)   Rp ~ {vet.get('rp_rearth', np.nan):.2f} R_earth",
        f"SNR = {vet.get('snr_red', np.nan):.1f} (white-noise BLS {det['snr']:.1f})   folded red-noise SNR = "
        f"{vet.get('red_snr', np.nan):.1f}   SDE = {det['sde']:.1f}   "
        f"transits = {vet.get('n_transits_measured', 0)}   max single-transit share = {vet.get('max_single_frac', np.nan):.2f}   "
        f"depth chi2/dof = {vet.get('chi2_depth', np.nan):.2f}",
        f"secondary: strongest other dip {vet.get('sec_snr', np.nan):.1f} sigma at phase {vet.get('sec_phase', np.nan):.2f}; "
        f"phase-0.5 depth {vet.get('phase05_depth_ppm', np.nan):.0f} ppm ({vet.get('phase05_sigma', np.nan):.1f} sigma)   "
        f"uniqueness {vet.get('uniqueness', np.nan):.1f}",
        f"centroid: chi2 z={vet.get('centroid_z', np.nan):.1f}, shift z={vet.get('centroid_shift_z', np.nan):.1f}, "
        f"implied source offset {vet.get('centroid_offset_arcsec', np.nan):.0f}\"   "
        f"background rise z={vet.get('bkg_signed_z', np.nan):.1f}   red-noise beta={vet.get('red_noise_beta', np.nan):.2f}   "
        f"SAP depth {vet.get('sap_depth_ppm', np.nan):.0f} +- {vet.get('sap_depth_err_ppm', np.nan):.0f} ppm   "
        f"edge fraction {vet.get('edge_frac', np.nan):.2f}",
        f"star: Prot ~ {meta.get('prot', np.nan):.2f} d (LS power {meta.get('prot_power', 0):.2f}), "
        f"noise {meta.get('noise_ppm', np.nan):.0f} ppm per 2 min, detrend window {meta.get('window', np.nan):.2f} d, "
        f"flares removed {100 * meta.get('flare_fraction', 0):.2f}%",
        "catalogue matches: " + (
            ", ".join(f"{h['name']} ({h['relation']})" for h in matches[:4])
            + (f" + {len(matches) - 4} more" if len(matches) > 4 else "") if matches else "none"),
    ]
    lines = [l if len(l) <= 190 else l[:187] + "..." for l in lines]
    ax.text(0, 1, "\n".join(lines), va="top", family="monospace", fontsize=8.5)

    fig.suptitle(f"TIC {lc.tic}  {title_extra}", fontsize=14, x=0.01, ha="left")
    if path:
        fig.savefig(path, dpi=90, bbox_inches="tight")
    plt.close(fig)
