"""Step 15: final candidate dossiers, summary table and a draft ExoFOP CTOI upload file.

Pulls together, for each of the five candidates:
  results/followup/          detection, vetting, Gaia neighbours, half-data searches
  results/hardening/localize pixel-level source position (scripts/12_localize.py)
  results/hardening/mcmc     transit parameters with full uncertainties (scripts/14_mcmc.py)
  results/hardening/triceratops  false-positive probabilities (scripts/13_*.py)
  results/hardening/spoc_*   SPOC DV cross-check and period-drift test (scripts/11_*.py)

Writes results/candidates/TIC<tic>.md, results/candidates/summary.csv and
results/hardening/exofop/params_planet_DRAFT.txt (NOT submitted: the tag and paper
fields are placeholders that need the submitter's ExoFOP username and a public URL).
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from tess_search import DATA, RESULTS

OUT = RESULTS / "candidates"
HARD = RESULTS / "hardening"
CANDIDATES = [(32090583, 5, "TOI-218"), (229689348, 1, None), (198412174, 1, None), (149390648, 1, None),
              (294053492, 1, None)]
PREDICT_BJD = 2461557.5   # 2027 June 1: timing uncertainty quoted at this date
from tess_search.assessments import SUBMIT, VERDICT, assessment  # noqa: E402  hand-written verdicts


def load(path):
    return json.loads(Path(path).read_text()) if Path(path).exists() else None


def q(d, fmt="{:.2f}"):
    """median +plus -minus from an mcmc summary entry."""
    if d is None:
        return "n/a"
    m, p, mi = d["median"], d["plus"], d["minus"]
    if abs(p - mi) / max(p, mi, 1e-12) < 0.15:
        return f"{fmt.format(m)} ± {fmt.format(0.5 * (p + mi))}"
    return f"{fmt.format(m)} +{fmt.format(p)} −{fmt.format(mi)}"


def ephemeris(fu, mc):
    """Period from the least-squares fit (errors inflated by the red-noise factor), epoch from MCMC."""
    f = fu["fit"]
    beta = mc["prior"]["red_noise_beta"]
    p, p_err = f["period"], f["period_err"] * max(beta, 1.0)
    t0 = mc["prior"]["t0_btjd"]["median"] + 2457000.0
    t0_err = 0.5 * (mc["prior"]["t0_btjd"]["plus"] + mc["prior"]["t0_btjd"]["minus"])
    n = round((PREDICT_BJD - t0) / p)
    sig_min = float(np.hypot(t0_err, n * p_err) * 1440)
    return p, p_err, t0, t0_err, sig_min


def fmt_err(v, e, digits=None):
    if digits is None:
        digits = int(np.clip(1 - np.floor(np.log10(e)), 0, 8)) if e > 0 else 6
    return f"{v:.{digits}f} ± {e:.{digits}f}"


def dossier(tic, sig, name, fu, loc, mc, tri, dv, pchk, star, tri_c=None):
    p, p_err, t0, t0_err, sig_min = ephemeris(fu, mc)
    s = mc["prior"]
    fr = mc["free_density"]
    v = fu["vet"]
    g = fu.get("gaia", {})
    title = f"TIC {tic}" + (f" ({name})" if name else "")
    L = [f"# {title}: {VERDICT.get(tic, '')}", ""]
    text = assessment(tic)
    if text:
        L += [text, ""]
    L += ["## Ephemeris and transit parameters",
          "MCMC fit (emcee, batman; Rp/R*, b, stellar density with TIC prior, quadratic limb darkening with priors; "
          "1-min folded bins with red-noise-scaled errors). Planet radius includes the TIC stellar-radius error.", "",
          "| parameter | value |", "|---|---|",
          f"| period (d) | {fmt_err(p, p_err)} |",
          f"| mid-transit (BJD_TDB) | {fmt_err(t0, t0_err, 5)} |",
          f"| timing uncertainty on 2027 June 1 | ±{sig_min:.0f} min |",
          f"| depth (ppm) | {q(s['depth_ppm'], '{:.0f}')} |",
          f"| duration T14 (h) | {q(s['t14_h'])} |",
          f"| Rp/R* | {q(s['rp_rs'], '{:.4f}')} |",
          f"| impact parameter b | {q(s['b'])} |",
          f"| a/R* | {q(s['a_rs'], '{:.1f}')} |",
          f"| inclination (deg) | {q(s['inc_deg'], '{:.1f}')} |",
          f"| **planet radius (R_earth)** | **{q(s['rp_rearth'])}** |",
          f"| semi-major axis (AU) | {q(s['a_au'], '{:.4f}')} |",
          f"| equilibrium temperature (K, zero albedo) | {q(s['teq_k'], '{:.0f}')} |",
          f"| insolation (Earth = 1) | {q(s['insolation_earth'], '{:.0f}')} |",
          f"| stellar density, transit shape alone / TIC | {q(fr['rho_ratio_to_tic'])} |",
          ""]
    L += ["## Star", f"TESS mag {star['Tmag']:.2f}; Teff {star['Teff']:.0f} ± {star['e_Teff']:.0f} K; "
          f"R* {star['rad']:.3f} ± {star['e_rad']:.3f} R_sun; M* {star['mass']:.3f} M_sun (TIC v8.2); "
          f"distance {star['d']:.1f} pc; Gaia DR3 {star['GAIA']} (RUWE {g.get('target_ruwe', float('nan')):.2f}). "
          f"{len(fu['sectors'])} sectors of SPOC 2-min data.", ""]
    hs = [h for h in fu["halves"] if h.get("snr") is not None]
    L += ["## Detection and light-curve vetting",
          f"Red-noise-aware SNR {v.get('snr_red', float('nan')):.1f}; {v.get('n_transits_measured')} transits with data; "
          f"odd/even depths differ by {v.get('oddeven_sigma', float('nan')):.1f} sigma; depth at phase 0.5: "
          f"{v.get('phase05_depth_ppm', float('nan')):.0f} ppm ({v.get('phase05_sigma', float('nan')):.1f} sigma); "
          f"largest single-transit share {v.get('max_single_frac', float('nan')):.2f}; depth in SAP flux "
          f"{v.get('sap_depth_ppm', float('nan')):.0f} ppm vs {v.get('depth_ppm', float('nan')):.0f} ppm in PDCSAP. "
          "Independent searches of each half of the data: "
          + ", ".join(f"{h['half']} SNR {h['snr']:.1f}" + (" at the same period" if h.get("period_match") else " at a different period")
                      for h in hs) + ".", ""]
    if dv is not None and len(dv):
        L += ["## NASA SPOC pipeline history",
              "| SPOC multi-sector run | SPOC period (d) | SPOC SNR | this light curve's SNR at SPOC's period | transit drift over the baseline (h) | SPOC difference-image offset |",
              "|---|---|---|---|---|---|"]
        runs = {r["run"]: r for r in (pchk or {}).get("runs", [])}
        for _, r in dv.iterrows():
            pr = runs.get(r.run, {})
            L.append(f"| {r.run} | {r.period:.6f} | {r.fit_snr:.1f} | {pr.get('our_snr_at_spoc_period', float('nan')):.1f} | "
                     f"{pr.get('drift_over_baseline_h', float('nan')):.1f} | {r.ms_tic_offset_arcsec:.1f} ± "
                     f"{r.ms_tic_offset_err_arcsec:.1f}\" ({int(r.diff_images_good)}/{int(r.diff_images_attempted)} "
                     "difference images passed SPOC's quality test) |")
        if pchk:
            L += ["", f"At this work's period ({pchk['our_period']:.6f} d) the same light curve gives box SNR "
                  f"{pchk['our_snr_box']:.1f}. SPOC's SNR changes between runs follow its period estimate: a period "
                  "error of a few 1e-5 d smears a sub-hour transit by hours over the ~2,000-d baseline."]
        L += [""]
    if loc:
        alive = [r for r in loc["stars"] if not r["is_target"] and r["excluded_sigma_total"] < 3]
        comp = sorted([r for r in loc["stars"] if not r["is_target"]], key=lambda r: r["sep_arcsec"])[:3]
        L += ["## Pixel-level localization (where the light goes missing)",
              f"Difference images from {loc['n_sectors']} sectors ({loc['n_transits']} transits) fitted jointly with the "
              "SPOC PRF (scripts/12_localize.py). Errors include a 1.5\" systematic floor measured on 46 cases with "
              "known sources.", "",
              f"* best-fitting source position: {loc['offset_arcsec']:.1f} ± {loc['offset_err_total_arcsec']:.1f}\" from the "
              f"target (E {loc['offset_east']:+.0f}\", N {loc['offset_north']:+.0f}\"); the target is "
              f"{loc['target_sigma_total']:.1f} sigma from the best position",
              f"* light lost in transit: {loc['amplitude']:.2f} ± {loc['amplitude_err']:.2f} e-/s; "
              f"{loc['expected_amplitude']:.2f} e-/s expected if the target hosts the observed depth "
              f"(ratio {loc['amplitude_ratio']:.2f}; confirmed planets in the test set span 0.85-1.30)",
              "* nearest neighbours: " + "; ".join(
                  f"Gaia DR3 {r['source_id']} at {r['sep_arcsec']:.1f}\" (T = {r['tmag']:.1f}) excluded at "
                  f"{r['excluded_sigma_total']:.1f} sigma" for r in comp),
              "* neighbours not excluded at 3 sigma: " + ("none" if not alive else "; ".join(
                  f"Gaia DR3 {r['source_id']} at {r['sep_arcsec']:.1f}\" (T = {r['tmag']:.1f}; would need a "
                  f"{100 * r['implied_eclipse_depth']:.1f}% eclipse)" for r in alive))]
        if loc.get("split_halves"):
            L.append("* odd / even sectors separately: " + "; ".join(
                f"{h['half']} {np.hypot(h['offset_east'], h['offset_north']):.1f} ± {h['offset_err_total_arcsec']:.1f}\" "
                f"(target {h['target_sigma_total']:.1f} sigma)" for h in loc["split_halves"]))
        L += ["", f"![localization](../hardening/localize/TIC{tic}.png)", ""]
    if tri:
        top = ", ".join(f"{k.split(':')[0]} on TIC {k.split(':')[1]} {100 * v_:.1f}%" for k, v_ in tri["top_scenarios"][:4])
        L += ["## Statistical validation (TRICERATOPS)",
              f"FPP = {tri['FPP_mean']:.3f} ± {tri['FPP_std']:.3f}, NFPP = {tri['NFPP_mean']:.5f} ± {tri['NFPP_std']:.5f} "
              f"({tri['n_runs']} runs of {tri['n_draws']:,} draws; sectors {tri['sectors']}; Gaia DR3 field population). "
              f"Most probable scenarios: {top}. No high-resolution imaging was used, so unresolved companions are "
              "limited only by Gaia; imaging would lower the FPP.", ""]
        if tri_c:
            cl = tri_c.get("cleared_by_localization", [])
            L += [f"With the neighbours that the pixel localization excludes at more than 3 sigma treated as cleared "
                  f"({len(cl)} stars: " + (", ".join(f"TIC {c['tic']} at {c['sep']:.1f}\" ({c['excluded_sigma']:.1f} sigma)"
                                                  for c in cl[:4]) or "none") + (", ..." if len(cl) > 4 else "")
                  + f"), FPP = {tri_c['FPP_mean']:.3f} ± {tri_c['FPP_std']:.3f} and NFPP = {tri_c['NFPP_mean']:.5f} ± "
                  f"{tri_c['NFPP_std']:.5f}.", ""]
    L += ["## Files",
          f"* follow-up sheet and vetting: `results/followup/TIC{tic}_{sig}.png`, `.json`",
          f"* transit fit: `results/hardening/mcmc/TIC{tic}.png`, `.json`",
          f"* localization: `results/hardening/localize/TIC{tic}.png`, `.json`",
          f"* TRICERATOPS: `results/hardening/triceratops/TIC{tic}_result.json`", ""]
    return "\n".join(L), {"tic": tic, "name": name or "", "period": p, "period_err": p_err, "t0_bjd": t0, "t0_err": t0_err,
                          "timing_sigma_2027_min": sig_min}


def ctoi_row(tic, sig, row, mc, loc, tri, notes):
    s = mc["prior"]

    def pm(d):
        return d["median"], 0.5 * (d["plus"] + d["minus"])

    depth, depth_e = pm(s["depth_ppm"])
    dur, dur_e = pm(s["t14_h"])
    inc, inc_e = pm(s["inc_deg"])
    b, b_e = pm(s["b"])
    k, k_e = pm(s["rp_rs"])
    ar, ar_e = pm(s["a_rs"])
    rp, rp_e = pm(s["rp_rearth"])
    teq, teq_e = pm(s["teq_k"])
    ins, ins_e = pm(s["insolation_earth"])
    rho, rho_e = pm(s["rho_cgs"])
    a, a_e = pm(s["a_au"])
    f = [f"TIC{tic}.01", "newctoi", "PC", "TESS", "",
         f"{row['period']:.7f}", f"{row['period_err']:.7f}", f"{row['t0_bjd']:.5f}", f"{row['t0_err']:.5f}",
         f"{depth:.0f}", f"{depth_e:.0f}", f"{dur:.3f}", f"{dur_e:.3f}", f"{inc:.2f}", f"{inc_e:.2f}",
         f"{b:.3f}", f"{b_e:.3f}", f"{k:.4f}", f"{k_e:.4f}", f"{ar:.2f}", f"{ar_e:.2f}", f"{rp:.3f}", f"{rp_e:.3f}",
         "", "", f"{teq:.0f}", f"{teq_e:.0f}", f"{ins:.1f}", f"{ins_e:.1f}", f"{rho:.2f}", f"{rho_e:.2f}",
         f"{a:.5f}", f"{a_e:.5f}", "", "", "", "", "", "", "", "",
         "YYYYMMDD_EXOFOPUSERNAME_mdwarfdeepsearch_00001", "", "0", "REQUIRED_PUBLIC_URL", notes[:120]]
    return "|".join(f)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    tic_tab = pd.read_parquet(DATA / "tic_15plus.parquet").set_index("tic")
    dv_all = pd.read_csv(HARD / "spoc_dv_summary.csv") if (HARD / "spoc_dv_summary.csv").exists() else pd.DataFrame()
    pchk_all = load(HARD / "spoc_period_check.json") or {}
    rows, ctoi = [], []
    for tic, sig, name in CANDIDATES:
        fu = load(RESULTS / "followup" / f"TIC{tic}_{sig}.json")
        loc = load(HARD / "localize" / f"TIC{tic}.json")
        mc = load(HARD / "mcmc" / f"TIC{tic}.json")
        tri = load(HARD / "triceratops" / f"TIC{tic}_result.json")
        tri_c = load(HARD / "triceratops" / f"TIC{tic}_result_cleared.json")
        dv = dv_all[dv_all.tic == tic].sort_values("run") if len(dv_all) else None
        star = tic_tab.loc[tic].to_dict()
        text, row = dossier(tic, sig, name, fu, loc, mc, tri, dv, pchk_all.get(str(tic)), star, tri_c)
        (OUT / f"TIC{tic}.md").write_text(text)
        s = mc["prior"]
        row.update(verdict=VERDICT.get(tic, ""), rp_rearth=s["rp_rearth"]["median"],
                   rp_err=0.5 * (s["rp_rearth"]["plus"] + s["rp_rearth"]["minus"]),
                   depth_ppm=s["depth_ppm"]["median"], t14_h=s["t14_h"]["median"], b=s["b"]["median"],
                   teq_k=s["teq_k"]["median"], insolation=s["insolation_earth"]["median"],
                   snr=fu["vet"].get("snr_red"), loc_offset_arcsec=loc["offset_arcsec"] if loc else None,
                   loc_err_arcsec=loc["offset_err_total_arcsec"] if loc else None,
                   loc_target_sigma=loc["target_sigma_total"] if loc else None,
                   fpp=tri["FPP_mean"] if tri else None, nfpp=tri["NFPP_mean"] if tri else None,
                   fpp_cleared=tri_c["FPP_mean"] if tri_c else None, nfpp_cleared=tri_c["NFPP_mean"] if tri_c else None)
        rows.append(row)
        if SUBMIT.get(tic):
            ctoi.append(ctoi_row(tic, sig, row, mc, loc, tri, SUBMIT[tic] if isinstance(SUBMIT[tic], str) else ""))
        print(f"TIC {tic}: written", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False)
    template = (HARD / "exofop" / "exofop_template_params_planet.txt").read_text().splitlines()
    header = next(line for line in template if line.startswith("target|"))
    (HARD / "exofop" / "params_planet_DRAFT.txt").write_text(
        "\\ DRAFT - NOT SUBMITTED. Fill in the tag (your ExoFOP username) and paper (a public URL), check the next\n"
        "\\ free TIC<id>.NN candidate number on each target's ExoFOP overview page, then rename to\n"
        "\\ params_planet_YYYYMMDD_001.txt before uploading. Format: ExoFOP planet-parameter bulk upload template.\n"
        + header + "\n" + "\n".join(ctoi) + "\n")
    print(f"CTOI draft: {len(ctoi)} rows")
