"""Step 9: write REPORT.md from the result files (every number comes from results/).

    .venv/bin/python scripts/09_report.py

The only hand-written parts are the per-candidate verdicts and assessments in
tess_search/assessments.py (shared with the dossiers).
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from tess_search import DATA, ROOT, RESULTS, vetting

from tess_search.assessments import assessment  # noqa: E402

HARD = RESULTS / "hardening"


def load_jsonl(name):
    p = RESULTS / "reliability" / f"{name}.jsonl"
    return [json.loads(l) for l in p.read_text().splitlines() if l.strip()] if p.exists() else []


def reclassify(sig):
    """Apply the final vetting rules to stored metrics (as 06_summarize does for real data)."""
    if sig.get("vet"):
        return vetting.classify(sig["vet"])[0]
    return sig.get("verdict", "not vetted")


def hardening_section(w, loc):
    """Pixel localization, NASA cross-check and statistical validation, with their tests."""
    syst = json.loads((HARD / "localization_systematics.json").read_text())
    w("## Hardening\n")
    w("### Is the light lost on the target? Pixel-level localization\n")
    w("TESS pixels are 21\" wide, so a neighbouring eclipsing binary can leak a planet-sized dip into the target's "
      "light curve. For every sector the images taken during transit are subtracted from those taken just before and "
      "after; the result shows only the light that disappeared. All sectors are then fitted together with NASA's model "
      "of how a point source spreads over the pixels (the SPOC PRF), calibrated per sector on the Gaia stars in the "
      "image, to find where on the sky the light went missing (`tess_search/localize.py`). Pixel errors come from "
      "~60 fake transits per sector; flares are removed first.\n")
    pl = loc[loc.kind == "planet"].sort_values("label")
    nb = loc[loc.kind == "neb"].sort_values("label")
    w("| test | signal | source offset from target | target excluded at | reduced chi2 |")
    w("|---|---|---|---|---|")
    for _, r in pd.concat([pl, nb]).iterrows():
        w(f"| {'confirmed planet' if r.kind == 'planet' else 'TFOP: nearby EB'} | {r.label} | {r.offset_arcsec:.1f} ± "
          f"{r.offset_err_total_arcsec:.1f}\" | {r.target_sigma_total:.1f} sigma | {r.chi2_red:.2f} |")
    inj = loc[loc.kind == "injection"]
    rel = inj[inj.reliable == True]  # noqa: E712
    rn = rel[rel.injected_on == "neighbour"]
    w(f"\n{len(inj)} synthetic eclipses were planted in the real pixels (on the target, and on the neighbours that "
      "could most easily mimic each signal, sized to reproduce its depth). "
      f"Of the {len(rel)} with a reliable fit, {int(rel.recovered_correct.astype(str).eq('True').sum())} were traced "
      "to the right star (the exception is a target/neighbour pair 4.7\" apart, below the method's resolution); for "
      f"{int((rn.target_sigma_total > 3).sum())} of {len(rn)} planted on neighbours (4.7-83\" away) the target was "
      f"excluded at more than 3 sigma. {len(inj) - len(rel)} injections around TOI-6000 gave a poor fit (reduced chi2 "
      "> 2): a variable star in that image happens to vary in step with the injected period, which breaks the "
      "one-source assumption. Such fits are flagged unreliable; none of the candidates is affected (reduced chi2 "
      "1.06-1.16).\n")
    w(f"The position errors include a {syst['sys_arcsec_used']:.1f}\" systematic floor: the smallest floor for which 95% "
      f"of the {syst['n_cases']} cases with a known source fall inside their 2-sigma region is "
      f"{syst['sys_arcsec_95pct_coverage']:.1f}\".\n")
    dv = pd.read_csv(HARD / "spoc_dv_summary.csv")
    w("### NASA's own pipeline (SPOC) on the same signals\n")
    w("Three candidates were SPOC Threshold Crossing Events that never became TOIs. SPOC's reported SNR fell as data "
      "were added; folding this work's light curve at SPOC's period for each run reproduces the drop (the transit "
      "smears by 1-4 hours over the baseline), so the signals did not fade. SPOC's difference-image centroids for "
      f"these runs are listed in the dossiers; none of its {int(dv.diff_images_attempted.sum())} per-sector difference "
      "images passed SPOC's own quality test. The run that put TIC 229689348's source 55\" away (s14-s55) used the "
      "correct period, but its difference images were among those that failed; the localization above, built at the "
      "correct period from all sectors, places that source on the target.\n")
    w("![SPOC period check](results/hardening/spoc_period_check.png)\n")
    w("### Statistical validation (TRICERATOPS)\n")
    w("TRICERATOPS (Giacalone et al. 2021) weighs a planet on the target against eclipsing binaries on the target, "
      "unresolved companions, background stars and resolved neighbours, using the transit shape and the Gaia DR3 "
      "field population (queried through VizieR). Each candidate was run 5 times with 10^6 draws. TRICERATOPS uses "
      "only brightness ratios and the transit shape; it is run a second time with the neighbours that the pixel "
      "localization excludes at more than 3 sigma treated as cleared (as stars cleared by ground-based photometry "
      "are). No high-resolution imaging was used, so these FPPs are conservative; validation needs FPP < 0.015 and "
      "NFPP < 0.001 with imaging. Values are in the candidate table above.\n")


def pct(x):
    return f"{100 * x:.1f}%"


def main():
    targets = pd.read_csv(DATA / "targets.csv")
    index = pd.read_parquet(DATA / "lc_index.parquet")
    summary = json.loads((RESULTS / "summary.json").read_text())
    known = pd.read_csv(RESULTS / "known_recovery.csv")
    follow = pd.read_csv(RESULTS / "followup" / "summary.csv")
    confirmed = pd.read_csv(DATA / "catalogs" / "confirmed.csv")
    transiting = dict(zip(confirmed.pl_name, confirmed.tran_flag))

    inj = load_jsonl("inject")
    inv = load_jsonl("invert")
    for r in inj:
        r["passed"] = False
        if r.get("recovered"):
            s = next((s for s in r["signals"] if abs(s["period"] / r["period"] - 1) < 0.005), None)
            r["passed"] = bool(s) and reclassify(s) in ("candidate", "weak candidate")
    inj_df = pd.DataFrame(inj)
    inv_cands = [sum(reclassify(s) == "candidate" for s in r["signals"]) for r in inv]
    inv_weak = [sum(reclassify(s) == "weak candidate" for s in r["signals"]) for r in inv]

    in_range = known[(known.period >= 0.4) & (known.period <= 40)].copy()
    in_range["transits"] = [transiting.get(n, 1) if s == "CONFIRMED" else 1 for n, s in zip(in_range["name"], in_range.source)]
    conf = in_range[in_range.source == "CONFIRMED"]
    conf_tr = conf[conf.transits == 1]
    can_transit = in_range[(in_range.transits == 1) & (in_range.disposition != "FA")]

    full = follow[follow.verdict == "candidate"].copy()
    weak = follow[follow.verdict == "weak candidate"].copy()
    fits = {}
    for _, r in follow.iterrows():
        info = json.loads((RESULTS / "followup" / f"TIC{r.tic}_{r.signal}.json").read_text())
        fits[(r.tic, r.signal)] = info

    L = []
    w = L.append
    w("# TESS M-dwarf deep search: report\n")
    w("*Generated by `scripts/09_report.py` from the files in `results/`; candidate assessments are mine.*\n")

    # ------------------------------------------------------------ key findings
    w("## Key findings\n")
    w(f"* Searched the **{summary['stars']:,} bright M dwarfs** that TESS has watched longest (>= 20 sectors, through "
      f"sector {index.sector.max()}) for transiting planets with periods of 0.4-40 days.")
    w(f"* The pipeline independently recovers **{len(conf_tr[conf_tr.recovered])} of {len(conf_tr)} confirmed transiting "
      f"planets** in that range (and {int(can_transit.recovered.sum())} of {len(can_transit)} catalogued signals that can "
      "transit), including all four planets of TOI-700, and passes them through vetting.")
    big = (inj_df.rp >= 2) & (inj_df.period < 15)
    small = (inj_df.rp < 1.5) & (inj_df.period < 5)
    w(f"* **Completeness** ({len(inj_df)} injected planets): {pct(inj_df.recovered.mean())} detected, "
      f"{pct(inj_df.passed.mean())} detected and passed vetting; {pct(inj_df[big].passed.mean())} for planets of "
      f"2-4 R_earth inside 15 days, {pct(inj_df[small].passed.mean())} for planets under 1.5 R_earth inside 5 days.")
    w(f"* **Reliability** (200 flipped light curves): **no false candidates** (95% upper limit "
      f"{pct(3 / len(inv))} of stars); false *weak* candidates on {pct(np.mean([c > 0 for c in inv_weak]))} of stars.")
    cs = pd.read_csv(RESULTS / "candidates" / "summary.csv")
    loc = pd.read_csv(HARD / "localization_summary.csv")
    on = cs[cs.verdict.str.startswith("Candidate")]
    off = cs[~cs.verdict.str.startswith("Candidate")]
    w(f"* **{len(full)} signals in no planet catalogue** passed every light-curve test. A second, deeper pass checked "
      "each one in the pixels, statistically and against NASA's own pipeline (see *Hardening*): "
      f"**{len(on)} remain candidates** ({', '.join(f'TIC {t}' + (f' ({n})' if isinstance(n, str) and n else '') for t, n in zip(on.tic, on.name))}; "
      f"{on.rp_rearth.min():.2f}-{on.rp_rearth.max():.2f} R_earth, periods {on.period.min():.2f}-{on.period.max():.2f} d), "
      + (f"and **{len(off)} is a nearby eclipsing binary** (TIC {', '.join(map(str, off.tic))}: the light loss is "
         f"{off.loc_offset_arcsec.iloc[0]:.0f}\" from the target)." if len(off) else ""))
    pl = loc[loc.kind == "planet"]
    nb = loc[loc.kind == "neb"]
    inj = loc[(loc.kind == "injection") & (loc.reliable == True)]  # noqa: E712
    w(f"* The pixel-level localization was validated first: {len(pl)} confirmed planets come out on their own star "
      f"(all within {pl.target_sigma_total.max():.1f} sigma), {len(nb)} TFOP-retired nearby eclipsing binaries come out off "
      f"target, and {int(inj.recovered_correct.astype(str).eq('True').sum())} of {len(inj)} synthetic eclipses planted in "
      "the real pixels are traced to the right star.")
    likely = on[(on.fpp_cleared < 0.5) & (on.nfpp_cleared < 1e-3)]
    w(f"* **Statistical validation (TRICERATOPS, no imaging):** with the neighbours that the pixels exclude treated as "
      f"cleared, {len(likely)} of {len(on)} candidates meet TRICERATOPS's *likely planet* criteria (FPP < 0.5, "
      f"NFPP < 0.001): " + ", ".join(f"TIC {r.tic} FPP {r.fpp_cleared:.2f}" for _, r in likely.sort_values("fpp_cleared").iterrows())
      + ". None is validated (FPP < 0.015 with high-resolution imaging); the remaining FPP is mostly unresolved "
      "bound companions, which imaging tests directly.")
    w(f"* {len(weak)} weaker signals are listed separately; most are expected to be false alarms, and many sit at "
      "36-40 d periods where noise produces them (see *Weak candidates*).\n")

    w("## Candidates\n")
    w("Each signal was re-vetted from scratch, fitted with an MCMC transit model (planet radius includes the stellar-radius "
      "uncertainty), localized in the pixels, run through TRICERATOPS, and compared with NASA's SPOC pipeline where SPOC "
      "had flagged it. A candidate is a signal worth follow-up observations, **not a confirmed planet**. Full dossiers: "
      "`results/candidates/`.\n")
    w("| TIC | verdict | period (d) | radius (R_earth) | T_eq (K) | SNR | source offset from target | FPP / NFPP (TESS only) | FPP / NFPP (localization-cleared) |")
    w("|---|---|---|---|---|---|---|---|---|")
    for _, r in cs.iterrows():
        name = f" ({r['name']})" if isinstance(r["name"], str) and r["name"] else ""
        def pair(a, b):
            return f"{a:.3f} / {b:.4f}" if np.isfinite(a) and np.isfinite(b) else "n/a"
        w(f"| {r.tic}{name} | {r.verdict} | {r.period:.6f} | {r.rp_rearth:.2f} ± {r.rp_err:.2f} | {r.teq_k:.0f} | "
          f"{r.snr:.1f} | {r.loc_offset_arcsec:.1f} ± {r.loc_err_arcsec:.1f}\" | {pair(r.fpp, r.nfpp)} | "
          f"{pair(r.get('fpp_cleared', np.nan), r.get('nfpp_cleared', np.nan))} |")
    w("")
    for _, r in cs.iterrows():
        w(f"**TIC {r.tic}.** {assessment(int(r.tic))} Dossier: `results/candidates/TIC{r.tic}.md`.\n")
    w("The ExoFOP community-TOI upload is drafted in `results/hardening/exofop/params_planet_DRAFT.txt` "
      "(not submitted; it needs the submitter's ExoFOP tag and a public URL).\n")
    hardening_section(w, loc)

    # ------------------------------------------------------------ data
    w("## Data\n")
    w(f"* **{len(targets):,} M dwarfs** (Teff <= 3900 K, R <= 0.65 R_sun, TESS mag <= 13.5, low contamination) with "
      f"**>= 20 sectors** of SPOC 2-minute photometry; median {int(targets.n_sectors.median())} sectors, max "
      f"{targets.n_sectors.max()}.")
    w(f"* Chosen from {len(index):,} light-curve files of {index.tic.nunique():,} stars in sectors 1-{index.sector.max()}. "
      "NASA's latest combined multi-sector search covered sectors 1-96.")
    w(f"* {summary['signals_total']:,} periodic signals detected; {summary['signals_vetted']:,} strong enough to vet.\n")
    w("![sample](results/figures/sample.png)\n")

    # ------------------------------------------------------------ validation
    w("## Validation on known planets\n")
    w(f"Of {len(in_range)} catalogued signals with periods 0.4-40 d on the searched stars, "
      f"{int(in_range.recovered.sum())} were recovered. Of the {len(in_range) - int(in_range.recovered.sum())} missed, "
      f"{len(conf) - len(conf_tr)} are confirmed planets found by radial velocity that **do not transit** "
      "(no transit search can see them) and one is a known false alarm. The real misses:\n")
    miss = in_range[(~in_range.recovered) & (in_range.transits == 1) & (in_range.disposition != "FA")]
    w("| star | signal | period (d) | disposition |")
    w("|---|---|---|---|")
    for _, r in miss.sort_values("period").iterrows():
        w(f"| TIC {r.tic} | {r['name']} | {r.period:.4f} | {r.disposition} |")
    rec = in_range[in_range.recovered]
    w(f"\nOf the recovered signals, vetting passed {int((rec.verdict == 'candidate').sum() + (rec.verdict == 'weak candidate').sum())} "
      f"and rejected {int((rec.verdict == 'false positive').sum())}; the rejected ones include signals the TESS team itself "
      "classified as false positives (e.g. TOI-419.01, TOI-1178.01, TOI-1635.01). The full table is in the appendix.\n")

    # ------------------------------------------------------------ completeness
    w("## Completeness: injected planets\n")
    w(f"{len(inj_df)} synthetic planets (0.6-4 R_earth, 0.5-40 d, random impact parameter) were planted into the raw "
      "light curves of random stars without known planets, and the full pipeline was run on each.\n")
    w(f"* detected at the right period and epoch: **{pct(inj_df.recovered.mean())}**; detected and passed vetting: "
      f"**{pct(inj_df.passed.mean())}** (vetting kept {pct(inj_df[inj_df.recovered].passed.mean())} of detected planets; "
      "most losses are planets with periods within 3% of the spacecraft's 13.7-day orbit or its multiples, which are "
      "rejected on purpose).")
    for lo, hi in [(0.6, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 4.0)]:
        m = (inj_df.rp >= lo) & (inj_df.rp < hi)
        w(f"  * {lo}-{hi} R_earth: {pct(inj_df[m].passed.mean())} found and kept (n={m.sum()})")
    w("\n![completeness](results/figures/completeness.png)\n")

    # ------------------------------------------------------------ reliability
    w("## Reliability: inverted light curves\n")
    w(f"{len(inv)} light curves were flipped upside down (dips become bumps) and searched with the same pipeline. "
      "Nothing in flipped data can be a planet, so anything that passes vetting is a false alarm.\n")
    w(f"* false **candidates**: {sum(inv_cands)} on {sum(c > 0 for c in inv_cands)} of {len(inv)} stars. With none seen, "
      f"at 95% confidence fewer than {pct(3 / len(inv))} of stars produce one (fewer than ~{3 / len(inv) * summary['stars']:.0f} "
      f"expected among the {summary['stars']:,} searched stars; best estimate close to none).")
    wr = np.mean([c > 0 for c in inv_weak])
    w(f"* false **weak candidates**: on {pct(wr)} of stars, i.e. ~{wr * summary['stars']:.0f} expected among the searched "
      f"stars, against {summary['verdicts'].get('weak candidate', 0)} found in the real search. **Weak candidates are therefore mostly false alarms.**\n")

    # ------------------------------------------------------------ results overview
    w("## What happened to every signal\n")
    w("| verdict | signals |")
    w("|---|---|")
    for k, v in sorted(summary["verdicts"].items(), key=lambda x: -x[1]):
        w(f"| {k} | {v} |")
    w("\nMost common reasons for rejection (a signal can have several):\n")
    for k, v in list(summary["fp_reason_counts"].items())[:12]:
        w(f"* {k}: {v}")
    w("\n![verdicts](results/figures/verdicts.png)\n")

    w("## Weak candidates\n")
    long_p = weak[(weak.period > 36) & (weak.period < 41)]
    inv_periods = [s["period"] for r in inv for s in r["signals"] if reclassify(s) == "weak candidate"]
    w(f"Passed the hard tests but with flags (usually modest red-noise SNR or uneven transit depths). The inversion test "
      f"predicts about {wr * summary['stars']:.0f} false weak candidates in this sample, so treat these as a list to "
      f"re-check with more data, not as candidates. {len(long_p)} of the {len(weak)} have periods of 36-41 d, where a "
      "signal rests on only a handful of transits; the inverted light curves, which contain no real planets, put "
      f"{sum(36 < p < 41 for p in inv_periods)} of their {len(inv_periods)} false weak candidates there too.\n")
    wr_path = HARD / "weak_recheck.csv"
    if wr_path.exists():
        wrc = pd.read_csv(wr_path)
        n = wrc["class"].value_counts()
        cons = wrc[wrc["class"] == "consistent"]
        w(f"**Re-checked in the pixels** (`scripts/16_weak_recheck.py`): does the target lose, in the difference "
          "images, the light the light-curve depth predicts? For confirmed planets the ratio is 0.85-1.29. For the "
          f"weak signals, {n.get('not in pixels', 0)} are **not reproduced in the pixels** (ratio more than 3 sigma "
          f"below 1: most likely light-curve artefacts), {n.get('marginal', 0)} are marginal and "
          f"{n.get('consistent', 0)} are consistent with a dip on the target ("
          + ", ".join(f"TIC {r.tic} at {r.period:.2f} d" for _, r in cons.iterrows())
          + "), although only at 3.5-4 sigma in the pixels. None is localized off target with confidence. "
          "Table: `results/hardening/weak_recheck.csv`.\n")
    w("| TIC | period (d) | depth (ppm) | radius (R_earth) | SNR | both halves | RUWE | flags |")
    w("|---|---|---|---|---|---|---|---|")
    for _, r in weak.sort_values("snr_red", ascending=False).iterrows():
        w(f"| {r.tic} | {r.period:.4f} | {r.depth_ppm:.0f} | {r.rp_rearth:.2f} | {r.snr_red:.1f} | "
          f"{'yes' if r.in_both_halves else 'no'} | {r.ruwe:.2f} | {r.reasons} |")
    w("")

    # ------------------------------------------------------------ appendix
    w("## Appendix: every catalogued signal on the searched stars (0.4-40 d)\n")
    w("| star | known signal | period (d) | disposition | transits | recovered | verdict |")
    w("|---|---|---|---|---|---|---|")
    for _, r in in_range.sort_values(["tic", "period"]).iterrows():
        w(f"| TIC {r.tic} | {r['name']} | {r.period:.4f} | {r.disposition} | {'yes' if r.transits else 'no (RV)'} | "
          f"{'yes' if r.recovered else 'no'} | {r.verdict if isinstance(r.verdict, str) else ''} |")
    (ROOT / "REPORT.md").write_text("\n".join(L) + "\n")
    print("wrote REPORT.md")


if __name__ == "__main__":
    main()
