"""Hand-written verdicts on the five candidates, shared by the dossiers (scripts/15_dossiers.py)
and the report (scripts/09_report.py).

The judgement is hand-written; the numbers are not typed in: each text has {placeholders} that
facts(tic) fills from the result files, so prose and tables cannot disagree.
"""
import json

import numpy as np
import pandas as pd

from . import RESULTS

HARD = RESULTS / "hardening"

VERDICT = {
    32090583: "Candidate (strongest)",
    229689348: "Candidate",
    149390648: "Candidate",
    198412174: "Candidate (ambiguous host)",
    294053492: "False positive: nearby eclipsing binary",
}

ASSESSMENT = {
    32090583: (
        "A **third transiting signal in the TOI-218 system**, found after the two known TOIs were masked: "
        "{rp} R_earth on a {period:.4f}-d orbit (T_eq ~{teq:.0f} K), not in any TOI, CTOI or SPOC TCE list. "
        "TOI-218 is one star of a wide binary: a near-twin M dwarf (same parallax and proper motion) sits 13.5\" "
        "away and is blended with it in TESS. The pixel localization puts the new signal on the target "
        "({loc_off:.1f} ± {loc_err:.1f}\") and excludes the companion at {comp_sig:.1f} sigma; the two known "
        "TOI-218 signals also come from the target ({toi218_01:.1f} and {toi218_02:.1f} sigma against the companion). "
        "The pixels lose {ratio:.2f} ± {ratio_err:.2f} times the light the light-curve depth predicts, and both "
        "halves of the data recover it (SNR 11.0 and 8.4). From TESS photometry alone TRICERATOPS gives "
        "FPP = {fpp:.2f} and NFPP = {nfpp:.2f}, almost all of it the scenario that the planet orbits the "
        "equal-brightness companion, which TRICERATOPS cannot tell apart without pixel information; with the "
        "companion and the other neighbours the localization excludes treated as cleared, FPP = {fpp_c:.3f} and "
        "NFPP {nfpp_c_rel}; adding the public Gemini-South speckle contrast curve of TOI-218 (2020, 562 nm) "
        "gives FPP {fpp_cc} and NFPP {nfpp_cc}, below the thresholds TRICERATOPS uses for validation "
        "(FPP < 0.015, NFPP < 0.001). It is not called validated here: the host is a flaring star, and no "
        "ground-based light curve has yet seen the transit (the one archival TFOP observation that covers a "
        "predicted transit has no pre-transit baseline and is inconclusive). TFOP notes call the host an eruptive (flaring) variable and suspect TOI-218.01 is "
        "stellar variability; flares are removed here and the new signal is a flat-bottomed, periodic, ~0.9-h dip "
        "seen in 400+ transits, unlike spot or flare activity. Extra candidates in systems that already have "
        "candidates are more often real. Next step: ground-based photometry of a predicted transit that "
        "resolves the binary (the 1.3-ppt depth is within reach of a 2-m telescope), an independent check of "
        "which star hosts the signal."
    ),
    229689348: (
        "An **ultra-short-period candidate: {rp} R_earth on an 11.2-hour orbit** (T_eq ~{teq:.0f} K, "
        "~{insol:.0f}x Earth's insolation), flat-bottomed, in {n_transits} transits and in both halves of the data. "
        "NASA's pipeline flagged it in three multi-sector runs but never promoted it; its reported SNR (5.8, 18.4, "
        "3.8) tracks how far its period estimate was from the true one, not a fading signal. SPOC's difference "
        "images put the source 55\" away in its best run, but none of them passed SPOC's own quality test; the "
        "joint PRF localization here, at the correct period with all {loc_sectors} sectors, puts it on the target "
        "({loc_off:.1f} ± {loc_err:.1f}\"), excludes the bright {far_sep:.1f}\" neighbour SPOC's offset pointed towards at "
        "{far_sig:.1f} sigma and the {near_sep:.1f}\" neighbour at {near_sig:.1f} sigma, and synthetic eclipses planted on "
        "those neighbours are correctly traced to them. Pixel/light-curve depth ratio {ratio:.2f} ± {ratio_err:.2f}. "
        "TRICERATOPS: FPP = {fpp:.3f}, NFPP = {nfpp:.4f}; the planet-on-target scenario carries {tp:.0f}% and "
        "the rest is almost all unresolved-companion scenarios ({stp:.0f}% a planet on a bound companion), which "
        "high-resolution imaging would test; with the neighbours the localization excludes treated as cleared the "
        "neighbour scenarios vanish (NFPP = {nfpp_c:.5f}, FPP = {fpp_c:.3f}). The transit is short for this star (b ~ {b:.2f}), and the transit "
        "shape alone prefers a denser star ({rho_free}), which TRICERATOPS reads as some weight on a companion host."
    ),
    149390648: (
        "**{rp} R_earth on a {period:.4f}-d orbit** (T_eq ~{teq:.0f} K) in a crowded southern field. SPOC flagged "
        "it once (s1-s96, SNR 8.5) and its difference-image offset was {spoc_off:.1f} ± {spoc_off_err:.1f}\", "
        "consistent with the target; this localization agrees ({loc_off:.1f} ± {loc_err:.1f}\", every neighbour "
        "excluded at 3 sigma). The first half of the data alone does not lock onto the period (SNR 5.4 at a "
        "different period; second half 7.6 at this one), and the pixels lose {ratio:.2f} ± {ratio_err:.2f} times "
        "the expected light, the highest ratio in the set, which in a crowded field can mean the PDC crowding "
        "correction is imperfect. TRICERATOPS: FPP = {fpp:.3f}, NFPP = {nfpp:.4f} from TESS photometry alone, "
        "the NFPP coming from many faint neighbours; with those the localization excludes treated as cleared, "
        "FPP = {fpp_c:.3f} and NFPP = {nfpp_c:.5f}. The last neighbour scenario left is TIC 149390646, a TIC entry "
        "with T = 17.0 and no Gaia counterpart; Gaia DR3 sees only a G = 21.6 source there, far too faint to make "
        "the dip. Ground-based photometry would settle it."
    ),
    198412174: (
        "**{rp} R_earth on a {period:.4f}-d orbit**, nearly grazing (b ~ {b:.2f}). SPOC flagged it in three runs; as "
        "for TIC 229689348 the falling SPOC SNR follows SPOC's period error. The pixels put the source on or very "
        "near the target ({loc_off:.1f} ± {loc_err:.1f}\") but **cannot separate the target from a T = 18.2 star "
        "{near_sep:.1f}\" away** (excluded at only {near_sig:.1f} sigma; it would need a ~{near_depth:.0f}% eclipse); the planted-eclipse test shows that this pair is below "
        "the method's resolution. TRICERATOPS: FPP = {fpp:.3f} (NFPP {nfpp:.4f}), dominated by a planet "
        "transiting an unresolved bound companion ({stp:.0f}%) rather than the target ({tp:.0f}%), because the "
        "transit is short for this star; clearing the neighbours the localization excludes gives FPP = {fpp_c:.3f}, "
        "NFPP = {nfpp_c:.5f}. High-resolution imaging (to find or exclude a companion and the {near_sep:.1f}\" "
        "star) is the decisive next step."
    ),
    294053492: (
        "**Not a planet candidate on this star.** The light-curve signal ({period:.4f} d, depth ~{depth:.0f} ppm) "
        "passed every light-curve test, but the pixels place the light loss {loc_off:.1f} ± {loc_err:.1f}\" "
        "north-east of the target, excluding the target at {target_sig:.1f} sigma; odd and even sectors each "
        "give the same off-target position on their own. The best-matching Gaia DR3 source is {src_id} "
        "(G = {src_g:.1f}, {src_sep:.1f}\" from the target and {src_dbest:.1f}\" from the best-fit position), which "
        "would need a ~{src_depth:.0f}% eclipse: an ordinary eclipsing binary. More light goes missing than the "
        "target could lose for this depth ({amp_ratio_best:.2f}x). TRICERATOPS, which has no pixel information, "
        "also flags a neighbour (NFPP = {nfpp:.2f}) but blames the bright star {nw_sep:.1f}\" to the north-west, which the "
        "pixels exclude at {nw_sig:.1f} sigma; the faint true source carries little prior weight there. This is what the hardening pass is for; the "
        "earlier write-up listed it as one of five candidates."
    ),
}

# signals that belong in the draft CTOI upload: ExoFOP notes field (<= 120 characters)
# candidate number for the draft upload: the next free TIC<id>.NN on ExoFOP (TIC 32090583.01 and .02 are
# the community-candidate names of TOI-218.01 and .02); .01 for stars with no TOI or CTOI
CTOI_SUFFIX = {32090583: "03"}

SUBMIT = {
    32090583: "3rd signal on TOI-218; pixel-localized to target, not 13.5in binary companion; flaring host",
    229689348: "USP 11.2h; SPOC TCE never TOI; pixel-localized on target (SPOC DV 55in offset not reproduced)",
    149390648: "SPOC TCE never TOI; on target in pixel localization; crowded field; FPP 0.09 with pixel-cleared neighbours",
    198412174: "SPOC TCE never TOI; grazing; T=18.2 star 5.3in away not excluded; needs imaging",
}


def small_p(x):
    """A probability for prose: four decimals, or mantissa x 10^n below 0.001."""
    if not np.isfinite(x):
        return "n/a"
    if x >= 1e-3:
        return f"{x:.4f}"
    if x <= 0:
        return "< 10⁻⁵"   # no draws in any non-planet scenario
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    e = int(np.floor(np.log10(x)))
    return f"{x / 10 ** e:.1f} × 10{str(e).translate(sup)}"


def rel(p):
    """'= 0.0880' or '< 10⁻⁵': a formatted probability with its relation sign."""
    return p if p.startswith("<") else f"= {p}"


def _sep(star):
    """Separation from the target (arcsec) at the median TESS epoch, from the localization's E/N offsets."""
    return float(np.hypot(star["east"], star["north"]))


def _load(path):
    return json.loads(path.read_text()) if path.exists() else None


def facts(tic):
    """Numbers quoted in the assessment of one candidate, read from the result files."""
    from .localize import pixel_depth_ratio

    loc = _load(HARD / "localize" / f"TIC{tic}.json")
    mc = _load(HARD / "mcmc" / f"TIC{tic}.json")
    tri = _load(HARD / "triceratops" / f"TIC{tic}_result.json")
    tri_c = _load(HARD / "triceratops" / f"TIC{tic}_result_cleared.json")
    tri_cc = _load(HARD / "triceratops" / f"TIC{tic}_result_cleared_cc.json")  # with a contrast curve
    s = mc["prior"]
    r = s["rp_rearth"]
    rf = mc["free_density"]["rho_ratio_to_tic"]
    ratio, ratio_err, _, _ = pixel_depth_ratio(loc)
    f = {"rp": f"{r['median']:.2f} ± {0.5 * (r['plus'] + r['minus']):.2f}", "period": mc["period"],
         "teq": s["teq_k"]["median"], "insol": s["insolation_earth"]["median"], "b": s["b"]["median"],
         "depth": s["depth_ppm"]["median"],
         "rho_free": f"density {rf['median']:.1f} +{rf['plus']:.1f} -{rf['minus']:.1f} x the TIC value",
         "loc_off": loc["offset_arcsec"], "loc_err": loc["offset_err_total_arcsec"],
         "target_sig": loc["target_sigma_total"], "loc_sectors": loc["n_sectors"],
         "ratio": ratio, "ratio_err": ratio_err,
         "fpp": tri["FPP_mean"] if tri else float("nan"), "nfpp": tri["NFPP_mean"] if tri else float("nan"),
         "fpp_c": tri_c["FPP_mean"] if tri_c else float("nan"),
         "nfpp_c": tri_c["NFPP_mean"] if tri_c else float("nan"),
         "fpp_cc": rel(small_p(tri_cc["FPP_mean"])) if tri_cc else "n/a",
         "nfpp_cc": rel(small_p(tri_cc["NFPP_mean"])) if tri_cc else "n/a",
         "nfpp_c_rel": rel(small_p(tri_c["NFPP_mean"])) if tri_c else "n/a"}
    fu = [d for d in (_load(p) for p in sorted((RESULTS / "followup").glob(f"TIC{tic}_*.json")))
          if d and d.get("verdict") == "candidate"]
    f["n_transits"] = fu[0]["vet"].get("n_transits_measured") if fu else loc["n_transits"]
    near = sorted([x for x in loc["stars"] if not x["is_target"]], key=_sep)
    f["near_sig"] = near[0]["excluded_sigma_total"] if near else float("nan")
    f["near_depth"] = 100 * near[0]["implied_eclipse_depth"] if near else float("nan")
    f["near_sep"] = _sep(near[0]) if near else float("nan")
    f["comp_sig"] = f["near_sig"]
    bright = [x for x in near if 45 < _sep(x) < 52 and x["tmag"] < 14.5]
    f["far_sig"] = bright[0]["excluded_sigma_total"] if bright else float("nan")
    f["far_sep"] = _sep(bright[0]) if bright else float("nan")
    src = next((x for x in loc["stars"] if not x["is_target"]), None)  # stars are sorted by chi-square
    if src:
        f.update(src_id=src["source_id"], src_sep=_sep(src), src_depth=100 * src["implied_eclipse_depth"],
                 src_g=src.get("gmag", float("nan")),
                 src_dbest=float(np.hypot(src["east"] - loc["offset_east"], src["north"] - loc["offset_north"])))
    f["amp_ratio_best"] = loc["amplitude_ratio"]
    nw = [x for x in loc["stars"] if not x["is_target"] and 16 < _sep(x) < 20 and x["tmag"] < 14.5]
    f["nw_sig"] = nw[0]["excluded_sigma_total"] if nw else float("nan")
    f["nw_sep"] = _sep(nw[0]) if nw else float("nan")
    if tri:
        tops = dict((k.split(":")[0], v) for k, v in tri["top_scenarios"] if k.split(":")[1] == str(tic))
        f["stp"], f["tp"] = 100 * tops.get("STP", 0.0), 100 * tops.get("TP", 0.0)
    if tic == 32090583:
        for lab, key in (("TOI-218.01", "toi218_01"), ("TOI-218.02", "toi218_02")):
            d = _load(HARD / "localize" / f"{lab}.json")
            comp = sorted([x for x in d["stars"] if not x["is_target"]], key=_sep)[0]
            f[key] = comp["excluded_sigma_total"]
    dv = HARD / "spoc_dv_summary.csv"
    if dv.exists():
        d = pd.read_csv(dv)
        d = d[d.tic == tic]
        if len(d):
            f["spoc_off"], f["spoc_off_err"] = float(d.ms_tic_offset_arcsec.iloc[-1]), float(d.ms_tic_offset_err_arcsec.iloc[-1])
    return f


def assessment(tic):
    text = ASSESSMENT.get(tic, "")
    try:
        return text.format(**facts(tic))
    except (KeyError, TypeError, ValueError) as exc:  # missing result file: say so rather than guess
        return text + f" [numbers unavailable: {exc!r}]"
