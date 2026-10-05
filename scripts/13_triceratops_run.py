"""Step 13b: TRICERATOPS false-positive probabilities (run with .venv-tri/bin/python).

For each candidate TRICERATOPS (Giacalone et al. 2021) weighs every way the dip
could arise: a planet on the target (TP), an eclipsing binary on the target (EB),
an unresolved bound or background companion (PTP/PEB/STP/SEB/DTP/DEB/BTP/BEB),
or a resolved neighbouring star (NTP/NEB), each with priors from the field's
star population and a fit to the transit shape. It reports
  FPP  = probability the signal is NOT a planet on the target
  NFPP = probability it comes from a resolved neighbour.
With --cleared, neighbours that the pixel-level localization (scripts/12_localize.py)
excludes as the source at more than 3 sigma are treated as cleared, the way stars cleared by
ground-based photometry are: their required depth is set to zero so TRICERATOPS no longer
considers them, and the result goes to <label>_result_cleared.json.

Giacalone et al. call a candidate statistically validated at FPP < 0.015 and
NFPP < 0.001 (with high-resolution imaging), likely planet at FPP < 0.5 and
NFPP < 0.001. No imaging contrast curves are used here, so unresolved companions
are only limited by Gaia; these FPPs are therefore upper-end values.

The calculation is Monte Carlo; it is repeated N_RUNS (5) times and the mean and
scatter reported. Results: results/hardening/triceratops/<label>_result.json.
"""
import json
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np

import triceratops.triceratops as tr

ROOT = Path(__file__).resolve().parent.parent
BG_DIR = ROOT / "results" / "hardening" / "triceratops" / "background"


def gaia_background_vizier(ra, dec, tic, field=0.1, mag_limit=21.0):
    """TRICERATOPS's field-star population (its query_gaia_background), with the Gaia DR3
    query sent to VizieR instead of the ESA archive, which hung for >10 min in Oct 2026.
    The conversion from Gaia photometry to stellar properties is TRICERATOPS's own code
    (dwarf sequence of Pecaut & Mamajek 2013), repeated here unchanged."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier
    from math import pi
    from pandas import DataFrame
    from triceratops.funcs import G, Msun, Rsun, gaia_to_Tmag, mamajek_column

    fname = BG_DIR / f"{tic}_gaia_background.csv"
    if fname.exists():
        return str(fname)
    radius = np.sqrt(field / pi)
    v = Vizier(columns=["Gmag", "BPmag", "RPmag", "Plx", "e_Plx"], row_limit=-1,
               column_filters={"Gmag": f"<{mag_limit}"})
    t = v.query_region(SkyCoord(ra * u.deg, dec * u.deg), radius=radius * u.deg, catalog="I/355/gaiadr3")[0]
    Gmag = np.ma.filled(t["Gmag"], np.nan).astype(float)
    BPmag = np.ma.filled(t["BPmag"], np.nan).astype(float)
    RPmag = np.ma.filled(t["RPmag"], np.nan).astype(float)
    plx = np.ma.filled(t["Plx"], np.nan).astype(float)
    plx_snr = plx / np.ma.filled(t["e_Plx"], np.nan).astype(float)
    Tmag = gaia_to_Tmag(Gmag, BPmag, RPmag)

    def seq_table(xcol, ycol):
        x, y = mamajek_column(xcol), mamajek_column(ycol)
        keep = np.isfinite(x) & np.isfinite(y)
        x, y = x[keep], y[keep]
        order = np.argsort(x)
        x, y = x[order], y[order]
        keep = np.concatenate(([True], np.diff(x) > 1e-9))
        return x[keep], y[keep]

    cx, cy = seq_table("BpRp", "Teff")
    Teff = np.interp(BPmag - RPmag, cx, cy, left=cy[0], right=cy[-1])
    vals = {}
    for col in ("Msun", "R_Rsun", "M_G", "M_J", "M_Ks", "JH"):
        tx, ty = seq_table("Teff", col)
        vals[col] = np.interp(Teff, tx, ty)
    logg = np.log10(G * (vals["Msun"] * Msun) / (vals["R_Rsun"] * Rsun) ** 2)
    use_plx = np.isfinite(plx) & (plx > 0.0) & (plx_snr > 5.0)
    dmod = np.where(use_plx, 5.0 * np.log10(1000.0 / np.where(plx > 0.0, plx, np.nan)) - 5.0, Gmag - vals["M_G"])
    good = np.isfinite(Teff) & np.isfinite(vals["Msun"]) & np.isfinite(dmod) & np.isfinite(Tmag)
    df = DataFrame(np.column_stack([
        vals["Msun"][good], logg[good], np.log10(Teff[good]), np.zeros(good.sum()), Tmag[good],
        (vals["M_J"] + dmod)[good], (vals["M_J"] - vals["JH"] + dmod)[good], (vals["M_Ks"] + dmod)[good]]),
        columns=["Mact", "logg", "logTe", "[M/H]", "TESS", "J", "H", "Ks"])
    BG_DIR.mkdir(parents=True, exist_ok=True)
    fname = BG_DIR / f"{tic}_gaia_background.csv"
    df.to_csv(fname, index=False)
    print(f"TIC {tic}: {len(df)} Gaia DR3 field stars (VizieR, G < {mag_limit})", flush=True)
    return str(fname)
OUT = ROOT / "results" / "hardening" / "triceratops"
CACHE = ROOT / "data" / "tess_cache"
N_RUNS = int(os.environ.get("TRI_RUNS", 5))  # run-to-run FPP scatter is ~1%, so 5 runs suffice
N_DRAWS = 1_000_000
LOC_DIR = ROOT / "results" / "hardening" / "localize"
CLEAR_SIGMA = 3.0


def clear_neighbours(target, tic):
    """Zero the required depth of TRICERATOPS stars that the pixel localization excludes.

    TRICERATOPS stars are TIC entries (offsets from the target from their TIC positions);
    localization stars are Gaia DR3 sources. Match by sky offset (within 4") and brightness."""
    loc = json.loads((LOC_DIR / f"TIC{tic}.json").read_text())
    gstars = [s_ for s_ in loc["stars"] if not s_["is_target"]]
    st = target.stars
    cleared = []
    for i in range(1, len(st)):
        if st["tdepth"].iloc[i] <= 0:
            continue
        sep, pa = float(st["sep (arcsec)"].iloc[i]), float(st["PA (E of N)"].iloc[i])
        e, n = sep * np.sin(np.radians(pa)), sep * np.cos(np.radians(pa))
        best = min(gstars, key=lambda g: np.hypot(g["east"] - e, g["north"] - n), default=None)
        if best is None or np.hypot(best["east"] - e, best["north"] - n) > 4.0:
            continue
        if abs(best["tmag"] - float(st["Tmag"].iloc[i])) > 1.5:
            continue
        if best["excluded_sigma_total"] > CLEAR_SIGMA:
            cleared.append({"tic": int(st["ID"].iloc[i]), "gaia": best["source_id"], "sep": sep,
                            "tmag": float(st["Tmag"].iloc[i]), "excluded_sigma": best["excluded_sigma_total"]})
            target.stars.iloc[i, target.stars.columns.get_loc("tdepth")] = 0.0
    return cleared

if __name__ == "__main__":
    cleared_mode = "--cleared" in sys.argv
    order = [int(a) for a in sys.argv[1:] if a != "--cleared"]  # optional: TIC IDs to run, in this order
    CACHE.mkdir(parents=True, exist_ok=True)
    paths = sorted(OUT.glob("*_input.json"))
    if order:
        by_tic = {json.loads(p.read_text())["tic"]: p for p in paths}
        paths = [by_tic[t] for t in order if t in by_tic]
    for path in paths:
        d = json.loads(path.read_text())
        out_path = OUT / f"{d['label']}_result{'_cleared' if cleared_mode else ''}.json"
        if out_path.exists():
            print(f"{d['label']}: done already", flush=True)
            continue
        t_start = time.time()
        try:
            sectors = np.array([int(s) for s in d["apertures"]])
            aps = [np.array(d["apertures"][str(s)]) for s in sectors]
            bg = gaia_background_vizier(d["ra"], d["dec"], d["tic"])
            target = tr.target(ID=d["tic"], sectors=sectors, search_radius=10,
                               lightkurve_cache_dir=str(CACHE), trilegal_fname=bg)
            target.calc_depths(tdepth=d["depth"], all_ap_pixels=aps)
            cleared = clear_neighbours(target, d["tic"]) if cleared_mode else []
            if cleared_mode:
                print(f"{d['label']}: cleared by pixel localization: {cleared}", flush=True)
            runs = []
            for i in range(N_RUNS):
                target.calc_probs(time=np.array(d["time"]), flux_0=np.array(d["flux"]), flux_err_0=d["flux_err"],
                                  P_orb=d["period"], N=N_DRAWS, parallel=True, verbose=0)
                runs.append({"FPP": float(target.FPP), "NFPP": float(target.NFPP),
                             "degenerate": bool(getattr(target, "FPP_degenerate", False)),
                             "probs": target.probs[["ID", "scenario", "prob"]].to_dict("records")})
                print(f"{d['label']} run {i + 1}: FPP {target.FPP:.4f} NFPP {target.NFPP:.5f}", flush=True)
            fpp = np.array([r["FPP"] for r in runs])
            nfpp = np.array([r["NFPP"] for r in runs])
            # scenario probabilities averaged over runs
            scen = {}
            for r in runs:
                for p in r["probs"]:
                    key = f"{p['scenario']}:{p['ID']}"
                    scen[key] = scen.get(key, 0.0) + p["prob"] / len(runs)
            top = sorted(scen.items(), key=lambda kv: -kv[1])[:8]
            res = {"tic": d["tic"], "label": d["label"], "sectors": sectors.tolist(), "n_runs": N_RUNS, "n_draws": N_DRAWS,
                   "FPP_mean": float(fpp.mean()), "FPP_std": float(fpp.std(ddof=1)),
                   "NFPP_mean": float(nfpp.mean()), "NFPP_std": float(nfpp.std(ddof=1)),
                   "top_scenarios": top, "runs": [{k: v for k, v in r.items() if k != "probs"} for r in runs],
                   "stars": [],
                   "background_population": str(target.trilegal_fname), "runtime_s": round(time.time() - t_start),
                   "cleared_by_localization": cleared}
            out_path.write_text(json.dumps(res, indent=1))  # save before anything else can fail
            try:
                cols = [c for c in ("ID", "Tmag", "sep (arcsec)", "PA (E of N)", "fluxratio", "tdepth")
                        if c in target.stars.columns]
                res["stars"] = json.loads(target.stars[cols].to_json(orient="records"))
                out_path.write_text(json.dumps(res, indent=1))
            except Exception:
                print(f"{d['label']}: star table not saved:\n{traceback.format_exc()}", flush=True)
            print(f"{d['label']}: FPP {fpp.mean():.4f} +- {fpp.std(ddof=1):.4f}, NFPP {nfpp.mean():.5f} "
                  f"({time.time() - t_start:.0f} s)", flush=True)
        except Exception:
            print(f"{d['label']} failed:\n{traceback.format_exc()}", flush=True)
