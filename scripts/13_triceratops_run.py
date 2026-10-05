"""Step 13b: TRICERATOPS false-positive probabilities (run with .venv-tri/bin/python).

For each candidate TRICERATOPS (Giacalone et al. 2021) weighs every way the dip
could arise: a planet on the target (TP), an eclipsing binary on the target (EB),
an unresolved bound or background companion (PTP/PEB/STP/SEB/DTP/DEB/BTP/BEB),
or a resolved neighbouring star (NTP/NEB), each with priors from the field's
star population and a fit to the transit shape. It reports
  FPP  = probability the signal is NOT a planet on the target
  NFPP = probability it comes from a resolved neighbour.
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

if __name__ == "__main__":
    only = {int(a) for a in sys.argv[1:]}
    CACHE.mkdir(parents=True, exist_ok=True)
    for path in sorted(OUT.glob("*_input.json")):
        d = json.loads(path.read_text())
        if only and d["tic"] not in only:
            continue
        out_path = OUT / f"{d['label']}_result.json"
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
            runs = []
            for i in range(N_RUNS):
                target.calc_probs(time=np.array(d["time"]), flux_0=np.array(d["flux"]), flux_err_0=d["flux_err"],
                                  P_orb=d["period"], N=N_DRAWS, parallel=True, verbose=0)
                runs.append({"FPP": float(target.FPP), "NFPP": float(target.NFPP),
                             "degenerate": bool(getattr(target, "FPP_degenerate", False)),
                             "probs": target.probs[["ID", "scenario", "prob"]].to_dict("records")})
                print(f"{d['label']} run {i + 1}: FPP {target.FPP:.4f} NFPP {target.NFPP:.5f}", flush=True)
            stars = target.stars[["ID", "Tmag", "sep (arcsec)", "PA (deg)", "fluxratio", "tdepth"]].copy()
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
                   "stars": json.loads(stars.to_json(orient="records")),
                   "background_population": str(target.trilegal_fname), "runtime_s": round(time.time() - t_start)}
            out_path.write_text(json.dumps(res, indent=1))
            print(f"{d['label']}: FPP {fpp.mean():.4f} +- {fpp.std(ddof=1):.4f}, NFPP {nfpp.mean():.5f} "
                  f"({time.time() - t_start:.0f} s)", flush=True)
        except Exception:
            print(f"{d['label']} failed:\n{traceback.format_exc()}", flush=True)
