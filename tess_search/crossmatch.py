"""Is a detection already known? Compare against every public list we can get.

Sources (all downloaded into data/catalogs/ by scripts/03_catalogs.py):
    TOI        TESS Objects of Interest (ExoFOP)
    CTOI       Community TOIs (ExoFOP)
    CONFIRMED  NASA Exoplanet Archive confirmed planets (pscomppars)
    TCE        every SPOC Threshold Crossing Event, single- and multi-sector
               (signals NASA's pipeline found, including ones never promoted to TOIs)
    EB         TESS eclipsing binary catalogue (Prsa et al. 2022, VizieR J/ApJS/258/16)

A detection "matches" a known signal if the periods agree within 0.3%, or are
in a simple ratio (1:2, 2:1, 1:3, 3:1, 2:3, 3:2) within 0.3%, since pipelines
often report half or double the true period.
"""
import glob
import re

import numpy as np
import pandas as pd

from . import DATA

CAT = DATA / "catalogs"
RATIOS = {1.0: "same", 0.5: "half", 2.0: "double", 1 / 3: "third", 3.0: "triple",
          2 / 3: "2/3", 1.5: "3/2"}
TOL = 0.003


def _tic_int(series):
    return pd.to_numeric(series.astype(str).str.replace("TIC", "").str.strip(), errors="coerce")


def load_tois():
    df = pd.read_csv(CAT / "toi.csv")
    return pd.DataFrame({
        "tic": _tic_int(df["TIC ID"]), "period": pd.to_numeric(df["Period (days)"], errors="coerce"),
        "epoch_btjd": pd.to_numeric(df["Epoch (BJD)"], errors="coerce") - 2457000,
        "name": "TOI-" + df["TOI"].astype(str),
        "disposition": df["TFOPWG Disposition"].fillna(df["TESS Disposition"]).fillna(""),
        "depth_ppm": pd.to_numeric(df["Depth (ppm)"], errors="coerce"), "source": "TOI",
    })


def load_ctois():
    df = pd.read_csv(CAT / "ctoi.csv")
    return pd.DataFrame({
        "tic": _tic_int(df["TIC ID"]), "period": pd.to_numeric(df["Period (days)"], errors="coerce"),
        "epoch_btjd": pd.to_numeric(df["Transit Epoch (BJD)"], errors="coerce") - 2457000,
        "name": "CTOI-" + df["CTOI"].astype(str),
        "disposition": df["User Disposition"].fillna(""),
        "depth_ppm": pd.to_numeric(df["Depth ppm"], errors="coerce"), "source": "CTOI",
    })


def load_confirmed():
    df = pd.read_csv(CAT / "confirmed.csv")
    return pd.DataFrame({
        "tic": _tic_int(df["tic_id"]), "period": pd.to_numeric(df["pl_orbper"], errors="coerce"),
        "epoch_btjd": pd.to_numeric(df["pl_tranmid"], errors="coerce") - 2457000,
        "name": df["pl_name"], "disposition": "CONFIRMED", "depth_ppm": np.nan, "source": "CONFIRMED",
    })


def load_tces(tics=None, cache=CAT / "tce_targets.parquet"):
    if cache.exists():
        df = pd.read_parquet(cache)
    else:
        parts = []
        for path in sorted(glob.glob(str(CAT / "tce" / "*_dvr-tcestats.csv"))):
            m = re.search(r"-(s\d{4}-s\d{4})_", path)
            t = pd.read_csv(path, comment="#", usecols=lambda c: c in {
                "ticid", "tce_plnt_num", "tce_period", "tce_time0bt", "tce_duration",
                "tce_depth", "tce_model_snr", "tce_num_transits"})
            t["run"] = m.group(1) if m else path
            if tics is not None:
                t = t[t.ticid.isin(set(tics))]
            parts.append(t)
        df = pd.concat(parts, ignore_index=True)
        df.to_parquet(cache)
    return pd.DataFrame({
        "tic": df.ticid, "period": df.tce_period, "epoch_btjd": df.tce_time0bt,
        "name": "TCE " + df.run + "-" + df.tce_plnt_num.astype(str),
        "disposition": "snr=" + df.tce_model_snr.round(1).astype(str),
        "depth_ppm": df.tce_depth, "source": "TCE",
    })


def load_ebs():
    path = CAT / "tess_ebs.csv"
    if not path.exists():
        from astroquery.vizier import Vizier

        v = Vizier(columns=["TIC", "Per", "BJD0", "Morph"], row_limit=-1)
        tab = v.get_catalogs("J/ApJS/258/16")[0]
        tab.to_pandas().to_csv(path, index=False)
    df = pd.read_csv(path)
    return pd.DataFrame({
        "tic": pd.to_numeric(df["TIC"], errors="coerce"), "period": pd.to_numeric(df["Per"], errors="coerce"),
        "epoch_btjd": pd.to_numeric(df["BJD0"], errors="coerce") - 2457000,
        "name": "EB TIC " + df["TIC"].astype(str), "disposition": "morph=" + df["Morph"].astype(str),
        "depth_ppm": np.nan, "source": "EB",
    })


def load_known(tics=None):
    frames = [load_tois(), load_ctois(), load_confirmed(), load_tces(tics), load_ebs()]
    known = pd.concat(frames, ignore_index=True)
    known = known.dropna(subset=["tic"])
    known["tic"] = known.tic.astype(np.int64)
    if tics is not None:
        known = known[known.tic.isin(set(int(t) for t in tics))]
    return known.reset_index(drop=True)


def period_relation(p_det, p_known, tol=TOL):
    """Name of the period relation ('same', 'half', ...) or None."""
    if not (np.isfinite(p_det) and np.isfinite(p_known)) or p_known <= 0:
        return None
    r = p_det / p_known
    for ratio, name in RATIOS.items():
        if abs(r / ratio - 1) < tol:
            return name
    return None


def match(tic, period, known):
    """All known signals on this star whose period relates to `period`."""
    rows = known[known.tic == int(tic)]
    hits = []
    for _, k in rows.iterrows():
        rel = period_relation(period, k.period)
        if rel:
            hits.append({"source": k.source, "name": k["name"], "relation": rel,
                         "known_period": k.period, "disposition": k.disposition})
    return hits


def novelty(hits):
    """Summarise matches: 'known planet/TOI', 'known EB', 'SPOC TCE only', or 'new'."""
    srcs = {h["source"] for h in hits}
    if srcs & {"CONFIRMED", "TOI", "CTOI"}:
        return "known (planet/TOI/CTOI)"
    if "EB" in srcs:
        return "known EB"
    if "TCE" in srcs:
        return "SPOC TCE only (never promoted)"
    return "new"
