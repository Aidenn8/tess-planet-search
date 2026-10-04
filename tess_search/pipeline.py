"""Everything for one star, start to finish: prepare -> search -> vet -> crossmatch -> plot."""
import json
import time

import numpy as np

from . import RESULTS, crossmatch, lightcurve, report, search, vetting

PLOT_SNR = 6.0    # make a diagnostic sheet for every detection at least this strong
VET_SNR = 6.0     # vet everything at least this strong (weaker ones are recorded but not vetted)

_KNOWN = None


def known_table(tics=None):
    global _KNOWN
    if _KNOWN is None:
        _KNOWN = crossmatch.load_known(tics)
    return _KNOWN


def process_star(row, known=None, raw=None, out_dir=None, plots=True, max_signals=5):
    """Run the full pipeline on one target row (from data/targets.csv). Returns the record."""
    t_start = time.time()
    tic = int(row["tic"])
    known = known if known is not None else known_table()
    lc = lightcurve.prepare(tic, raw=raw)
    r_star = float(row["rad"])
    m_star = float(row["mass"]) if np.isfinite(row.get("mass", np.nan)) else r_star
    dets, periods, first_sde = search.search_star(lc, r_star, m_star, max_signals=max_signals)

    signals = []
    for k, det in enumerate(dets):
        d = det.as_dict()
        entry = {"detection": d}
        if d["snr"] >= VET_SNR:
            v = vetting.vet(lc, d, r_star, m_star)
            verdict, reasons = vetting.classify(v)
            hits = crossmatch.match(tic, d["period"], known)
            entry.update(vet=v, verdict=verdict, reasons=reasons, matches=hits,
                         novelty=crossmatch.novelty(hits))
            if plots and d["snr"] >= PLOT_SNR and out_dir is not None:
                png = out_dir / "plots" / f"TIC{tic}_{k + 1}.png"
                png.parent.mkdir(parents=True, exist_ok=True)
                report.plot_detection(lc, d, v, verdict, reasons, hits, periods, first_sde, path=png,
                                      title_extra=f"signal {k + 1}: {verdict}, {entry['novelty']}")
                entry["plot"] = str(png.relative_to(out_dir))
        signals.append(entry)

    return {
        "tic": tic, "r_star": r_star, "m_star": m_star, "teff": float(row["Teff"]),
        "tmag": float(row["Tmag"]), "contratio": float(row.get("contratio", np.nan)),
        "n_points": int(len(lc.time)), "n_sectors": int(len(np.unique(lc.sector))),
        "sectors": sorted(int(s) for s in np.unique(lc.sector)),
        "n_seasons": int(len(np.unique(lc.season))), "time_span": float(np.ptp(lc.time)),
        "n_trial_periods": int(len(periods)), "diagnostics": lc.star_meta,
        "signals": signals, "runtime_s": round(time.time() - t_start, 1),
    }


def save(record, out_dir):
    path = out_dir / "search" / f"{record['tic']}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(record, default=_json_default))
    tmp.rename(path)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


DEFAULT_OUT = RESULTS
