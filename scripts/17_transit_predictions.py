"""Step 17: upcoming transits of the four candidates, and when they can be observed from the ground.

    .venv/bin/python scripts/17_transit_predictions.py [start YYYY-MM-DD] [end YYYY-MM-DD]

For every predicted transit between the start and end dates (default 2026-10-06 to 2027-03-31):
mid-transit time (BJD_TDB and UTC), its 1-sigma uncertainty from the ephemeris (results/candidates/
summary.csv), and, for each site of the Las Cumbres Observatory network (the telescopes most used for
TESS follow-up), whether the whole transit plus 30 minutes of baseline on each side happens with the
target above 30 degrees altitude and the Sun more than 12 degrees below the horizon.

Writes results/followup_planning/transits.csv (observable events only) and summary.json.
These are planning aids; observers should confirm with their own scheduling tools (e.g. TAPIR).
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time

from tess_search import DATA, RESULTS

OUT = RESULTS / "followup_planning"
SITES = {  # Las Cumbres Observatory network (approximate coordinates)
    "LCO Cerro Tololo (Chile)": (-30.1674, -70.8048, 2198),
    "LCO SAAO Sutherland (South Africa)": (-32.3797, 20.8107, 1798),
    "LCO Siding Spring (Australia)": (-31.2733, 149.0711, 1165),
    "LCO Teide (Tenerife)": (28.3003, -16.5117, 2390),
    "LCO McDonald (Texas)": (30.6797, -104.0247, 2070),
    "LCO Haleakala (Hawaii)": (20.7069, -156.2581, 3055),
}
MIN_ALT = 30.0
MAX_SUN_ALT = -12.0
BASELINE_H = 0.5


def bjd_to_utc(bjd, coord):
    """BJD_TDB -> UTC at the geocentre (two iterations of the barycentric light-travel time)."""
    t = Time(bjd, format="jd", scale="tdb")
    obs = t
    for _ in range(2):
        obs = Time(t.jd - obs.light_travel_time(coord, kind="barycentric",
                                                  location=EarthLocation.from_geocentric(0, 0, 0, u.m)).to_value(u.day),
                   format="jd", scale="tdb")
    return obs.utc


if __name__ == "__main__":
    start = Time(sys.argv[1] if len(sys.argv) > 1 else "2026-10-06", scale="utc")
    end = Time(sys.argv[2] if len(sys.argv) > 2 else "2027-03-31", scale="utc")
    OUT.mkdir(parents=True, exist_ok=True)
    cs = pd.read_csv(RESULTS / "candidates" / "summary.csv")
    cs = cs[cs.verdict.str.startswith("Candidate")]
    tic_tab = pd.read_parquet(DATA / "tic_15plus.parquet").set_index("tic")
    rows, summary = [], {}
    for _, c in cs.iterrows():
        t = tic_tab.loc[int(c.tic)]
        coord = SkyCoord(float(t.ra) * u.deg, float(t.dec) * u.deg)
        name = f"TIC {int(c.tic)}" + (f" ({c['name']})" if isinstance(c["name"], str) and c["name"] else "")
        n0 = int(np.ceil((start.tdb.jd - c.t0_bjd) / c.period))
        n1 = int(np.floor((end.tdb.jd - c.t0_bjd) / c.period))
        n = np.arange(n0, n1 + 1)
        tc_bjd = c.t0_bjd + n * c.period
        sig_min = np.hypot(c.t0_err, n * c.period_err) * 1440
        tc_utc = bjd_to_utc(tc_bjd, coord)
        half = (c.t14_h / 2 + BASELINE_H) / 24
        counts = {}
        for site, (lat, lon, h) in SITES.items():
            loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=h * u.m)
            # sample the window (ingress-baseline .. egress+baseline) at 9 points per event
            offs = np.linspace(-half, half, 9)
            times = Time(tc_utc.jd[:, None] + offs[None, :], format="jd", scale="utc")
            frame = AltAz(obstime=times.ravel(), location=loc)
            alt = coord.transform_to(frame).alt.deg.reshape(times.shape)
            sun = get_sun(times.ravel()).transform_to(frame).alt.deg.reshape(times.shape)
            ok = (alt.min(1) > MIN_ALT) & (sun.max(1) < MAX_SUN_ALT)
            counts[site] = int(ok.sum())
            for i in np.flatnonzero(ok):
                rows.append({"candidate": name, "tic": int(c.tic), "epoch": int(n[i]),
                             "mid_transit_bjd_tdb": round(float(tc_bjd[i]), 5),
                             "mid_transit_utc": tc_utc[i].iso[:16], "sigma_min": round(float(sig_min[i]), 1),
                             "t14_h": round(float(c.t14_h), 2), "depth_ppm": round(float(c.depth_ppm)),
                             "site": site,
                             "window_start_utc": Time(tc_utc.jd[i] - half, format="jd").iso[:16],
                             "window_end_utc": Time(tc_utc.jd[i] + half, format="jd").iso[:16],
                             "min_altitude_deg": round(float(alt[i].min()), 1)})
        summary[name] = {"period_d": float(c.period), "transits_in_window": int(len(n)),
                         "timing_sigma_min_range": [round(float(sig_min.min()), 1), round(float(sig_min.max()), 1)],
                         "observable_per_site": counts}
        print(name, len(n), "transits;", {k.split(" (")[0].replace("LCO ", ""): v for k, v in counts.items()}, flush=True)
    df = pd.DataFrame(rows).sort_values(["candidate", "mid_transit_utc", "site"])
    df.to_csv(OUT / "transits.csv", index=False)
    (OUT / "summary.json").write_text(json.dumps({"start_utc": start.iso[:10], "end_utc": end.iso[:10],
                                                  "min_altitude_deg": MIN_ALT, "max_sun_altitude_deg": MAX_SUN_ALT,
                                                  "baseline_h": BASELINE_H, "candidates": summary}, indent=1))
    print(f"{len(df)} observable (event, site) pairs -> {OUT / 'transits.csv'}")
