"""Compare pipeline results on known-planet hosts with the catalogues.

    .venv/bin/python scripts/check_known.py [results_dir]

For every processed star that hosts a confirmed planet, TOI or CTOI, prints the
known periods, what the search found, whether each known signal was recovered,
and the vetting verdict it got. Known planets should mostly pass vetting;
signals the TESS team already rejected (FP/FA) should mostly fail it.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from tess_search import DATA, RESULTS, crossmatch

out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else RESULTS
targets = pd.read_csv(DATA / "targets.csv")
known = crossmatch.load_known(targets.tic.tolist())
planets = known[known.source.isin(["TOI", "CTOI", "CONFIRMED"])]

rows = []
for path in sorted((out_dir / "search").glob("*.json")):
    r = json.loads(path.read_text())
    k = planets[planets.tic == r["tic"]]
    if k.empty:
        continue
    sigs = r["signals"]
    print(f"TIC {r['tic']} ({r['n_sectors']} sectors, noise {r['diagnostics']['noise_ppm']:.0f} ppm/2-min)")
    for _, kk in k.drop_duplicates("period").iterrows():
        found, verdict = "MISSED", ""
        for i, s in enumerate(sigs):
            rel = crossmatch.period_relation(s["detection"]["period"], kk.period)
            if rel:
                found = f"{rel} (signal {i + 1}, SNR {s['detection']['snr']:.1f})"
                verdict = s.get("verdict", "not vetted")
                break
        print(f"   known {kk['name']:<22} P={kk.period:9.4f} d  {kk.disposition:<10} -> {found:<28} {verdict}")
        rows.append({"tic": r["tic"], "name": kk["name"], "period": kk.period, "source": kk.source,
                     "disposition": kk.disposition, "recovered": found != "MISSED", "verdict": verdict})
    for s in sigs:
        d = s["detection"]
        print(f"   found P={d['period']:9.5f} d  SNR={d['snr']:5.1f}  depth={d['depth'] * 1e6:6.0f} ppm  "
              f"-> {s.get('verdict', '-'):<15} {s.get('novelty', ''):<26} {'; '.join(s.get('reasons', []))[:110]}")
if rows:
    df = pd.DataFrame(rows).drop_duplicates(["tic", "period"])
    in_range = df[(df.period >= 0.4) & (df.period <= 40)]
    print(f"\nrecovered {in_range.recovered.sum()}/{len(in_range)} known signals with 0.4 < P < 40 d")
    df.to_csv(out_dir / "known_recovery.csv", index=False)
