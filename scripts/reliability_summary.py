"""Summarise the injection and inversion tests under the final vetting rules.

    .venv/bin/python scripts/reliability_summary.py
"""
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from tess_search import RESULTS, vetting

PASS = ("candidate", "weak candidate")


def load(name):
    p = RESULTS / "reliability" / f"{name}.jsonl"
    return [json.loads(l) for l in p.read_text().splitlines() if l.strip()]


def verdict(s):
    return vetting.classify(s["vet"])[0] if s.get("vet") else s.get("verdict", "not vetted")


inj = load("inject")
rows = []
reasons = Counter()
for r in inj:
    match = next((s for s in r["signals"] if abs(s["period"] / r["period"] - 1) < 0.005), None) if r["recovered"] else None
    v = verdict(match) if match else None
    if match and v not in PASS and match.get("vet"):
        for x in vetting.classify(match["vet"])[1]:
            reasons[x.split(" (")[0][:70]] += 1
    rows.append({"period": r["period"], "rp": r["rp"], "b": r["b"], "depth_ppm": r["depth_true"] * 1e6,
                 "noise_ppm": r["noise_ppm"], "recovered": r["recovered"], "verdict": v,
                 "passed": v in PASS, "candidate": v == "candidate"})
df = pd.DataFrame(rows)
print(f"INJECTION: {len(df)} planets")
print(f"  detected at right period+epoch: {df.recovered.mean():.1%}")
print(f"  detected and passed vetting (candidate or weak): {df.passed.mean():.1%}")
print(f"  detected and full candidate: {df.candidate.mean():.1%}")
print(f"  vetting kept {df[df.recovered].passed.mean():.1%} of detected injected planets")
print("  why vetting rejected detected injected planets:", reasons.most_common(6))
for lo, hi in [(0.6, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 4.0)]:
    m = (df.rp >= lo) & (df.rp < hi)
    print(f"  Rp {lo}-{hi} R_earth: detected {df[m].recovered.mean():.0%}, passed {df[m].passed.mean():.0%} (n={m.sum()})")
for lo, hi in [(0.5, 2), (2, 7), (7, 20), (20, 40)]:
    m = (df.period >= lo) & (df.period < hi)
    print(f"  P {lo}-{hi} d: detected {df[m].recovered.mean():.0%}, passed {df[m].passed.mean():.0%} (n={m.sum()})")

inv = load("invert")
fc = [sum(verdict(s) == "candidate" for s in r["signals"]) for r in inv]
fw = [sum(verdict(s) == "weak candidate" for s in r["signals"]) for r in inv]
print(f"\nINVERSION: {len(inv)} flipped light curves")
print(f"  false candidates: {sum(fc)} on {sum(c > 0 for c in fc)} stars ({np.mean([c > 0 for c in fc]):.1%} of stars)")
print(f"  false weak candidates: {sum(fw)} on {sum(c > 0 for c in fw)} stars ({np.mean([c > 0 for c in fw]):.1%} of stars)")
for r, c, w in zip(inv, fc, fw):
    for s in r["signals"]:
        if verdict(s) in PASS:
            print(f"    TIC {r['tic']}: P={s['period']:.4f} snr_red={s.get('snr_red')} -> {verdict(s)}")
summary = {"inject_n": len(df), "inject_detected": df.recovered.mean(), "inject_passed": df.passed.mean(),
           "inject_candidate": df.candidate.mean(), "inject_vetting_keep": df[df.recovered].passed.mean(),
           "invert_n": len(inv), "invert_false_candidates": int(sum(fc)), "invert_false_weak": int(sum(fw)),
           "invert_stars_with_false_candidate": int(sum(c > 0 for c in fc))}
(RESULTS / "reliability" / "summary.json").write_text(json.dumps(summary, indent=1, default=float))
