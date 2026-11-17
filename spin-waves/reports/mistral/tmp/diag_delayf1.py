"""Diagnose anomalies in delay_f1.h5 phase sweeps before trusting any tau fit.

Suspects:
  A) keys with readoff suffix other than "1" (e.g. "1EX") = different SA settings
  B) baseline-recalibration steps (every ~90 s = ~180 sweeps) landing inside a
     (delay,f1) phase block -> random subset of phases offset -> fake |c1|
  C) genuine physics
"""
import numpy as np
import h5py
from collections import defaultdict

BASE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/delay_extract"

with h5py.File(f"{BASE}/delay_f1.h5", "r") as f:
    all_keys = sorted(f.keys())
    suffixes = defaultdict(int)
    for k in all_keys:
        suffixes[k.split("_")[-1]] += 1
    print("readoff suffixes:", dict(suffixes))

    # group per (f1, delay) for f1=4069 and f1=4074 (hot band, clean-ish)
    for F1 in (4069, 4074):
        rec = defaultdict(dict)
        for k in all_keys:
            d, f1, ph, ro = k.split("_")
            if ro != "1" or float(f1) != F1:
                continue
            rec[float(d)][float(ph)] = k
        ds = sorted(rec)
        print(f"\n### f1={F1}: {len(ds)} delays")
        # noise scale: per-bin noise from a flat region, estimate from a few sweeps
        def I_of(k):
            tr = f[k]["trace"][:]
            s = tr[:, 1, :].astype(np.float64)
            return s.sum() / s.shape[0], s
        # print sorted table with phase-resolved I for selected delays
        sel = [-100.0, -68.4211, -11.7188, -4.8438, -2.0312, -0.1562, 0.0, 1.0938,
               2.0312, 2.5, 3.125, 4.5, 4.8438, 5.7812, 6.9062, 7.9688, 13.2812]
        print("t_gap |   A0   |   A1   |  phi0_deg | I(16 phases sorted by phi)")
        for d in sel:
            if d not in rec:
                continue
            phs = np.array(sorted(rec[d]))
            Is, one = [], None
            for p in phs:
                v, s = I_of(rec[d][p])
                Is.append(v)
                one = s
            I = np.array(Is)
            u = np.exp(-1j * phs)
            c1 = np.mean(I * u)
            A0, A1, p0 = I.mean(), 2 * abs(c1), np.degrees(np.angle(c1))
            # line vs floor split: top-8-bin sum around max bin (robust to floor steps?)
            # also step diagnostic: max |I_i - median| / (spread of smooth cos)
            mx = one[0].argmax()
            print(f"{d:+8.3f} | {A0*1e6:6.1f} | {A1*1e6:6.1f} | {p0:8.1f} | "
                  + " ".join(f"{v*1e6:5.1f}" for v in I))
