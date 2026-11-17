"""Phase-table diagnostic on delay_delta.h5: is persistent A1 at low Delta real smooth cos?"""
import numpy as np
import h5py
from collections import defaultdict

BASE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/delay_extract"
SEL = [50.0, 58.0, 62.0, 66.0, 70.0, 74.0, 88.0]
TD = [0.0, 2.0312, 3.125, 4.5312, 5.8594, 6.9531]

with h5py.File(f"{BASE}/delay_delta.h5", "r") as f:
    rec = defaultdict(dict)
    for k in sorted(f.keys()):
        pr = k.split("_")
        if len(pr) != 4 or pr[3] != "1":
            continue
        rec[(float(pr[0]), float(pr[1]))][float(pr[2])] = k
    for D in SEL:
        f1 = 4106 - D / 2
        # find delays present
        ds = sorted(d for (d, b) in rec if b == D)
        if not ds:
            print(f"Delta={D}: no keys"); continue
        print(f"\n=== Delta={D} (f1={f1:.0f}, f2={4106 + D/2:.0f}) ===")
        print("t_gap  |  A0  | depth% | phi0deg |  I (16 phases)")
        for t in TD:
            d = min(ds, key=lambda x: abs(x - t))
            if abs(d - t) > 0.2 or (d, D) not in rec:
                continue
            pdict = rec[(d, D)]
            phs = np.array(sorted(pdict))
            tr0 = f[pdict[phs[0]]]["trace"][:]
            I = np.zeros(len(phs))
            for i, p in enumerate(phs):
                s = f[pdict[p]]["trace"][:, 1, :].astype(np.float64)
                I[i] = s.sum() / s.shape[0]
            u = np.exp(-1j * phs)
            c1 = np.mean(I * u)
            A0, A1 = I.mean(), 2 * abs(c1)
            depth = (I.max() - I.min()) / (I.max() + I.min()) * 100 if A0 > 0 else 0
            # residual after best cos: unstructured?
            fit = A0 + A1 * np.cos(phs - np.angle(c1))
            rms = np.sqrt(np.mean((I - fit) ** 2))
            print(f"{d:+7.3f} | {A0*1e6:5.0f} | {depth:5.1f} | {np.degrees(np.angle(c1)):7.1f} | rms/resid={rms*1e6:5.1f} " + " ".join(f"{v*1e6:5.0f}" for v in I))
