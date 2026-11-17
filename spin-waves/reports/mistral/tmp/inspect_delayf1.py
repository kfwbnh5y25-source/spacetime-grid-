import h5py, numpy as np
from collections import defaultdict
p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/delay_extract/delay_f1.h5"
with h5py.File(p, "r") as f:
    keys = sorted(f.keys())
    print("n_keys:", len(keys), "attrs:", list(f.attrs.keys()))
    rec = defaultdict(set)
    shapes = defaultdict(int)
    bad = 0
    for k in keys:
        pr = k.split("_")
        if len(pr) == 4:
            delay, f1, ph, ro = pr
            rec[(float(delay), float(f1))].add(float(ph))
            shapes[f[k]["trace"].shape] += 1
        else:
            bad += 1
    print("bad keys:", bad)
    print("shapes:", dict(shapes))
    delays = sorted(set(d for d, _ in rec)); f1s = sorted(set(f for _, f in rec))
    print("n delays:", len(delays))
    print("delays:", delays)
    print("n f1:", len(f1s), "f1:", f1s)
    pos = [d for d in delays if d >= 0]
    print("positive delays:", pos)
    # completeness: phases per (delay,f1) for a hot f1
    for f1q in (4072.0, 4091.0):
        rows = []
        for d in pos:
            key = (d, float(f1q))
            if key in rec:
                rows.append(f"{d:.3f}:{len(rec[key])}")
        print(f"f1={f1q} positive-delay phase counts:", " ".join(rows))
    k0 = keys[len(keys)//2]
    tr = f[k0]["trace"][:]
    print("example", k0, tr.shape, "ch0:", tr[0,0,0], tr[0,0,500], "ch1 max:", tr[:,1,:].max())
