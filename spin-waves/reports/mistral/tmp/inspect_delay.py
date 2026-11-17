import h5py, numpy as np
from collections import defaultdict
p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/delay_fine.h5"
with h5py.File(p, "r") as f:
    keys = sorted(f.keys())
    print("n_keys:", len(keys))
    print("first 8:", keys[:8])
    print("last 8:", keys[-8:])
    print("attrs:", list(f.attrs.keys()))
    # key field counts
    nfields = defaultdict(int)
    for k in keys:
        nfields[len(k.split("_"))] += 1
    print("field counts:", dict(nfields))
    if f.attrs.get("script_content"):
        sc = f.attrs["script_content"]
        if isinstance(sc, bytes): sc = sc.decode()
        open("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/delay_fine_script.py","w").write(sc)
        print("script saved, len", len(sc))
    tr = f[keys[0]]["trace"]
    print("trace shape:", tr.shape, tr.dtype)
    a = tr[:]
    print("ch0 (freq) head/mid/tail:", a[0,0,:2], a[0,0,500], a[0,0,-2:])
    print("ch1 max/min/mean:", a[0,1].max(), a[0,1].min(), a[0,1].mean())
    shapes = set()
    for k in keys[:200]: shapes.add(f[k]["trace"].shape)
    print("shapes sample:", shapes)
