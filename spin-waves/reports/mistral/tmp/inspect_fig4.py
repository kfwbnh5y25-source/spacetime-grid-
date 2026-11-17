import h5py, numpy as np, sys

p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/detune_fine031126.h5"
with h5py.File(p, "r") as f:
    keys = sorted(f.keys())
    print("n_keys:", len(keys))
    print("first 5:", keys[:5])
    print("last 5:", keys[-5:])
    # parse keys
    f1s, f2s, phs, ro = set(), set(), set(), set()
    bad = 0
    for k in keys:
        parts = k.split("_")
        if len(parts) == 4:
            f1s.add(float(parts[0])); f2s.add(float(parts[1])); phs.add(float(parts[2])); ro.add(parts[3])
        else:
            bad += 1
    print("bad keys:", bad)
    print("f1 vals:", sorted(f1s))
    print("f2 vals:", sorted(f2s))
    print("phases n:", len(phs), "range:", min(phs), max(phs))
    print("readoffs:", ro)
    # attrs
    print("file attrs:", dict(f.attrs) if f.attrs else "none")
    # inspect one trace
    k0 = keys[len(keys)//2]
    tr = f[k0]["trace"]
    print("example key:", k0, "trace shape:", tr.shape, tr.dtype)
    a = tr[:]
    print("trace attrs:", dict(tr.attrs))
    # group attrs
    g = f[k0]
    print("group attrs:", dict(g.attrs) if g.attrs else "none")
    print("dim0 (passes) first-vs-last diff:", np.abs(a[0,1,:]-a[-1,1,:]).mean(), "vs within-pass noise scale:", a[:,1,:].std())
    print("ch0 values (head/mid/tail):", a[0,0,:3], a[0,0,500], a[0,0,-3:])
    # check all traces same shape
    shp = set()
    for k in keys[:20]+keys[-20:]:
        shp.add(f[k]["trace"].shape)
    print("shapes seen (sample):", shp)
