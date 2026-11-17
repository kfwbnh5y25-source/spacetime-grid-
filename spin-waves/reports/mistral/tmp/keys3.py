import h5py
with h5py.File("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/3/3.h5", "r") as f:
    ks = sorted(f.keys())
    print(len(ks), ks[:15], ks[-8:])
    import collections
    print(collections.Counter(len(k.split("_")) for k in ks))
