import h5py, collections
with h5py.File("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/3/3.h5", "r") as f:
    ks = sorted(f.keys())
    sec = collections.Counter(int(k.split("_")[1]) for k in ks)
    fir = sorted(set(int(k.split("_")[0]) for k in ks))
    print("distinct first-field (amplitude?):", len(fir), fir[:5], fir[-3:])
    weird = [k for k in ks if int(k.split("_")[1]) < 20]
    print("keys with second field < 20:", weird[:20], len(weird))
    print("distinct second-field l:", len(sec), sorted(sec.keys())[:10], sorted(sec.keys())[-5:])
