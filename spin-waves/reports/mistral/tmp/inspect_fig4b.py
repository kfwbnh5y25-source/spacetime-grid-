import h5py, numpy as np
from collections import defaultdict

p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/detune_fine031126.h5"
with h5py.File(p, "r") as f:
    keys = sorted(f.keys())
    recs = defaultdict(list)   # (f1,f2) -> [phases]
    shapes = defaultdict(int)
    for k in keys:
        a,b,c,d = k.split("_")
        recs[(float(a),float(b))].append(float(c))
        shapes[f[k]["trace"].shape] += 1
    print("shapes:", dict(shapes))
    print("n (f1,f2) pairs:", len(recs))
    f2s = sorted(set(k[1] for k in recs))
    f1s = sorted(set(k[0] for k in recs))
    print("f2s:", f2s)
    comp = {}
    for (f1,f2),ph in recs.items():
        comp[(f1,f2)] = len(set(ph))
    # completeness per f2
    for f2 in f2s:
        f1s_here = sorted(k[0] for k in recs if k[1]==f2)
        nph = [comp[(f1_,f2)] for f1_ in f1s_here]
        print(f"f2={f2}: {len(f1s_here)} f1 pts, f1 range {min(f1s_here)}-{max(f1s_here)}, phases counts min/max {min(nph)}/{max(nph)}")
    # check freq-axis differ per key & tail behavior of spectra
    kA = "4072.00_4132.00_0.000_1"
    kB = "4060.00_4122.00_0.000_1"
    for kk in [kA,kB]:
        tr = f[kk]["trace"][:]
        print(kk, "center:", tr[0,0,500]/1e6, "MHz; spec max*1e8:", tr[0,1].max()*1e8, "mean*1e8:", tr[0,1].mean()*1e8)
    # spectrum profile of hot band point
    tr = f[kA]["trace"][:][0,1]
    fr = f[kA]["trace"][:][0,0]
    top = np.argsort(tr)[-10:]
    print("top bins freq (MHz):", fr[top]/1e6)
    print("profile (dB rel max, every 50 bins):", np.round(10*np.log10(np.maximum(tr,1e-14)/tr.max()),1)[::50])
