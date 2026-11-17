import h5py, numpy as np
from collections import defaultdict
p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/detune_fine031126.h5"
with h5py.File(p, "r") as f:
    G = defaultdict(dict)
    for k in f.keys():
        a,b,c,d = k.split("_")
        tr = f[k]["trace"][:][0]
        G[(float(a),float(b))][float(c)] = tr[1].astype(np.float64).mean()
def fit(ph,I):
    X = np.column_stack([np.ones_like(ph),np.cos(ph),np.sin(ph),np.cos(2*ph),np.sin(2*ph)])
    c,*_ = np.linalg.lstsq(X,I,rcond=None)
    pred=X@c; ssr=((I-pred)**2).sum(); sst=((I-I.mean())**2).sum()
    A1c,A1s,A2c,A2s=c[1],c[2],c[3],c[4]
    return dict(A1=np.hypot(A1c,A1s),A2=np.hypot(A2c,A2s),phi0=np.degrees(np.arctan2(-A1s,A1c)),
                phi2=np.degrees(np.arctan2(-A2s,A2c)),r2=1-ssr/sst)
res={k:fit(np.array(sorted(d)),np.array([d[pp] for pp in sorted(d)])) for k,d in G.items()}
print("phi0 (deg) map:")
for f2 in (4122.0,4125.0,4132.0,4135.0,4139.0,4142.0,4148.0):
    row=[]
    for f1 in range(4060,4085,2):
        if (float(f1),f2) in res:
            r=res[(float(f1),f2)]
            row.append(f"{f1}:{r['phi0']:+6.0f}")
    print(f"f2={f2}:", " ".join(row))
print()
print("A2/A1 ratio map:")
for f2 in (4122.0,4125.0,4132.0,4135.0,4139.0,4142.0,4148.0):
    row=[]
    for f1 in range(4060,4085,2):
        if (float(f1),f2) in res:
            r=res[(float(f1),f2)]
            row.append(f"{f1}:{r['A2']/r['A1']:.2f}")
    print(f"f2={f2}:", " ".join(row))
