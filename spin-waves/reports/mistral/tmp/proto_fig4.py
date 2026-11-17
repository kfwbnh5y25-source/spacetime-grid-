# Exploratory prototype (tmp/) for Fig 4 phase-Fourier analysis; not a figure script.
import h5py, numpy as np
from collections import defaultdict

p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/detune_fine031126.h5"

def fourier_fit(ph, I):
    X = np.column_stack([np.ones_like(ph), np.cos(ph), np.sin(ph), np.cos(2*ph), np.sin(2*ph)])
    c, *_ = np.linalg.lstsq(X, I, rcond=None)
    pred = X @ c
    ss_res = np.sum((I-pred)**2); ss_tot = np.sum((I-I.mean())**2)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
    return dict(A0=c[0], A1c=c[1], A1s=c[2], A2c=c[3], A2s=c[4],
                A1=np.hypot(c[1], c[2]), A2=np.hypot(c[3], c[4]), r2=r2,
                A1rel=np.hypot(c[1], c[2])/max(abs(c[0]), 1e-30))

with h5py.File(p, "r") as f:
    groups = defaultdict(dict)
    for k in f.keys():
        a,b,c,d = k.split("_")
        tr = f[k]["trace"][:][0]
        groups[(float(a),float(b))][float(c)] = tr[1].astype(np.float64).mean()

results = {}
for (f1,f2), pd_ in sorted(groups.items()):
    ph = np.array(sorted(pd_))
    I = np.array([pd_[pp] for pp in ph])
    results[(f1,f2)] = fourier_fit(ph, I)

f1s = sorted(set(k[0] for k in results)); f2s = sorted(set(k[1] for k in results))
A1 = np.full((len(f2s), len(f1s)), np.nan); A2 = np.full_like(A1, np.nan)
R2 = np.full_like(A1, np.nan)
for (f1,f2), r in results.items():
    A1[f2s.index(f2), f1s.index(f1)] = r["A1"]
    A2[f2s.index(f2), f1s.index(f1)] = r["A2"]
    R2[f2s.index(f2), f1s.index(f1)] = r["r2"]

print("A1: min %.3e max %.3e median %.3e" % (np.nanmin(A1), np.nanmax(A1), np.nanmedian(A1)))
print("A2 max %.3e | R2: median %.3f max %.3f, frac>0.8 %.2f" % (np.nanmax(A2), np.nanmedian(R2), np.nanmax(R2), np.nanmean(R2>0.8)))
amax = np.unravel_index(np.nanargmax(A1), A1.shape)
print("A1 max at f1=%d f2=%d A1=%.3e R2=%.3f Delta=%d Sigma=%.0f" %
      (f1s[amax[1]], f2s[amax[0]], A1[amax], R2[amax], f2s[amax[0]]-f1s[amax[1]], f2s[amax[0]]+f1s[amax[1]]))
for f2 in f2s:
    row = [results[(f1,f2)] for f1 in f1s if (f1,f2) in results]
    print(f"f2={f2}: n={len(row)} medianA1={np.nanmedian([r['A1'] for r in row]):.3e} maxA1={np.nanmax([r['A1'] for r in row]):.3e} maxR2={np.nanmax([r['r2'] for r in row]):.3f}")
offmed = np.nanmedian([results[k]['A1'] for k in results if k[1] in (4139.0,4142.0,4143.0,4148.0)])
print("A1 contrast max/off-band-median: %.1f dB (off-median=%.3e)" % (10*np.log10(np.nanmax(A1)/offmed), offmed))
print("--- hot band (f2=4132,4135, f1 4066..4078) phase-sign check ---")
for (f1,f2) in sorted([(a,b) for (a,b) in results if b in (4132.0,4135.0) and 4066<=a<=4078]):
    r = results[(f1,f2)]
    print(f"f1={f1} f2={f2}: A1c={r['A1c']:+.3e} A1rel={r['A1rel']:.3f} R2={r['r2']:.3f}")
print("--- low band (f2=4122,4125) f1 4070..4074 for comparison ---")
for (f1,f2) in sorted([(a,b) for (a,b) in results if b in (4122.0,4125.0) and 4070<=a<=4074]):
    r = results[(f1,f2)]
    print(f"f1={f1} f2={f2}: A1c={r['A1c']:+.3e} A1rel={r['A1rel']:.3f} R2={r['r2']:.3f}")
