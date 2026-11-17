import h5py, numpy as np, re

p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
out_dir = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp"

f = h5py.File(p, 'r')
keys = list(f.keys())
pat = re.compile(r'(?P<a1>[-\d.]+)_(?P<a2>[-\d.]+)_(?P<f1>[-\d.]+)_(?P<f2>[-\d.]+)_(?P<phase>[-\d.]+)_(?P<label>\w+)')

spec = {}
spec_ref = {}
npass = {}
for k in keys:
    d = pat.match(k).groupdict()
    a1, a2, f1, f2, ph = (float(d['a1']), float(d['a2']), float(d['f1']), float(d['f2']), float(d['phase']))
    tr = f[k]['trace'][:]          # (passes, 2, 1001); [.., 1, :] = linear power
    pw = tr[:, 1, :].astype(np.float64).mean(axis=0)   # average over passes
    if a1 == 0 and f1 == 0:
        spec_ref[(a2, f2)] = pw
    else:
        spec[(a2, f2, ph)] = pw
        npass[(a2, f2, ph)] = tr.shape[0]

a1_vals = sorted(set(float(pat.match(k).groupdict()['a1']) for k in keys if float(pat.match(k).groupdict()['a1']) != 0))
print("a1 values present (nonzero):", a1_vals)

a2s = sorted(set(k[0] for k in spec))
f2s = sorted(set(k[1] for k in spec))
phs = sorted(set(k[2] for k in spec))
nmax = max(npass.values())
print("a2 range", a2s[0], "-", a2s[-1], len(a2s))
print("f2", f2s[0], "-", f2s[-1], len(f2s), "phases", phs, "max passes", nmax)

def tot(x):
    return float(np.sum(np.clip(x, 0, None)))

C = np.full((len(a2s), len(f2s)), np.nan)
Cd = np.full((len(a2s), len(f2s)), np.nan)
Sref = np.full((len(a2s), len(f2s)), np.nan)
for ia, a2 in enumerate(a2s):
    for jf, f2 in enumerate(f2s):
        k1, k2 = (a2, f2, 1.31), (a2, f2, 4.55)
        if k1 not in spec or k2 not in spec:
            continue
        r = spec_ref.get((a2, f2))
        if r is None:
            continue
        s1, s2, sr = tot(spec[k1]), tot(spec[k2]), tot(r)
        if sr <= 0:
            continue
        C[ia, jf] = 10*np.log10(s2/s1)
        Cd[ia, jf] = 10*np.log10(spec[k2].max()/max(spec[k1].max(), 1e-300))
        Sref[ia, jf] = 10*np.log10(s1/sr)

np.save(out_dir+"/contrast_fullband_avg.npy", C)
np.save(out_dir+"/contrast_peak_avg.npy", Cd)
np.save(out_dir+"/gain_p131_vs_ref_avg.npy", Sref)

idx = np.unravel_index(np.nanargmax(np.abs(C)), C.shape)
print("max |C| %.2f dB at a2=%d f2=%d" % (C[idx], a2s[idx[0]], f2s[idx[1]]))
idx = np.unravel_index(np.nanargmax(np.abs(Cd)), Cd.shape)
print("max |Cpeak| %.2f dB at a2=%d f2=%d" % (Cd[idx], a2s[idx[0]], f2s[idx[1]]))
idx = np.unravel_index(np.nanargmax(Sref), Sref.shape)
print("max gain(ph=1.31 vs f2-only ref) %.2f dB at a2=%d f2=%d" % (Sref[idx], a2s[idx[0]], f2s[idx[1]]))

print("median C (all):", np.nanmedian(C))
hb = slice(f2s.index(4120), f2s.index(4140)+1)
print("median C (f2 4120-4140):", np.nanmedian(C[:, hb]))
print("median Cd (f2 4120-4140):", np.nanmedian(Cd[:, hb]))
print("median C (a2>=240):", np.nanmedian(C[[i for i,a in enumerate(a2s) if a>=240], :]))

print("\nC (full-band dB, rows a2 250..275 step5, f2 4105..4150 step5):")
print("a2\\f2 " + " ".join("%6d" % fv for fv in range(4105, 4151, 5)))
for a2 in range(250, 276, 5):
    ia = a2s.index(float(a2))
    row = "%5d " % a2 + " ".join(
        ("%6.1f" % C[ia, f2s.index(float(fv))]) if fv in f2s and np.isfinite(C[ia, f2s.index(float(fv))]) else "   ---"
        for fv in range(4105, 4151, 5))
    print(row)

f.close()
print("\ncontrast = 10log10(P(4.55)/P(1.31)) full band, pass-averaged")
