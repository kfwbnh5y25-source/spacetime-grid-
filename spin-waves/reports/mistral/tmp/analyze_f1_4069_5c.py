import h5py, numpy as np, re

p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
out_dir = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp"

f = h5py.File(p, 'r')
keys = list(f.keys())
pat = re.compile(r'(?P<a1>[-\d.]+)_(?P<a2>[-\d.]+)_(?P<f1>[-\d.]+)_(?P<f2>[-\d.]+)_(?P<phase>[-\d.]+)_(?P<label>\w+)')

spec, spec_ref = {}, {}
for k in keys:
    d = pat.match(k).groupdict()
    a1, a2, f1, f2, ph = (float(d['a1']), float(d['a2']), float(d['f1']), float(d['f2']), float(d['phase']))
    pw = f[k]['trace'][:, 1, :].astype(np.float64).mean(axis=0)
    if a1 == 0 and f1 == 0:
        spec_ref[(a2, f2)] = pw
    else:
        spec[(a1, a2, f1, f2, ph)] = pw

combos = sorted(set((a1, f1) for (a1, a2, f1, f2, ph) in spec))
print("(a1, f1) combos:", combos)

def tot(x):
    return float(np.sum(np.clip(x, 0, None)))

for a1, f1 in combos:
    sel = [k for k in spec if k[0] == a1 and k[2] == f1]
    a2s = sorted(set(k[1] for k in sel))
    f2s = sorted(set(k[3] for k in sel))
    C = np.full((len(a2s), len(f2s)), np.nan)
    G = np.full((len(a2s), len(f2s)), np.nan)
    for ia, a2 in enumerate(a2s):
        for jf, f2 in enumerate(f2s):
            k1, k2 = (a1, a2, f1, f2, 1.31), (a1, a2, f1, f2, 4.55)
            if k1 in spec and k2 in spec:
                C[ia, jf] = 10*np.log10(tot(spec[k2])/tot(spec[k1]))
            r = spec_ref.get((a2, f2))
            if r is not None and k2 in spec:
                G[ia, jf] = 10*np.log10(tot(spec[k2])/tot(r))
    np.save(out_dir+f"/f1grp_a{int(a1)}_contrast.npy", C)
    np.save(out_dir+f"/f1grp_a{int(a1)}_gain_vs_ref.npy", G)
    fin = np.isfinite(C)
    print(f"\na1={a1:.0f} f1={f1:.0f}: {fin.sum()} pts, a2 {a2s[0]:.0f}-{a2s[-1]:.0f} ({len(a2s)}), "
          f"max|C| {np.nanmax(np.abs(C)):.2f} dB, median C {np.nanmedian(C):+.2f}, "
          f"max gain {np.nanmax(G):.1f} dB")
    # top-5 contrast rows
    Cf = np.abs(C)
    if fin.any():
        idx = np.unravel_index(np.nanargmax(Cf), Cf.shape)
        print("  max|C| at a2=%d f2=%d (C=%+.2f dB); gain there %+.1f dB" %
              (a2s[idx[0]], f2s[idx[1]], C[idx], G[idx]))
f.close()
