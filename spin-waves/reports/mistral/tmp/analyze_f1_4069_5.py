import h5py, numpy as np, re, json

p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
out_dir = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp"

f = h5py.File(p, 'r')
keys = list(f.keys())
pat = re.compile(r'(?P<a1>[-\d.]+)_(?P<a2>[-\d.]+)_(?P<f1>[-\d.]+)_(?P<f2>[-\d.]+)_(?P<phase>[-\d.]+)_(?P<label>\w+)')

data = {}
ref = {}
for k in keys:
    m = pat.match(k)
    d = m.groupdict()
    a1, a2, f1, f2, ph = (float(d['a1']), float(d['a2']), float(d['f1']), float(d['f2']), float(d['phase']))
    tr = f[k]['trace'][:]  # (2, 2, 1001)
    if a1 == 0 and f1 == 0:
        ref[(a2, f2)] = tr
    else:
        data[(a2, f2)] = tr  # dict key a2/f2 unique since a1,f1 fixed

a2s = sorted(set(k[0] for k in data))
f2s = sorted(set(k[1] for k in data))
print("a2 grid", a2s)
print("f2 grid", f2s[0], "..", f2s[-1], len(f2s))

# trace shape (2, 2, 1001): first dim = phase index (0: ph=1.310?), second = ? (maybe 2 repeated traces?), third = SA points
tr0 = data[(a2s[0], f2s[0])]
print("trace shape", tr0.shape)

# Check second dim: correlation/means
for i in range(2):
    for j in range(2):
        seg = np.power(10.0, tr0[i, j].astype(np.float64)/10.0) if tr0.min() < 50 else tr0[i, j]
        print(f"dim[{i},{j}] min {tr0[i,j].min():.1f} max {tr0[i,j].max():.1f} mean {tr0[i,j].mean():.1f}")

# Assume trace[p, j, :] : p indexes the 2 phases; j maybe two consecutive SA sweeps -> average them
def power_db(tr):
    # returns per-phase mean power (linear) averaged over sweeps
    lin = tr.astype(np.float64)
    if lin.max() < 100:  # dBm stored
        lin = np.power(10.0, lin/10.0)
    return lin.mean(axis=1)  # (2, 1001)

# Contrast 4.55 vs 1.31 relative to reference
# Need phase order: check against script phases [4.55, 1.31]; try both orderings by correlation with ref
# Reference trace: ref[(a2,f2)] shape (2,2,1001) as well (phase 0 only -> likely both slots same)
res = np.full((len(a2s), len(f2s)), np.nan)
res_abs = np.full((len(a2s), len(f2s), 2), np.nan)
for ia, a2 in enumerate(a2s):
    for jf, f2 in enumerate(f2s):
        if (a2, f2) not in data:
            continue
        tr = data[(a2, f2)]
        r = ref.get((a2, f2))
        P = power_db(tr)  # (2, 1001) linear, per phase
        if r is not None:
            Pr = power_db(r)
            # normalize: subtract mean reference power (scalar)
            Pref = Pr.mean()
        else:
            Pref = P.min() * 0 + np.nan
        # phase index guess: find which slot is closer to reference (phase closer to ideal 0 would be arbitrary)
        # instead compute contrast = 10log10(P[0].sum()/P[1].sum()) and record; ordering verdict separate
        s = P.sum(axis=1)
        res[ia, jf] = 10*np.log10(s[0]/s[1])
        res_abs[ia, jf, :] = s / 1001  # mean linear power per phase

order_guess = {}
# Determine phase ordering: response at phase 1.31 vs 4.55. Use a2=275 (max amp) f2 where interaction strong.
np.save(out_dir+"/contrast_a2xf2.npy", res)
np.save(out_dir+"/absmean_a2xf2xph.npy", res_abs)

# summary: best |contrast| tuples
idx = np.unravel_index(np.nanargmax(np.abs(res)), res.shape)
print("max |contrast| %.2f dB at a2=%d f2=%d" % (res[idx], a2s[idx[0]], f2s[idx[1]]))
idxn = np.unravel_index(np.nanargmin(res), res.shape)
print("min contrast %.2f dB at a2=%d f2=%d" % (res[idxn], a2s[idxn[0]], f2s[idxn[1]]))

# Contrast averaged over f2 for hot band rows
hotband_f2 = np.array([4110, 4115, 4120, 4125, 4130])
for a2 in a2s[-6:]:
    row = [f"{a2}mV:"] + ["%.1f" % res[a2s.index(a2), f2s.index(fv)] for fv in hotband_f2 if fv in f2s]
    print(" ".join(row))

# absolute on/off style: compare each phase vs reference
print("\nPhase responses vs reference (mean over trace, dB rel ref):")
for ph_i in range(2):
    mat = np.full((1,1), np.nan)
    vals = []
    for ia, a2 in enumerate(a2s):
        for jf, f2 in enumerate(f2s):
            if (a2, f2) in ref and not np.isnan(res_abs[ia,jf,ph_i]):
                r = power_db(ref[(a2,f2)]).mean()
                vals.append(10*np.log10(res_abs[ia,jf,ph_i]/r))
    vals = np.array(vals)
    print(f"ph slot {ph_i}: median {np.nanmedian(vals):.2f} dB rel ref, p95 {np.nanpercentile(vals,95):.2f}, max {np.nanmax(vals):.2f}")

f.close()
print("saved npy to", out_dir)
