import h5py, numpy as np, time

H5 = '/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5'
t0 = time.time()
f = h5py.File(H5, 'r')
ks = np.array(sorted(f.keys()))
a2s = np.arange(200, 276)
f2s = np.arange(4075, 4151)

def metrics(k):
    tr = f[k]['trace'][:].astype(np.float64)
    p = tr[:, 1, :].mean(0)
    pc = np.maximum(p, 0)
    sm = np.convolve(pc, np.ones(5)/5, mode='same')
    return dict(total=float(pc.sum()), peak=float(sm.max()),
                rawpeak=float(pc.max()))

# references a1=0 f1=0
ref = {}
for k in ks:
    if k.startswith('0_'):
        p = k.split('_')
        ref[(int(p[1]), int(p[3]))] = metrics(k)
print('ref loaded', len(ref), time.time()-t0)

# a1=324 both phases
dat = {}
for k in ks:
    if k.startswith('324_'):
        p = k.split('_')
        dat[(int(p[1]), int(p[3]), p[4])] = metrics(k)
print('324 loaded', len(dat), time.time()-t0)

gmap = np.full((len(a2s), len(f2s)), np.nan)   # total gain phase 1.31 / ref
cmap = np.full((len(a2s), len(f2s)), np.nan)   # 10log10 total(4.55)/total(1.31)
gpk  = np.full((len(a2s), len(f2s)), np.nan)   # peak gain
cpk  = np.full((len(a2s), len(f2s)), np.nan)
for i, a2 in enumerate(a2s):
    for j, f2 in enumerate(f2s):
        r = ref.get((a2, f2)); a = dat.get((a2, f2, '1.310')); b = dat.get((a2, f2, '4.550'))
        if r and a and b:
            gmap[i, j] = a['total']/r['total']
            cmap[i, j] = 10*np.log10(b['total']/a['total'])
            gpk[i, j] = a['peak']/r['peak']
            cpk[i, j] = 10*np.log10(b['peak']/a['peak'])

for name, m in [('gain_total', gmap), ('contrast_total', cmap), ('gain_peak', gpk), ('contrast_peak', cpk)]:
    print(name, 'min %.3f max %.3f  argmax (a2=%d, f2=%d)' % (
        np.nanmin(m), np.nanmax(m), a2s[np.unravel_index(np.nanargmax(m), m.shape)[0]],
        f2s[np.unravel_index(np.nanargmax(m), m.shape)[1]]))

# gain vs a2 at a few f2
for f2 in [4095, 4111, 4120]:
    j = list(f2s).index(f2)
    row = np.array([dat[(int(a2), f2, ph)]['total'] for a2 in a2s for ph in ['1.310','4.550']]).reshape(-1,2)
    rowr = np.array([ref[(int(a2), f2)]['total'] for a2 in a2s])
    print('\nf2=%d  gain(1.31) dB: ' % f2, np.array2string(10*np.log10(row[:,0]/rowr), precision=1))
np.save('/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/gmap.npy', gmap)
np.save('/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/cmap.npy', cmap)
print('elapsed', time.time()-t0)
