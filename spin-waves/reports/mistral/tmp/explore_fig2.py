import h5py, numpy as np
f = h5py.File('/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5','r')
for k in ['0_200_0_4075_0.000_1','0_200_0_4111_0.000_1','0_200_0_4150_0.000_1']:
    tr = f[k]['trace'][:]
    fr = tr[0,0,:]
    print(k, 'center %.4f MHz'%(fr.mean()*1e-3), 'span %.3f MHz'%((fr[-1]-fr[0])*1e-3))

def load(k):
    tr = f[k]['trace'][:].astype(np.float64)
    return tr[:,0,:].mean(0), tr[:,1,:].mean(0)

fr, p311 = load('324_275_4069_4111_1.310_1')
_, p455 = load('324_275_4069_4111_4.550_1')
_, pref  = load('0_275_0_4111_0.000_1')
print('freq edges MHz:', fr[0]*1e-3, fr[-1]*1e-3)
print('Ptot clip0 sum: phase1.31 %.4g  phase4.55 %.4g  ref %.4g'%(p311.clip(0).sum(), p455.clip(0).sum(), pref.clip(0).sum()))
print('peak max: p1.31 %.4g at %.3f MHz | p4.55 %.4g at %.3f MHz | ref %.4g at %.3f MHz'%(
    p311.max(), fr[p311.argmax()]*1e-3, p455.max(), fr[p455.argmax()]*1e-3, pref.max(), fr[pref.argmax()]*1e-3))
d = p311 - pref
print('top |diff| points (phase1.31 - ref):')
for i in np.argsort(np.abs(d))[-10:]:
    print('  %+.3g at %.3f MHz'%(d[i], fr[i]*1e-3))
print('spectrum overview phase1.31 (every 50th pt, uW-ish linear):')
print(np.array2string(p311[::100], precision=3))
print('ref overview:')
print(np.array2string(pref[::100], precision=3))
