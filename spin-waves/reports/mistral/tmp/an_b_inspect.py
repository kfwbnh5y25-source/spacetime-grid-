import h5py, numpy as np, json
fn = 'phase gating/delay_f1_JULY_varyf1a1.h5'
with h5py.File(fn,'r') as f:
    print('attrs:', list(f.attrs.keys()))
    for k,v in f.attrs.items():
        s = str(v)
        print('---', k, len(s))
        print(s[:3000])
    ks = list(f.keys())
    print(len(ks), ks[:3], ks[-3:])
    g = f[ks[0]]
    print(dict(g.attrs))
    for k,v in g.items():
        print(' ', k, v.shape, v.dtype, np.array(v).ravel()[:4] if v.size<10 else '')
