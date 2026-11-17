import h5py, numpy as np
fn = 'phase gating/delay_f1_JULY_varyf1a1.h5'
with h5py.File(fn,'r') as f:
    sc = f.attrs['script_content']
    print(sc[2900:6500])
    ks = list(f.keys())
    print(len(ks), ks[:3], ks[-3:])
    g = f[ks[0]]
    print('group attrs:', dict(g.attrs))
    for k,v in g.items():
        print(' ', k, v.shape, v.dtype)
