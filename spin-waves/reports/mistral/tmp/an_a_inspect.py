import h5py, numpy as np
with h5py.File('phase gating/3/3.h5','r') as f:
    print(dict(f.attrs))
    ks = list(f.keys())
    print(len(ks), ks[:5], ks[-5:])
    g = f[ks[0]]
    print(g.keys() if hasattr(g,'keys') else type(g), dict(g.attrs) if hasattr(g,'attrs') else '')
    for k,v in g.items():
        print(' ', k, v.shape, v.dtype)
