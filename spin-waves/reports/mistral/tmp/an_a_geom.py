import h5py, numpy as np
with h5py.File('phase gating/3/3.h5','r') as f:
    ks = list(f.keys())
    import re, collections
    pats = collections.Counter()
    for k in ks:
        a,l = k.split('_')
        pats[(int(l)//10 if int(l)>=20 else 'stray')] += 1
    # key l values present
    lv = sorted({int(k.split('_')[1]) for k in ks})
    strays = [l for l in lv if l<20 or l%10!=0]
    print('l values strays (<20 or not mult of 10):', strays)
    print('count strays:', sum(1 for k in ks if int(k.split("_")[1])<20))
    g = f[ks[0]]
    tr = g['trace'][0]
    print('field0 range MHz:', tr[0].min()/1e6, tr[0].max()/1e6, 'diff', (tr[0][1]-tr[0][0])/1e6)
    print('field1 stats: mean', tr[1].mean(), 'min', tr[1].min(), 'max', tr[1].max())
    # check a few keys field consistent
    ff = []
    for k in [ks[0], ks[len(ks)//2], ks[-1]]:
        ff.append(np.allclose(f[k]['trace'][0,0], tr[0]))
    print('freq axis identical:', ff)
    # timestamp range
    ts0 = f[ks[0]]['timestamp'][0,0]; ts1 = f[ks[-1]]['timestamp'][0,0]
    print('ts first key', ts0, 'ts last key', ts1)
