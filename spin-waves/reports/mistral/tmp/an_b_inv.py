import h5py, numpy as np, collections
fn='phase gating/delay_f1_JULY_varyf1a1.h5'
recs=[]
with h5py.File(fn,'r') as f:
    for k in f:
        try:
            a,f1,ph,ro=k.split('_')
            recs.append((float(a),int(f1),float(ph),ro,k,f[k]['timestamp'][0,0].decode() if isinstance(f[k]['timestamp'][0,0],bytes) else str(f[k]['timestamp'][0,0])))
        except Exception:
            print('unparsed',k)
print('total keys',len(recs))
a1s=sorted({r[0] for r in recs}); f1s=sorted({r[1] for r in recs})
phs=sorted({r[2] for r in recs}); ros=sorted({r[3] for r in recs})
print('a1 n=%d range %s..%s'%(len(a1s),a1s[0],a1s[-1]))
print('f1 n=%d range %d..%d'%(len(f1s),f1s[0],f1s[-1]))
print('f1 step set:',sorted({f1s[i+1]-f1s[i] for i in range(len(f1s)-1)}))
print('phases n=%d'%(len(phs)),'first/last',phs[0],phs[-1], 'step',phs[1]-phs[0])
print('readoff suffixes:',collections.Counter(r[3] for r in recs))
cnt=collections.Counter((r[0],r[1]) for r in recs)
carr=np.array([cnt[k] for k in cnt])
print('phases-per-block: ', {int(v):int((carr==v).sum()) for v in sorted(set(carr))})
ts=np.array([r[5] for r in recs])
print('ts min',ts.min(),'ts max',ts.max())
a1_even=[a for a in a1s]
# check duplicates in name space
print('example block counts for (256, 4075):', cnt.get((256.0,4075)))
