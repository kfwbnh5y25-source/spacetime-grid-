import h5py, numpy as np, collections
fn='phase gating/delay_f1_JULY_varyf1a1.h5'
blocks=collections.defaultdict(list)
with h5py.File(fn,'r') as f:
    for k in f:
        a,f1,ph,ro=k.split('_')
        blocks[(float(a),int(f1))].append(ph)
g16=g32=0
gridsets=collections.Counter()
for kk,v in blocks.items():
    p=np.array(sorted(v),dtype=float); d=np.diff(np.round(p,3))
    gridsets[(len(v), round(float(np.median(d)),3))]+=1
print('block (nphases, median dphase):', dict(gridsets))
with h5py.File(fn,'r') as f:
    k0=[k for k in f if k.startswith('328.0000_4118')][0]
    tr=f[k0]['trace'][0]
    s=tr[1].astype(np.float64); fr=tr[0]/1e6
peak=collections.Counter()
with h5py.File(fn,'r') as f:
    kk=[k for k in f if k.split('_')[1]=='4075']
    print('n keys at f1=4075:', len(kk))
    acc=np.zeros(1001); n=0
    big=np.zeros(1001); nb=0
    for k in kk:
        s=f[k]['trace'][0,1].astype(np.float64)
        a=float(k.split('_')[0])
        if a>320: big+=s; nb+=1
        acc+=s; n+=1
    fm=np.convolve(acc/n,np.ones(9)/9,'same')
    print('mean over f1=4075 keys: peak at %.3f MHz, val %.3g, background median %.3g'%(fr[np.argmax(fm)],fm.max(),np.median(fm)))
    fmb=np.convolve(big/nb,np.ones(9)/9,'same')
    print('a1>320 subset: peak %.3f MHz val %.3g bg %.3g'%(fr[np.argmax(fmb)],fmb.max(),np.median(fmb)))
    print('profile +-3MHz:', np.array2string(fmb[(np.abs(fr-2075)<3)]*1e12, precision=1))
