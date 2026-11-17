import h5py, numpy as np, collections, pickle
fn='phase gating/delay_f1_JULY_varyf1a1.h5'
OUT='reports/mistral/tmp/'
F2=4150.0
with h5py.File(fn,'r') as f:
    allk=[k for k in f]
a1s=sorted({float(k.split('_')[0]) for k in allk})
f1s=sorted({int(k.split('_')[1]) for k in allk})
NA,NF=len(a1s),len(f1s)
ai={v:i for i,v in enumerate(a1s)}; fi={v:i for i,v in enumerate(f1s)}
A1=np.full((NA,NF),np.nan); A2=np.full_like(A1,np.nan); R21=np.full_like(A1,np.nan)
A0=np.full_like(A1,np.nan); ARG2=np.full_like(A1,np.nan); PHI0=np.full_like(A1,np.nan)
CORR=np.zeros((NA,NF),bool); NP=np.zeros((NA,NF),int)
blockdatalist=[]
SKIPPED=0
with h5py.File(fn,'r') as f:
    blocks=collections.defaultdict(dict)
    for k in allk:
        a,f1,ph,ro=k.split('_')
        blocks[(float(a),int(f1))][float(ph)]=k
    freq_axis=None
    for (a,f1),pd in sorted(blocks.items()):
        phs=sorted(pd); ph=np.array(phs,dtype=float); n=len(ph)
        if n<9:
            SKIPPED+=n; continue
        dphi=np.diff(ph); dp=np.median(dphi)
        if np.max(np.abs(dphi-dp))>0.02 or dp>0.5:
            SKIPPED+=n; continue
        I=np.zeros(n); O=np.zeros(n); T=[]; keys=[]
        for j,p in enumerate(phs):
            k=pd[p]
            tr=f[k]['trace'][0]
            if freq_axis is None:
                freq_axis=tr[0].copy(); fr=freq_axis/1e6
            s=tr[1].astype(np.float64)
            lw=np.abs(fr-2075.0)<=1.75
            I[j]=s[lw].sum(); O[j]=np.median(s[~lw])
            ts=f[k]['timestamp'][0,0]
            T.append(ts.decode() if isinstance(ts,bytes) else str(ts)); keys.append(k)
        order=np.argsort(T); Iq=I[order]; Oq=O[order]
        corrupt=False
        if n>=8:
            best=None
            for cut in range(3,n-2):
                s1=Oq[:cut]-Oq[:cut].mean(); s2=Oq[cut:]-Oq[cut:].mean()
                v=(np.sum(s1**2)/max(cut-1,1)+np.sum(s2**2)/max(n-cut-1,1))
                dlev=Oq[cut:].mean()-Oq[:cut].mean()
                se=np.sqrt(v/cut+v/(n-cut))
                t=dlev/se if se>0 else 0
                if best is None or abs(t)>abs(best[0]): best=(abs(t),cut,dlev)
            if best and best[0]>6: corrupt=True
        e=np.exp(-1j*ph)
        c1=2*np.sum(I*e)/n
        c2=2*np.sum(I*np.exp(-2j*ph))/n
        ia,ib=ai[a],fi[f1]
        A1[ia,ib]=np.abs(c1); A2[ia,ib]=np.abs(c2)
        R21[ia,ib]=abs(c2/c1) if abs(c1)>0 else np.nan
        A0[ia,ib]=I.mean(); ARG2[ia,ib]=np.angle(c2); PHI0[ia,ib]=np.angle(c1)
        CORR[ia,ib]=corrupt; NP[ia,ib]=n
        blockdatalist.append(dict(a1=a,f1=f1,n=n,phi=ph,I=I,off=O,ts=np.array(T),
                                  keys=np.array(keys),corrupt=corrupt))
np.save(OUT+'B_A1.npy',A1); np.save(OUT+'B_A2.npy',A2); np.save(OUT+'B_R21.npy',R21)
np.save(OUT+'B_A0.npy',A0); np.save(OUT+'B_ARG2.npy',ARG2); np.save(OUT+'B_PHI0.npy',PHI0)
np.save(OUT+'B_CORRUPT.npy',CORR); np.save(OUT+'B_NPH.npy',NP)
np.save(OUT+'B_a1s.npy',np.array(a1s)); np.save(OUT+'B_f1s.npy',np.array(f1s,dtype=int))
pickle.dump(blockdatalist,open(OUT+'B_blocks.pkl','wb'))
print('blocks used:',len(blockdatalist),'skipped keys:',SKIPPED)
print('corrupt blocks:',int(CORR.sum()),'of',int(np.isfinite(A1).sum()))
ok=~CORR&np.isfinite(A1)
ix=np.unravel_index(np.nanargmax(np.where(ok,A1,np.nan)),A1.shape)
print('A1 max %.4g at (a1=%.0f, f1=%d, Delta=%.1f MHz); min %.4g'%(A1[ix], a1s[ix[0]], f1s[ix[1]], F2-f1s[ix[1]], np.nanmin(A1[ok])))
print('A0 max %.4g min %.4g'%(np.nanmax(A0[ok]),np.nanmin(A0[ok])))
r=np.where(ok,R21,np.nan)
print('A2/A1: max %.3f min %.3f median %.3f'%(np.nanmax(r),np.nanmin(r),np.nanmedian(r)))
g2=ARG2[ok]
print('ARG2 over ok: frac |arg(c2)|>pi/2: %.3f (%d pts)'%(np.mean(np.abs(g2)>np.pi/2),len(g2)))
