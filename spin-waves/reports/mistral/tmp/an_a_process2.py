import h5py, numpy as np
from scipy.signal import medfilt

OUT='reports/mistral/tmp/'
amps=np.arange(200,1001,10); ls=np.arange(20,600,10)
nf_a, nf_l = amps.size, ls.size
Pwin=np.full((nf_a,nf_l),np.nan); Nseg=np.zeros((nf_a,nf_l),int)
Pkpos=np.full_like(Pwin,np.nan); Fwhm=np.full_like(Pwin,np.nan)
Bg=np.full_like(Pwin,np.nan); Sig=np.full_like(Pwin,np.nan)
Ptot=np.full_like(Pwin,np.nan); freqs=None
with h5py.File('phase gating/3/3.h5','r') as f:
    for key in f:
        try: a,l=map(int,key.split('_'))
        except ValueError: continue
        if l<20 or l%10!=0: continue
        ia=np.where(amps==a)[0]; il=np.where(ls==l)[0]
        if ia.size!=1 or il.size!=1: continue
        ia,il=ia[0],il[0]
        tr=f[key]['trace'][0]
        if freqs is None: freqs=tr[0].copy()
        s=tr[1].astype(np.float64)
        fm=medfilt(s,15)
        bg=np.median(fm); pk=np.argmax(fm); h=fm[pk]-bg
        thr4=bg+0.25*h; thr2=bg+0.5*h
        lo=pk
        while lo>0 and fm[lo-1]>thr4: lo-=1
        hi=pk
        while hi<len(fm)-1 and fm[hi+1]>thr4: hi+=1
        lw=pk
        while lw>0 and fm[lw-1]>thr2: lw-=1
        hw=pk
        while hw<len(fm)-1 and fm[hw+1]>thr2: hw+=1
        seg=s[lo:hi+1]
        Pwin[ia,il]=seg.sum(); Nseg[ia,il]=len(seg)
        Pkpos[ia,il]=freqs[pk]/1e6; Fwhm[ia,il]=(hw-lw)*0.075
        Bg[ia,il]=bg
        std=1.4826*np.median(np.abs(s-np.median(s))); Sig[ia,il]=std
        Ptot[ia,il]=s.sum()
Pnet=Pwin-Bg*Nseg
Snr=Pnet/(Sig*np.sqrt(np.maximum(Nseg,1)))
np.save(OUT+'A_Pnet.npy',Pnet); np.save(OUT+'A_Nseg.npy',Nseg)
np.save(OUT+'A_Pkpos.npy',Pkpos); np.save(OUT+'A_fwhm.npy',Fwhm)
np.save(OUT+'A_Bg.npy',Bg); np.save(OUT+'A_noise.npy',Sig)
np.save(OUT+'A_Snr.npy',Snr); np.save(OUT+'A_Ptot.npy',Ptot)
print('threshold from NET window power, 5-sigma, persistence>=2 amps:')
det=Snr>5
a_th=np.full(nf_l,np.nan); a_th_strict=np.full(nf_l,np.nan)
for il in range(nf_l):
    d=det[:,il]
    idx=np.where(d)[0]
    if idx.size: a_th[il]=amps[idx[0]]
    # persistence: need 3 consecutive detections
    run=0; first=None
    for i in range(nf_a):
        run=run+1 if det[i,il] else 0
        if run>=3: first=amps[i-2]; break
    a_th_strict[il]=first if first else np.nan
for label,arr in [('naive',a_th),('persist3',a_th_strict)]:
    print(label, ' '.join(('%.0f'%v if np.isfinite(v) else '--') for v in arr))
np.save(OUT+'A_ath_naive.npy',a_th); np.save(OUT+'A_ath_persist.npy',a_th_strict)
