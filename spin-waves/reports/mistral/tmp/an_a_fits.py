import numpy as np, pickle
from scipy.optimize import curve_fit
OUT='reports/mistral/tmp/'
amps=np.load(OUT+'A_amps.npy'); ls=np.load(OUT+'A_ls.npy')
P=np.load(OUT+'A_Pnet.npy'); Snr=np.load(OUT+'A_Snr.npy')
Fw=np.load(OUT+'A_fwhm.npy'); Sig=np.load(OUT+'A_noise.npy'); Nseg=np.load(OUT+'A_Nseg.npy')
l_ns=ls/64.0
def f_exp(l,c0,Pi,tau): return c0+Pi*(1-np.exp(-l/tau))
res=[]
for ia,a in enumerate(amps):
    y=P[ia]; s=Snr[ia]
    sig=Sig[ia]*np.sqrt(np.maximum(Nseg[ia],1))
    m=np.isfinite(y)&np.isfinite(sig)&(sig>0)
    if m.sum()<10: continue
    d=dict(a=a)
    try:
        p,pcov=curve_fit(f_exp,l_ns[m],y[m],p0=[0,max(y[m].max(),1e-20),1.0],sigma=sig[m],
                      absolute_sigma=True,maxfev=20000,bounds=([-np.inf,0,0.02],[np.inf,np.inf,50]))
        perr=np.sqrt(np.diag(pcov))
        chi=np.mean((y[m]-f_exp(l_ns[m],*p))**2/sig[m]**2)
        d.update(tau=p[2],tau_err=perr[2],Pinf=p[1],c0=p[0],chi2=chi,
                 snr_end=float(s[-1]),maxsnr=float(np.nanmax(s)),npts=int(m.sum()))
        det=s>5; start=None
        for i in range(len(det)-2):
            if det[i] and det[i+1] and det[i+2]: start=i; break
        if start is not None and (len(ls)-start)>=8:
            mm=m.copy(); mm[:start]=False
            p3,_=curve_fit(f_exp,l_ns[mm],y[mm],p0=[0,y[m].max(),1.0],sigma=sig[mm],
                           absolute_sigma=True,maxfev=20000,bounds=([-np.inf,0,0.02],[np.inf,np.inf,50]))
            d.update(tau_ons=p3[2],Pinf_ons=p3[1],lstart=l_ns[start])
    except Exception as e:
        d['failed']=str(e)[:80]
    res.append(d)
pickle.dump(res,open(OUT+'A_ringup_fits.pkl','wb'))
print('=== ring-up fits ===')
for d in res:
    if 'failed' in d: print('a=%4d FAILED %s'%(d['a'],d['failed'])); continue
    print('a=%4d tau=%6.3f+-%5.3f ns Pinf=%.3g c0=%.3g chi2=%.2f snr(9.2ns)=%.0f | onset-restricted: tau=%s (l>=%s ns)'%(
        d['a'],d['tau'],d['tau_err'],d['Pinf'],d['c0'],d['chi2'],d['snr_end'],
        ('%.3f'%d['tau_ons'] if 'tau_ons' in d else '--'),('%.2f'%d['lstart'] if 'lstart' in d else '--')))
ath=np.load(OUT+'A_ath_persist.npy'); mm=np.isfinite(ath)&(ls>=30)
a_inf=np.nanmin(ath[mm])
def thr_model(l,tau): return a_inf/np.sqrt(1-np.exp(-l/tau))
p,pc=curve_fit(thr_model,l_ns[mm],ath[mm],p0=[1.5],maxfev=20000)
tau_b,tau_be=p[0],np.sqrt(pc[0,0])
rms=np.sqrt(np.mean((ath[mm]-thr_model(l_ns[mm],*p))**2))
print('\nthreshold build-up model: a_inf=%.0f mV tau_b=%.3f +- %.3f ns rms=%.1f mV'%(a_inf,tau_b,tau_be,rms))
print(' i, l_ns, data, model:', ' | '.join('%.2f:%.0f/%.0f'%(l_ns[i],ath[i],thr_model(l_ns[i],*p)) for i in range(len(ls)) if mm[i]))
print('\nraw Pnet (1e-12 W): l=1.562 ns vs 9.219 ns:')
for ia,a in enumerate(amps):
    if a>560: break
    print('  a=%3d: P(1.56ns)=%.3g P(9.22ns)=%.3g snr=(%.1f, %.1f)'%(a,P[ia,8],P[ia,-1],Snr[ia,8],Snr[ia,-1]))
