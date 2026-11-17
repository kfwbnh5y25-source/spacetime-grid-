import numpy as np, pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT='reports/mistral/tmp/'
amps=np.load(OUT+'A_amps.npy'); ls=np.load(OUT+'A_ls.npy'); l_ns=ls/64.0
P=np.load(OUT+'A_Pnet.npy'); Fw=np.load(OUT+'A_fwhm.npy')
ath=np.load(OUT+'A_ath_persist.npy'); res=pickle.load(open(OUT+'A_ringup_fits.pkl','rb'))
fig,ax=plt.subplots(2,2,figsize=(10,8))
im=ax[0,0].pcolormesh(l_ns,amps,np.log10(np.abs(P)+1e-13),shading='auto',cmap='viridis')
ax[0,0].plot(l_ns,ath,'r*',ms=4)
plt.colorbar(im,ax=ax[0,0])
ax[0,0].set_xlabel('pulse length (ns)'); ax[0,0].set_ylabel('amplitude (mV)')
ax[0,0].set_title('A: $P_{f/2}$ window power (log$_{10}$ W), star = $a_{th}(l)$')
a_i=np.nanmin(ath[np.isfinite(ath)&(ls>=30)])
tt=np.linspace(0.3,9.4,200)
ax[0,1].plot(l_ns,ath,'bo',label='data')
ax[0,1].plot(tt,a_i/np.sqrt(1-np.exp(-tt/2.153)),'r-',label='build-up $\\tau_b$=2.15 ns')
ax[0,1].axhline(a_i,color='gray',ls='--',label='$a_\\infty=320$ mV')
ax[0,1].legend(fontsize=8); ax[0,1].set_title('A: threshold vs pulse length')
ax[1,0].plot([d['a'] for d in res if 'tau' in d and 420<=d['a']<=900],[d['tau'] for d in res if 'tau' in d and 420<=d['a']<=900],'o-',label='full range')
ax[1,0].plot([d['a'] for d in res if 'tau_ons' in d and d['a']<=560],[d['tau_ons'] for d in res if 'tau_ons' in d and d['a']<=560],'s-',label='from onset')
ax[1,0].set_ylim(0,10); ax[1,0].set_xlabel('a (mV)'); ax[1,0].set_ylabel(r'$\tau_r$ (ns)'); ax[1,0].legend(fontsize=8)
ax[1,0].set_title('A: ring-up time vs amplitude')
for a,c in [(400,'C0'),(600,'C1'),(1000,'C2')]:
    ia=np.where(amps==a)[0][0]
    ax[1,1].plot(l_ns,Fw[ia],'o-',color=c,label='a=%d mV'%a)
ax[1,1].set_xlabel('pulse length (ns)'); ax[1,1].set_ylabel('FWHM (MHz)'); ax[1,1].legend(fontsize=8)
ax[1,1].set_title('A: line width vs pulse length')
plt.tight_layout(); plt.savefig(OUT+'fig_A_summary.png',dpi=110); plt.close()

a1s=np.load(OUT+'B_a1s.npy'); f1s=np.load(OUT+'B_f1s.npy')
A1=np.load(OUT+'B_A1.npy'); A2=np.load(OUT+'B_A2.npy'); R21=np.load(OUT+'B_R21.npy')
A0=np.load(OUT+'B_A0.npy'); ARG2=np.load(OUT+'B_ARG2.npy'); CORR=np.load(OUT+'B_CORRUPT.npy')
Ok=np.isfinite(A1); Use=~CORR&Ok
def msk(Z,log=False):
    Z=Z.copy().astype(float); Z[~Use]=np.nan
    if log: Z=np.log10(np.where(Z>0,Z,np.nan))
    return Z
D=4150-f1s
fig,ax=plt.subplots(2,2,figsize=(14,9))
panels=[(msk(A1,log=True),r'log$_{10}$ A1 (W bin-sum, cos$\phi$ amp)','viridis'),
        (msk(R21),'A2/A1','viridis'),
        (msk(A0,log=True),r'log$_{10}$ A0 (phase-mean power)','viridis'),
        (np.angle(np.exp(1j*msk(ARG2))),r'arg$(c_2)$ wrapped','twilight')]
for k,(Z,t,c) in enumerate(panels):
    im=ax.flat[k].pcolormesh(D,a1s,Z,shading='auto',cmap=c)
    plt.colorbar(im,ax=ax.flat[k])
    ax.flat[k].axvline(74.5,color='r',ls='--',lw=.8)
    ax.flat[k].set_xlabel(r'$\Delta=4150-f_1$ (MHz)'); ax.flat[k].set_ylabel('a1 (mV)'); ax.flat[k].set_title(t)
plt.tight_layout(); plt.savefig(OUT+'fig_B_maps.png',dpi=110); plt.close()
fig,ax=plt.subplots(1,3,figsize=(15,4))
Av=np.nanmean(np.where(Use,A1,np.nan),axis=0); Avs=np.nanstd(np.where(Use,A1,np.nan),axis=0)/np.sqrt(np.maximum(np.sum(Use,0),1))
ax[0].errorbar(D,Av,yerr=Avs,fmt='o-',ms=3,label='mean over a1')
i3=np.where(a1s==328)[0][0]
ax[0].plot(D,np.where(Use[i3],A1[i3],np.nan),'^-',ms=3,label='a1=328')
ax[0].axvspan(50,58,color='g',alpha=.15,label=r'$\Delta$ 50-58')
ax[0].axvline(74.5,color='r',ls='--')
ax[0].set_xlabel(r'$\Delta$ (MHz)'); ax[0].set_ylabel('A1 (W)'); ax[0].legend(fontsize=8); ax[0].set_title('B: A1 vs $\\Delta$')
gm=np.nanmedian(np.where(Use,np.angle(np.exp(1j*ARG2)),np.nan),axis=0)
ax[1].plot(D,gm/np.pi*180,'o-',ms=3); ax[1].axvline(74.5,color='r',ls='--')
ax[1].set_xlabel(r'$\Delta$ (MHz)'); ax[1].set_ylabel(r'median arg$(c_2)$ (deg)'); ax[1].set_title(r'B: arg$(c_2)$ vs $\Delta$')
fr=np.nansum(Use&(np.abs(ARG2)>np.pi/2),axis=0)/np.maximum(np.sum(Use,axis=0),1)
ax[2].plot(D,fr,'o-',ms=3); ax[2].axvline(74.5,color='r',ls='--')
ax[2].set_xlabel(r'$\Delta$ (MHz)'); ax[2].set_ylabel(r'frac |arg$(c_2)|>\pi/2$'); ax[2].set_title(r'B: $c_2$ $\pi$-jump fraction vs $\Delta$')
plt.tight_layout(); plt.savefig(OUT+'fig_B_delta.png',dpi=110); plt.close()
print('figures saved')
def stat(Z,name,wrap=False):
    v=Z.copy().astype(float); v[~Use]=np.nan
    if wrap: v=np.angle(np.exp(1j*v))
    ix=np.unravel_index(np.nanargmax(v),v.shape); jx=np.unravel_index(np.nanargmin(v),v.shape)
    print('B %-6s max %+.4g at (a1=%.0f, f1=%d, D=%.0f); min %+.4g at (a1=%.0f, f1=%d, D=%.0f); median %+.4g'%(
        name,v[ix],a1s[ix[0]],f1s[ix[1]],4150-f1s[ix[1]],v[jx],a1s[jx[0]],f1s[jx[1]],4150-f1s[jx[1]],np.nanmedian(v)))
stat(A1,'A1'); stat(A2,'A2'); stat(R21,'R21'); stat(A0,'A0')
stat(np.load(OUT+'B_PHI0.npy'),'phi0',wrap=True); stat(ARG2,'arg2',wrap=True)
