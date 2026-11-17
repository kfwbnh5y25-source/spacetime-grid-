import numpy as np, pickle
OUT='reports/mistral/tmp/'
amps=np.load(OUT+'A_amps.npy'); ls=np.load(OUT+'A_ls.npy')
P=np.load(OUT+'A_Pnet.npy'); Snr=np.load(OUT+'A_Snr.npy')
res=pickle.load(open(OUT+'A_ringup_fits.pkl','rb'))
fmap={d['a']:d for d in res}
lines=[]
lines.append('| a (mV) | P(l=0.938 ns) [W] | SNR | P(l=9.219 ns) [W] | SNR | $\\tau_r$ full [ns] | $\\tau_r$ from-onset [ns] | $\\chi^2_{\\rm red}$ |')
lines.append('|---|---|---|---|---|---|---|---|')
for ia,a in enumerate(amps):
    if a>560: break
    d=fmap.get(a,{})
    i0,i1=6,57  # l=60 -> 0.9375 ns ; l=590 -> 9.219
    tf='%.2f$\\pm$%.2f'%(d['tau'],d['tau_err']) if 'tau' in d else '--'
    to='%.2f'%d['tau_ons'] if 'tau_ons' in d else '--'
    chi='%.0f'%d['chi2'] if 'chi2' in d else '--'
    lines.append('| %d | %.3g | %.1f | %.3g | %.1f | %s | %s | %s |'%(a,P[ia,i0],Snr[ia,i0],P[ia,i1],Snr[ia,i1],tf,to,chi))
open(OUT+'tab_A_raw.md','w').write('\n'.join(lines))
print('\n'.join(lines[:8])); print('...rows:',len(lines)-2)
# threshold table compact
ath=np.load(OUT+'A_ath_persist.npy')
rows=[]
for i0 in range(0,58,2):
    seg=[]
    for i in range(i0,min(i0+2,58)):
        v=('%.0f'%ath[i]) if np.isfinite(ath[i]) else '--'
        seg.append('%d ns / l=%d: %s mV'%(ls[i]/64,ls[i],v))
    rows.append(' | '.join(seg))
open(OUT+'tab_A_ath.md','w').write('\n'.join(rows))
print('\n'.join(rows[:5]))
