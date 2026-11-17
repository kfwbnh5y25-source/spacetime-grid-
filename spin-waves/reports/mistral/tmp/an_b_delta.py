import numpy as np, pickle
OUT='reports/mistral/tmp/'
a1s=np.load(OUT+'B_a1s.npy'); f1s=np.load(OUT+'B_f1s.npy')
A1=np.load(OUT+'B_A1.npy'); A2=np.load(OUT+'B_A2.npy'); R21=np.load(OUT+'B_R21.npy')
A0=np.load(OUT+'B_A0.npy'); ARG2=np.load(OUT+'B_ARG2.npy'); PHI0=np.load(OUT+'B_PHI0.npy')
CORR=np.load(OUT+'B_CORRUPT.npy')
np.save=1
print('a1 grid:',a1s.astype(int)[::1], 'steps:',sorted(set(np.diff(a1s).astype(int)))[:5])
ok=~CORR&np.isfinite(A1)
D=4150.0-f1s  # per column
print('\nDelta coverage (ok blocks): %.0f .. %.0f MHz'%(D[ok.any(0)].min(),D[ok.any(0)].max()))
# Delta-binned means over a1
Ad=D
print('\nDelta-bin  n   A1_mean     A2/A1_med   frc|arg2|>pi/2  arg2_med(deg)')
bins=[(30,40),(40,50),(50,60),(60,70),(70,80),(80,90),(90,100),(100,110),(110,131)]
for lo,hi in bins:
    m=ok&((D>=lo)&(D<hi))[None,:]
    if not m.any(): print(f'{lo}-{hi}: none'); continue
    v1=A1[m]; v2=R21[m]; g=ARG2[m]
    print('%3d-%3d %4d  %.3g   %.3f     %.2f        %+.2f'%(lo,hi,m.sum(),np.nanmean(v1),np.nanmedian(v2),
          np.mean(np.abs(g[~np.isnan(g)])>np.pi/2),np.nanmedian(np.angle(np.exp(1j*g)))/np.pi*180))
# pi jump line near Delta=74.5: columns f1 = 4075 (D=75) and 4076 (D=74) + neighbors
print('\nf1 near 4075-4076 detail (Delta=74/75):')
for f1 in [4073,4074,4075,4076,4077,4078]:
    j=np.where(f1s==f1)[0]
    if j.size==0: print(' f1=%d: absent'%f1); continue
    j=j[0]; col=~CORR[:,j]&np.isfinite(A1[:,j])
    if not col.any(): print(' f1=%d: no ok blocks'%f1); continue
    g=ARG2[col,j]; r=R21[col,j]; v1=A1[col,j]
    print(' f1=%d (D=%.0f): n=%d arg2 mean=%+.2frad frac_pi/L=%.2f A1med=%.2g R21med=%.2f'%(
        f1,D[j],col.sum(),np.mean(np.angle(np.exp(1j*g)))/np.pi,np.mean(np.abs(g)>np.pi/2),np.median(v1),np.median(r)))
# A1 peak location vs Delta at highest a1:
for a in [328,320,300,280,256]:
    i=np.where(a1s==a)[0][0]
    row=np.where(~CORR[i]&np.isfinite(A1[i]))[0]
    if row.size==0: continue
    jbest=row[np.nanargmax(A1[i,row])]
    print('a1=%d: A1 peaks at f1=%d (Delta=%.0f MHz)'%(a,f1s[jbest],D[jbest]))
# A0 vs a1 at fixed f1=4075
j=np.where(f1s==4075)[0][0]
print('\nA0, A1, R21 vs a1 at f1=4075:')
for i,a in enumerate(a1s):
    if np.isfinite(A1[i,j]):
        print('  a1=%3d A0=%.3g A1=%.3g R21=%.2f arg2=%+.1f°%s'%(a,A0[i,j],A1[i,j],R21[i,j],ARG2[i,j]/np.pi*180,' [CORR]' if CORR[i,j] else ''))
