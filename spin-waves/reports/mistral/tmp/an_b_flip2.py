import numpy as np
OUT='reports/mistral/tmp/'
a1s=np.load(OUT+'B_a1s.npy'); f1s=np.load(OUT+'B_f1s.npy')
A1=np.load(OUT+'B_A1.npy'); A2=np.load(OUT+'B_A2.npy'); ARG2=np.load(OUT+'B_ARG2.npy'); CORR=np.load(OUT+'B_CORRUPT.npy')
Use=np.isfinite(A1)&~CORR
D=4150-f1s
print('locus of adjacent-D sign flips of cos(arg(c2)) weighted by A2:')
for i,a in enumerate(a1s):
    col=np.where(Use[i])[0]
    if col.size<10: continue
    col=np.sort(col)
    s=np.cos(ARG2[i,col])*A2[i,col]
    flips=[0.5*(D[col[k]]+D[col[k+1]]) for k in range(len(col)-1) if s[k]*s[k+1]<0]
    if flips: print('  a1=%3d: D~'%a, ' '.join('%.0f'%f for f in flips))
