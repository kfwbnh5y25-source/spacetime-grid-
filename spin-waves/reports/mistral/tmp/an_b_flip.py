import numpy as np
OUT='reports/mistral/tmp/'
a1s=np.load(OUT+'B_a1s.npy'); f1s=np.load(OUT+'B_f1s.npy')
A1=np.load(OUT+'B_A1.npy'); ARG2=np.load(OUT+'B_ARG2.npy'); CORR=np.load(OUT+'B_CORRUPT.npy')
Use=np.isfinite(A1)&~CORR
D=4150-f1s
print('per-a1: arg(c2) wrapped (deg) at selected D; and jump-locus search (|sign flip of cos(arg2)| between D bins):')
print('a1 :', ' '.join('D=%3d'%d for d in [34,50,58,66,74,75,82,90,100,108]))
for i,a in enumerate(a1s):
    row=[]
    for dt in [34,50,58,66,74,75,82,90,100,108]:
        j=np.argmin(np.abs(D-dt))
        if Use[i,j]: row.append('%4.0f'%(np.angle(np.exp(1j*ARG2[i,j]))/np.pi*180))
        else: row.append('  - ')
    print('%3d:'%a, ' '.join(row))
# jump locus: find adjacent-D pairs where arg2 near 0 vs near +-pi switches, weight by A2 significance
print('\nlocus (a1, D_edge): adjacent-D sign switches of cos(arg2)*A2:')
for i,a in enumerate(a1s):
    col=np.where(Use[i])[0]
    if col.size<10: continue
    col=np.sort(col)
    s=np.cos(ARG2[i,col])*A2[i,col]
    flips=[]
    for k in range(len(col)-1):
        if s[k]*s[k+1]<0:
            flips.append(0.5*(D[col[k]]+D[col[k+1]]))
    if flips: print('  a1=%3d flips at D~'%a, ' '.join('%.0f'%f for f in flips))
