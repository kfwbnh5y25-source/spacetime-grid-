import numpy as np, pickle
OUT='reports/mistral/tmp/'
amps=np.load(OUT+'A_amps.npy'); ls=np.load(OUT+'A_ls.npy')
P=np.load(OUT+'A_Pnet.npy'); Ptot=np.load(OUT+'A_Ptot.npy'); Fw=np.load(OUT+'A_fwhm.npy')
Bg=np.load(OUT+'A_Bg.npy'); Snr=np.load(OUT+'A_Snr.npy')
l_ns=ls/64.0
print('Ptot (1e-12 W·bin) per a at l=100(1.56ns) / 300(4.69) / 590(9.22):')
for a in [200,300,330,350,400,450,500,600,800,1000]:
    ia=np.where(amps==a)[0][0]
    print('  a=%4d: %.4g | %.4g | %.4g   (Bg median: %.2g/%.2g/%.2g)'%(a,Ptot[ia,8],Ptot[ia,28],Ptot[ia,-1],
          Bg[ia,8],Bg[ia,28],Bg[ia,-1]))
print()
print('FWHM(l) at fixed a (MHz):  a=400:')
print('  ', ' '.join('%.1f@%.1f'%(Fw[np.where(amps==400)[0][0],i],l_ns[i]) for i in range(0,58,6)))
ia=np.where(amps==600)[0][0]
print('a=600:', ' '.join('%.1f@%.1f'%(Fw[ia,i],l_ns[i]) for i in range(0,58,6)))
ia=np.where(amps==1000)[0][0]
print('a=1000:', ' '.join('%.1f@%.1f'%(Fw[ia,i],l_ns[i]) for i in range(0,58,6)))
# compare burst-limited 1/T widths
print('\nspectral dilution scale 1/T_burst (MHz):', ' '.join('%.0f@%.2f'%(64e9/(20*ls[i]*1.0)/1e6,l_ns[i]) for i in range(0,58,6)))
