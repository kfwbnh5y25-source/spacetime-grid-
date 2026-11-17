import numpy as np
OUT='reports/mistral/tmp/'
amps=np.load(OUT+'A_amps.npy'); ls=np.load(OUT+'A_ls.npy')
P=np.load(OUT+'A_Pwin.npy'); Snr=np.load(OUT+'A_Snrint.npy')
Pk=np.load(OUT+'A_Pkpos.npy'); Fw=np.load(OUT+'A_fwhm.npy'); Bg=np.load(OUT+'A_Bg.npy')
for a in [200,250,300,320,330,350,400,450,500,600,800,1000]:
    ia=np.where(amps==a)[0][0]
    row=P[ia]/np.array([max(wsio,1)*1 for wsio in range(1)])  # dummy
    print(f'a={a:4d}: ', ' '.join(f'{P[ia,il]/1e-12:7.1f}' for il in range(0,ls.size,4)))
print()
print('cols l =', ' '.join(f'{ls[i]}({ls[i]/64:.2f}ns)' for i in range(0,ls.size,4)))
print()
for a in [200,250,300,330,350,400,1000]:
    ia=np.where(amps==a)[0][0]
    print(f'a={a:4d} Snr: ', ' '.join(f'{Snr[ia,il]:6.1f}' for il in range(0,ls.size,4)))
print()
for a in [300,330,350,400,600,1000]:
    ia=np.where(amps==a)[0][0]
    print(f'a={a:4d} pkpos: ', ' '.join(f'{Pk[ia,il]:8.2f}' for il in range(0,ls.size,4)))
    print(f'         fwhm : ', ' '.join(f'{Fw[ia,il]:8.3f}' for il in range(0,ls.size,4)))
