import numpy as np
OUT='reports/mistral/tmp/'
amps=np.load(OUT+'A_amps.npy'); ls=np.load(OUT+'A_ls.npy'); P=np.load(OUT+'A_Pnet.npy'); Snr=np.load(OUT+'A_Snr.npy')
Fw=np.load(OUT+'A_fwhm.npy'); Bg=np.load(OUT+'A_Bg.npy')
for a in [400,420,430,450,500]:
    ia=np.where(amps==a)[0][0]
    print('a=%d: Pnet 0.938ns=%.4g (snr %.1f) | 3.44ns=%.4g (snr %.1f) | 9.22ns=%.4g (snr %.1f) | ratio9.2/3.4=%.2f | Bg=%.2g'%(
        a,P[ia,6],Snr[ia,6],P[ia,20],Snr[ia,20],P[ia,57],Snr[ia,57],P[ia,57]/P[ia,20],Bg[ia,57]))
