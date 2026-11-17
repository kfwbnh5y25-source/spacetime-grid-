import h5py, numpy as np
from scipy.signal import medfilt
with h5py.File('phase gating/3/3.h5','r') as f:
    for key in ['1000_590','400_590','200_590','1000_100','400_100']:
        tr=f[key]['trace'][0]
        s=tr[1].astype(np.float64)
        fm=medfilt(s,15)
        print(key,'mean=%.3g med=%.3g mad=%.3g max=%.3g min=%.3g  fm_max=%.3g fm_med=%.3g fm_argmax_f=%.2fMHz'%(
            s.mean(), np.median(s), 1.4826*np.median(np.abs(s-np.median(s))), s.max(), s.min(),
            fm.max(), np.median(fm), tr[0][np.argmax(fm)]/1e6))
        print('   raw sample [480:500]:', np.array2string(s[480:500]*1e12, precision=2))
