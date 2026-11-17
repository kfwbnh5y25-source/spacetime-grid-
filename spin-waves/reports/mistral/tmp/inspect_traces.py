import h5py, numpy as np
path = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
with h5py.File(path,'r') as f:
    for k in ['324_275_4069_4112_1.310_1','324_275_4069_4112_4.550_1','0_275_0_4112_0.000_1','0_275_0_4094_0.000_1']:
        tr = f[k]['trace'][:]
        print(k, "shape", tr.shape)
        fr = tr[:,0,:]; pw = tr[:,1,:]
        print("  freq0/1000:", fr[0,0]/1e9, fr[0,-1]/1e9, "GHz; n unique across passes:", np.unique(fr).size)
        print("  power min/med/max:", pw.min(), np.median(pw), pw.max())
        # peak freq
        pm = pw.mean(axis=0)
        ipk = pm.argmax()
        print("  peak freq:", fr[0,ipk]/1e9, "GHz, peak pow:", pm[ipk])
        print("  ts:", f[k]['timestamp'][()].isoformat() if hasattr(f[k]['timestamp'][()],'isoformat') else f[k]['timestamp'][()])
