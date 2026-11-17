import sys
sys.path.insert(0, "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/pipeline")
import numpy as np, h5py
H5="/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
def loader(k, tr):
    return tr[:,0,:].mean(axis=0), tr[:,1,:].mean(axis=0)
with h5py.File(H5,'r') as f:
    for k in ["324_265_4069_4103_1.310_1","0_265_0_4103_0.000_1"]:
        tr = f[k]["trace"][:]
        fr,pw = loader(k,tr)
        d = np.asarray(tr)
        print(k, "dtype", tr.dtype, "fr shape", fr.shape, "pw shape", pw.shape,
              "pw max", pw.max() if pw.size else "EMPTY", "n>0", int((pw>0).sum()),
              "n nan", int(np.isnan(pw).sum()), "n inf", int(np.isinf(pw).sum()))
        print("   fr range GHz", fr.min()/1e9, fr.max()/1e9)
