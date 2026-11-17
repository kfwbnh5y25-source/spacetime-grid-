import numpy as np, h5py
H5="/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
P0=np.load("../figures/fig3_overlap__P_on.npy"); P4=np.load("../figures/fig3_overlap__P_off.npy")
C=10*np.log10(P4/P0); A2=list(range(200,276)); F2=list(range(4075,4151))
line=P0>5e-7
rs=np.where(line.sum(axis=1)>=10, np.nanmean(np.abs(C),axis=1), -np.inf)
ri=int(np.argmax(rs)); ROW=A2[ri]
CR=C[ri].copy(); jm=int(np.nanargmax(CR)); ROW_F2=F2[jm]
print("ROW",ROW,"ROW_F2",ROW_F2,"C %.2f"%CR[jm], "line cells in row", int(line[ri].sum()))
with h5py.File(H5,'r') as f:
    for k in [f"324_{ROW}_4069_{ROW_F2}_1.310_1", f"324_{ROW}_4069_{ROW_F2}_4.550_1", f"0_{ROW}_0_{ROW_F2}_0.000_1"]:
        print(k, k in f.keys())
        if k in f.keys():
            tr=f[k]['trace'][:]
            pw=tr[:,1,:].mean(axis=0)
            print("   passes", tr.shape[0], "max %.3g"%pw.max(), "n>0", int((pw>0).sum()))
