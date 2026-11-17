import h5py, numpy as np, time
path = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
A2 = list(range(200,276)); F2 = list(range(4075,4151))
t0=time.time()
with h5py.File(path,'r') as f:
    grids = {}
    for ph in (1.310, 4.550):
        Ppk = np.full((76,76), np.nan); Pint = np.full((76,76), np.nan)
        for i,a2 in enumerate(A2):
            for j,f2 in enumerate(F2):
                tr = f[f'324_{a2}_4069_{f2}_{ph:.3f}_1']['trace'][:]
                pw = tr[:,1,:].mean(axis=0); fr = tr[:,0,:].mean(axis=0)
                df = np.diff(fr).mean()
                Ppk[i,j] = pw.max() - np.median(pw)
                Pint[i,j] = np.sum(np.clip(pw-np.median(pw),0,None))*df
        grids[ph]=(Ppk,Pint)
    P0pk,P0int = grids[1.310]; P4pk,P4int = grids[4.550]
    np.save("P_on_peak.npy",P0pk); np.save("P_off_peak.npy",P4pk)
    np.save("P_on_int.npy",P0int); np.save("P_off_int.npy",P4int)
print("load time", time.time()-t0)
for nm,(P0,P4) in [("peak",(P0pk,P4pk)),("int",(P0int,P4int))]:
    C = 10*np.log10(P4/P0)
    im,jm = np.unravel_index(np.nanargmax(C), C.shape)
    i2,j2 = np.unravel_index(np.nanargmin(C), C.shape)
    # row means of |C|
    rabs = np.nanmean(np.abs(C),axis=1)
    best_rows = np.argsort(rabs)[::-1][:5]
    print(nm, "C: min %.2f dB %s | max %.2f dB at a2=%d f2=%d | median %.2f" % (np.nanmin(C),(A2[i2],F2[j2]),C[im,jm],A2[im],F2[jm],np.nanmedian(C)))
    print("  top|C| rows:", [(A2[r], round(float(rabs[r]),2)) for r in best_rows])
