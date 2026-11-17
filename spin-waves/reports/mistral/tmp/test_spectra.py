import h5py, numpy as np
path = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
with h5py.File(path,'r') as f:
    def load(k):
        tr = f[k]['trace'][:]
        return tr[:,0,:].mean(axis=0), tr[:,1,:].mean(axis=0)
    for tag,keys in [("204/4137",["324_204_4069_4137_1.310_1","324_204_4069_4137_4.550_1","0_204_0_4137_0.000_1"]),
                     ("275/4112",["324_275_4069_4112_1.310_1","324_275_4069_4112_4.550_1","0_275_0_4112_0.000_1"]),
                     ("272/4107",["324_272_4069_4107_1.310_1","324_272_4069_4107_4.550_1","0_272_0_4107_0.000_1"])]:
        for k in keys:
            try: fr,pw = load(k)
            except KeyError: print(tag,k,"MISSING"); continue
            q = np.percentile
            print(f"{tag} {k}")
            print(f"   raw-sum={pw.sum():.4g}  median={np.median(pw):.3g}  max={pw.max():.3g}  min={pw.min():.3g}  p25={q(pw,25):.3g} p75={q(pw,75):.3g}")
            # top 10 bins
            idx = np.argsort(pw)[::-1][:8]
            print("   top freqs MHz:", np.round(fr[idx]/1e6,3), "pows:", np.array2string(pw[idx], precision=2))
    # ref line shape: is ref line at f2/2 too? top freqs
