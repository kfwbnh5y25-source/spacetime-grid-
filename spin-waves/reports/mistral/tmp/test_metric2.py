import h5py, numpy as np
path = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
with h5py.File(path,'r') as f:
    def load(k):
        tr = f[k]['trace'][:]
        return tr[:,0,:].mean(axis=0), tr[:,1,:].mean(axis=0)
    a2s = list(range(200,276,2))
    res = {}
    for name, ph in [("on",1.310),("off",4.550)]:
        Ppk = np.zeros((len(a2s),76)); Pint = np.zeros_like(Ppk); Praw = np.zeros_like(Ppk)
        for i,a2 in enumerate(a2s):
            for j,f2 in enumerate(range(4075,4151)):
                fr,pw = load(f'324_{a2}_4069_{f2}_{ph:.3f}_1')
                Ppk[i,j] = pw.max()
                Pint[i,j] = pw.sum()*np.diff(fr).mean()
                Praw[i,j] = np.sum(np.clip(pw,-0,None))  # raw sum no baseline
                if a2==200 and f2==4112:
                    print(k:="  debug", ph, "peak argmax idx:", pw.argmax(), "window center idx 500")
        res[name]=(Ppk,Pint)
    for mi,mname in enumerate(["peak","rawint"]):
        Po,Pf = res["on"][mi], res["off"][mi]
        with np.errstate(divide='ignore',invalid='ignore'):
            C = 10*np.log10(Pf/Po)
        im,jm = np.unravel_index(np.nanargmax(np.abs(C)), C.shape)
        print(mname, "C min/med/max: %.2f %.2f %.2f" % (np.nanmin(C),np.nanmedian(C),np.nanmax(C)),
              " max|C| %.2f dB at a2=%d f2=%d" % (C[im,jm], a2s[im], 4075+jm))
    # where is line center vs window center: on-grid sample
    offp=[]
    for i,a2 in enumerate([200,237,275]):
        for j,f2 in enumerate([4075,4112,4150]):
            fr,pw = load(f'324_{a2}_4069_{f2}_1.310_1')
            cen = fr[len(fr)//2]/1e6; pk = fr[pw.argmax()]/1e6
            offp.append(pk-cen)
            print(f"a2={a2} f2={f2} win-cent={cen:.3f}MHz peak={pk:.3f}MHz off={pk-cen:+.3f}MHz pk={pw.max():.3g}")
