import h5py, numpy as np, re, collections
path = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
pat = re.compile(r'^(\d+)_(\d+)_(\d+)_(\d+)_([\d.]+)_(.+)$')
with h5py.File(path,'r') as f:
    keys = sorted(f.keys())
    # pass counts on a sample
    pc = collections.Counter()
    tmin, tmax = None, None
    for k in keys[::997]:
        n = f[k]['trace'].shape[0]; pc[n]+=1
        ts = f[k]['timestamp'][()].ravel()[-1].decode()
        if tmin is None or ts < tmin: tmin = ts
        if tmax is None or ts > tmax: tmax = ts
    print("pass-count sample dist:", dict(pc), " t range:", tmin, tmax)
    # idler metric test: a1=324 grid subset rows
    def load(k):
        tr = f[k]['trace'][:]
        fr = tr[:,0,:].mean(axis=0); pw = tr[:,1,:].mean(axis=0)
        return fr, pw
    a2s = list(range(200,276,2))
    for name, ph in [("on",1.310),("off",4.550)]:
        powf = np.zeros((len(a2s), 76))
        for i,a2 in enumerate(a2s):
            for j,f2 in enumerate(range(4075,4151)):
                fr,pw = load(f'324_{a2}_4069_{f2}_{ph:.3f}_1')
                df = np.diff(fr).mean()
                base = np.median(pw)
                powf[i,j] = np.sum(np.clip(pw-base,0,None))*df
        print(name, "integ pmed  min/med/max:", powf.min(), np.median(powf), powf.max(), "pos frac:", (powf>0).mean())
        if name=="on": P_on=powf
        else:
            with np.errstate(divide='ignore', invalid='ignore'):
                C = 10*np.log10(P_on/powf)
            print("contrast min/med/max:", np.nanmin(C), np.nanmedian(C), np.nanmax(C))
    # row 275 detail
    row_f2 = range(4075,4151)
    con275=[]
    for f2 in row_f2:
        try:
            fo,po = load(f'324_275_4069_{f2}_1.310_1'); ff,pf = load(f'324_275_4069_{f2}_4.550_1')
            dfo = np.diff(fo).mean()
            Io = np.sum(np.clip(po-np.median(po),0,None))*dfo; If = np.sum(np.clip(pf-np.median(pf),0,None))*dfo
            con275.append((f2, Io, If, 10*np.log10(If/Io) if Io>0 and If>0 else np.nan))
        except Exception as e:
            print("missing", f2, e); con275.append((f2,0,0,np.nan))
    for c in con275: 
        if not np.isnan(c[3]) and abs(c[3])>3: print("big contrast f2=%d on=%.3g off=%.3g C=%.1f dB" % c)
