import h5py, numpy as np, collections, re, datetime, json
P = 'spin-waves/Figure!/arnold_tongue_zoom/data2.h5'
pat = re.compile(r'(\d+)_(\d+)_([a-z_]+?)_(\d{8}_\d{6})_(\w+?)_0$')
recs = []
with h5py.File(P, 'r') as f:
    for k in f.keys():
        m = pat.match(k)
        if not m: continue
        a, fr, lab, run, idx = int(m[1]), int(m[2]), m[3].replace('_', ''), m[4], m[5]
        g = f[k]
        tr = g['trace'][()]; ts = g['timestamp'][()]
        t = np.array([datetime.datetime.fromisoformat(x[0].decode()).timestamp() for x in ts])
        spec = tr[:, 1, :].astype(np.float64); fq = tr[0, 0, :]
        recs.append(dict(k=k, a=a, f=fr, lab=lab, run=run, idx=idx, t=t, spec=spec, fq=fq))
print('records', len(recs))
df = recs[0]['fq'][1] - recs[0]['fq'][0]
out = {}
for lab in ['halfaverage', 'full']:
    R = [r for r in recs if r['lab'] == lab]
    # integrated linear power per acquisition
    for r in R: r['P'] = r['spec'].sum(axis=1) * df
    # per-bin noise floor from a=0 nulls
    nulls = [r for r in R if r['a'] == 0]
    s = {}
    s['n_records'] = len(R)
    s['n_null'] = len(nulls)
    if nulls:
        allnull = np.concatenate([r['spec'] for r in nulls])
        s['null_bin_mean'] = float(allnull.mean()); s['null_bin_std'] = float(allnull.std())
        Pn = np.concatenate([r['P'] for r in nulls])
        s['null_P_mean'] = float(Pn.mean()); s['null_P_std'] = float(Pn.std())
    # within-key repeatability
    cv, mu = [], []
    for r in R:
        if len(r['P']) >= 2 and r['a'] > 0:
            mu.append(r['P'].mean()); cv.append(r['P'].std(ddof=1) / abs(r['P'].mean()) if r['P'].mean() != 0 else np.nan)
    mu = np.array(mu); cv = np.array(cv)
    ok = np.isfinite(cv)
    qs = np.quantile(mu[ok], [0, .25, .5, .75, .9, 1])
    s['within_cv_by_power_quantile'] = []
    for lo, hi in zip(qs[:-1], qs[1:]):
        sel = ok & (mu >= lo) & (mu <= hi)
        s['within_cv_by_power_quantile'].append(dict(lo=float(lo), hi=float(hi), n=int(sel.sum()), median_cv=float(np.median(cv[sel]))))
    # replicate groups across runs
    grp = collections.defaultdict(list)
    for r in R:
        if r['a'] > 0: grp[(r['a'], r['f'])].append(r)
    ratios, rel_between, rel_within, dts = [], [], [], []
    for key, rs in grp.items():
        if len(rs) < 2: continue
        means = np.array([r['P'].mean() for r in rs]); n_in = np.array([len(r['P']) for r in rs])
        within = np.mean([r['P'].var(ddof=1) / len(r['P']) for r in rs if len(r['P']) >= 2]) if any(len(r['P']) >= 2 for r in rs) else np.nan
        between = means.var(ddof=1)
        if within > 0 and np.mean(means) > 0:
            ratios.append(between / within)
            rel_between.append(np.sqrt(between) / np.mean(means))
            rel_within.append(np.sqrt(within) / np.mean(means))
            ts = np.array([r['t'].mean() for r in rs]); dts.append(ts.max() - ts.min())
    ratios = np.array(ratios)
    s['replicate_groups'] = len(ratios)
    s['between_over_within_var_ratio_median'] = float(np.median(ratios)) if len(ratios) else None
    s['between_over_within_var_ratio_mean_sqrt'] = float(np.sqrt(np.mean(ratios))) if len(ratios) else None
    s['rel_between_median'] = float(np.median(rel_between)) if len(rel_between) else None
    s['rel_within_sem_median'] = float(np.median(rel_within)) if len(rel_within) else None
    s['replicate_dt_s_median'] = float(np.median(dts)) if len(dts) else None
    # timing: seconds per repeat
    gaps = np.concatenate([np.diff(r['t']) for r in R if len(r['t']) >= 2])
    s['repeat_gap_s_median'] = float(np.median(gaps)); s['repeat_gap_s_p95'] = float(np.quantile(gaps, .95))
    # PCA effective rank of mean spectra (noise whitened by null std)
    M = np.array([r['spec'].mean(axis=0) for r in R if r['a'] > 0])
    sub = M[np.random.default_rng(0).choice(len(M), min(4000, len(M)), replace=False)]
    sc = s.get('null_bin_std', 1.0)
    sub = sub - sub.mean(axis=0)
    sv = np.linalg.svd(sub / sc, compute_uv=False)
    ev = sv ** 2; ev /= ev.sum()
    s['pca_cum_var_k'] = {str(k): float(ev[:k].sum()) for k in [1, 2, 3, 5, 10, 20]}
    s['participation_ratio'] = float(1 / (ev ** 2).sum())
    out[lab] = s
print(json.dumps(out, indent=1))
json.dump(out, open('/Users/davidzehner/.gemini/antigravity/brain/21fa05fd-fdbd-40bb-b413-bee9c27eaa97/scratch/repeatability.json', 'w'), indent=1)
