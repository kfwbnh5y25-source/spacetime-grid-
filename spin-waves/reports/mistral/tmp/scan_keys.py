import h5py, numpy as np, re

p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/detune_fine031126.h5"
with h5py.File(p, 'r') as f:
    keys = sorted(f.keys())
    prefixes = sorted(set(re.match(r'(\S+)_(\S+)_', k).group(0) for k in keys))
    print("f1_f2 combos:", len(prefixes))
    print(prefixes[:40])
    rng = np.random.default_rng(0)
    sample = list(rng.choice(keys, 24, replace=False))
    for k in sample:
        tr = f[k]['trace'][:]
        fr, pw = tr[0, 0, :].astype(np.float64), tr[0, 1, :].astype(np.float64)
        fc = fr.mean()
        c = int(np.argmin(np.abs(fr - fc)))
        band = pw[c - 25:c + 25]
        noise = np.median(pw)
        print(f"{k}: fc={fc/1e9:.5f}GHz peak_center={band.max():.3e} "
              f"median={noise:.3e} ratio={band.max()/max(noise,1e-30):.1f}")
