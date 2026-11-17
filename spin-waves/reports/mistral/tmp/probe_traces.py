import h5py, numpy as np

p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
f = h5py.File(p, 'r')
keys = list(f.keys())
rng = np.random.default_rng(0)
for k in [keys[0], keys[len(keys)//4], keys[len(keys)//2], keys[-1], '324_275_4069_4130_1.310_1', '324_275_4069_4130_4.550_1', '0_275_0_4130_0.000_1']:
    g = f[k]
    tr = g['trace'][:]
    ts = g['timestamp'][:]
    print(k, "shape", tr.shape, "ts", ts.ravel())
    for i in range(tr.shape[0]):
        for j in range(tr.shape[1]):
            seg = tr[i, j]
            print(f"  [{i},{j}] min {seg.min():.4g} max {seg.max():.4g} mean {seg.mean():.4g} std {seg.std():.4g}", end='')
        print()
f.close()
