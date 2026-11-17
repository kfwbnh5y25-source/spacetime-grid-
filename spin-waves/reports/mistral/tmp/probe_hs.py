import h5py, numpy as np

for p in [
    "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/detune_fine031126.h5",
    "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5",
]:
    print("=" * 20, p.split("/")[-1])
    try:
        with h5py.File(p, 'r') as f:
            keys = sorted(f.keys())
            print("n_keys", len(keys), "first", keys[:5], "last", keys[-3:])
            k = keys[0]
            print("key attrs:", dict(f[k].attrs))
            tr = f[k]["trace"]
            print("trace shape", tr.shape, tr.dtype)
            print("file attrs:", dict(f.attrs))
            a0 = f[k]["trace"][:, 0, :]
            print("ch0 (freq?) sample:", a0[0, :3], "...", a0[-1, -3:])
            a1 = f[k]["trace"][:, 1, :]
            print("ch1 stats mean", a1.mean(), "min", a1.min(), "max", a1.max())
            if tr.shape[1] >= 3:
                a2 = f[k]["trace"][:, 2, :]
                print("ch2 stats mean", a2.mean(), "min", a2.min(), "max", a2.max())
    except Exception as e:
        print("ERROR:", repr(e))
