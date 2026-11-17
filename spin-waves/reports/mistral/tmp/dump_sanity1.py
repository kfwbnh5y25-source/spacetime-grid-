import h5py
import numpy as np

path = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/measure/meas0_sanity.h5"
f = h5py.File(path, "r")
g = f["SAN_OVL_t0_p0.0000__SIG"]
print("SIG dataset members:", list(g.keys()))
for name in g:
    arr = np.array(g[name])
    try:
        print(f"  {name}: shape={arr.shape} dtype={arr.dtype} "
              f"min={arr.min():.6g} max={arr.max():.6g}")
    except TypeError:
        print(f"  {name}: shape={arr.shape} dtype={arr.dtype} value={arr[:60]}")

if "SAN_OVL_t0_p0.0000__BL" in f:
    b = f["SAN_OVL_t0_p0.0000__BL"]
    for name in b:
        arr = np.array(b[name])
        try:
            print(f"  BL {name}: shape={arr.shape} min={arr.min():.6g} max={arr.max():.6g}")
        except TypeError:
            print(f"  BL {name}: shape={arr.shape} dtype={arr.dtype} value={arr[:60]}")

s = np.array(g["trace"])[0] if "trace" in g else None
print("row0 head:", np.array(g["trace"][0, :5]) if "trace" in g else "-")
print("row1 head:", np.array(g["trace"][1, :5]) if "trace" in g else "-")

sg = np.array(g["trace"][:, 1], dtype=float)
sr = np.array(g["trace"][:, 0], dtype=float)
print("col0 like freq?", sr[:3], "... monotonic:", np.all(np.diff(sr) > 0))
print("SIG col1: mean", sg.mean(), "std", sg.std())

key_bl = "SAN_OVL_t0_p0.0000__BL"
if key_bl in f:
    bt = np.array(f[key_bl]["trace"][:, 1], dtype=float)
    print("BL col1: mean", bt.mean(), "std", bt.std())
    print("sig==bl exactly?", np.array_equal(sg, bt),
          " corr:", np.corrcoef(sg, bt)[0, 1])
    dmed = np.median(sg - bt)
    print("median(sig-bl):", dmed, "  max(sig-bl):", (sg - bt).max())
f.close()
