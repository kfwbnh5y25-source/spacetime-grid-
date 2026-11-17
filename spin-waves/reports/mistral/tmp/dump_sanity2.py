import h5py
import numpy as np

path = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/measure/meas0_sanity.h5"
f = h5py.File(path, "r")
for key in ["SAN_REF2_t0", "SAN_OVL_t0_p0.0000", "SAN_OVL_t1_p0.0000", "SAN_d-2.00_p0.0000"]:
    if key + "__BL" not in f:
        print(key, "-> no BL")
        continue
    s = np.mean(np.array(f[key + "__SIG"]["trace"])[:, 1, :], axis=0)
    b = np.mean(np.array(f[key + "__BL"]["trace"])[:, 1, :], axis=0)
    d = s - b
    print(f"{key}: mean SIG={s.mean():.3e} mean BL={b.mean():.3e} "
          f"bit-equal={np.array_equal(s, b)} corr={np.corrcoef(s, b)[0,1]:.4f} "
          f"|d|max={np.abs(d).max():.3e}")
f.close()
