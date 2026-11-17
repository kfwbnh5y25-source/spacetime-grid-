import numpy as np

base = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/"
for name in ["delay60_delay2_A1_t.npy", "delay60_delay2_A2_t.npy",
             "delay60_delay_A1_t.npy", "delay60_delay_A2_t.npy"]:
    try:
        a = np.load(base + name)
        print(name, a.shape)
        t, A, m = a[0], a[1], a[2].astype(bool)
        tt, AA = t[m], A[m]
        print("  t range:", tt.min(), tt.max(), " n=", len(tt))
        for probe in [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 6.0, 8.0]:
            i = np.argmin(np.abs(tt - probe))
            if abs(tt[i] - probe) <= 0.5:
                print(f"  A({tt[i]:+.1f}ns) = {AA[i]*1e6:8.2f} uW")
        print("  A max:", AA.max(), " mean:", AA.mean())
    except Exception as e:
        print(name, "ERR", e)
