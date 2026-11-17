"""June delivered-tone levels from delay2.h5 0EX/1EX sweeps vs today's probe.

Today (excitation path): f1(4075) peak -59.1 dBm, f2(4149) -43.3 dBm
-> ratio 15.8 dB, while (588/327)^2 = 3.23 -> expected ratio 5.1 dB.
So f1 delivery is ~10.7 dB below linear expectation TODAY.
Question: what was the 0EX(4075)/1EX(4149) ratio in June?
"""
import h5py
import numpy as np

BASE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/fig1_extract/"
F1, F2 = 4075e6, 4149e6
A1, A2 = 327.0, 588.0

for fn in ("delay2.h5", "delay.h5"):
    print("=" * 64)
    print(fn)
    peaks = {4075: [], 4149: []}
    n_key = {4075: 0, 4149: 0}
    with h5py.File(BASE + fn, "r") as f:
        for key in f.keys():
            excl = None
            for probe, center in (("0EX", F1), ("1EX", F2)):
                if key.endswith(probe):
                    excl = probe
                    tr = np.array(f[key]["trace"])
                    fr = np.mean(tr[:, 0, :], axis=0)
                    pw = np.mean(tr[:, 1, :], axis=0)
                    win = np.abs(fr - center) <= 2e6
                    if win.any():
                        p = float(pw[win].max())
                        n_key[center / 1e6] += 1
                        peaks[center / 1e6].append(p)
                    break
            if excl is None:
                # regular sweep (center f2/2) - not needed here
                continue
    print(f"  0EX(4075): n={len(peaks[4075])}, median peak = "
          f"{np.median(peaks[4075])*1e3:.3e} mW = {10*np.log10(max(np.median(peaks[4075]),1e-30)):.1f} dBm")
    print(f"  1EX(4149): n={len(peaks[4149])}, median peak = "
          f"{np.median(peaks[4149])*1e3:.3e} mW = {10*np.log10(max(np.median(peaks[4149]),1e-30)):.1f} dBm")
    if peaks[4075] and peaks[4149]:
        r = 10 * np.log10(np.median(peaks[4149]) / max(np.median(peaks[4075]), 1e-30))
        print(f"  June ratio 1EX/0EX = {r:.1f} dB  (linear expectation 5.1 dB)")
print()
print("today's probe: 1EX -43.3 dBm, 0EX -59.1 dBm, ratio 15.8 dB")
