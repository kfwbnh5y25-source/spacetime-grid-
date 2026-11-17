"""Analyze 'phase gating/3': single-tone (f=4072) pulse-length l x amplitude a sweep,
SA centered at f/2, RBW 300 kHz, avgs=10, x10 pulse train.

Observable: P(l, a) = baseline-subtracted power in the spectrum. Sub-threshold response
build-up with l measures the ring-up time (--> 2*Gamma_eff); threshold a_th(l) tests the
single-mode buildup model.
"""
import numpy as np
import h5py

BASE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/3/3.h5"

with h5py.File(BASE, "r") as f:
    keys = sorted(f.keys())
    # main grid: l in 20..590 step 10, a in 200..1000 step 10 (file also holds 90 stray l<20 probes)
    ls = list(range(20, 600, 10))
    amps = list(range(200, 1010, 10))
    rec = {}
    for k in keys:
        pr = k.split("_")
        a, l = int(pr[0]), int(pr[1])
        if a in amps and l in ls:
            rec[(l, a)] = k
    print(f"{len(rec)}/{len(keys)} keys on main grid (l 20-590 x10, a 200-1000 x10)")
    print(f"l in ns: {min(ls)/64:.3f}..{max(ls)/64:.3f}")
    g0 = f[rec[(min(ls), max(amps))]]
    freq0 = g0["trace"][0, 0, :].astype(np.float64)
    span_eff = (freq0[-1] - freq0[0]) / 1e6
    center = 0.5 * (freq0[0] + freq0[-1]) / 1e6
    print(f"SA center {center:.1f} MHz, span {span_eff:.1f} MHz, df = {span_eff/1000:.3f} MHz")
    # idler line: integrate +/- 3 MHz around center (line width unknown; RBW 300kHz)
    win = (np.abs(freq0 / 1e6 - center) <= 3.0)
    print(f"line window bins: {win.sum()}")

    def line_power(k):
        tr = f[rec[k]]["trace"][0, 1, :].astype(np.float64)
        return float(tr[win].sum())

    # matrix P(l, a)
    P = np.full((len(ls), len(amps)), np.nan)
    for i, l in enumerate(ls):
        for j, a in enumerate(amps):
            if (l, a) in rec:
                P[i, j] = line_power((l, a))

    np.save("p3_line_matrix.npy", P)
    lns = np.array(ls) / 64.0
    # print a coarse view: rows = selected a, cols = selected l
    jj = np.linspace(0, len(amps) - 1, 12).astype(int)
    ii = np.linspace(0, len(ls) - 1, 14).astype(int)
    print("\nP(f/2) [arb], rows=a (mV), cols=l (ns):")
    print("a\\l  " + " ".join(f"{lns[i]:6.2f}" for i in ii))
    for j in jj:
        print(f"{amps[j]:4d} " + " ".join(f"{P[i, j]*1e12:6.2f}" for i in ii))

    # ring-up fit for sub-threshold-saturating amplitudes: P(l) = P_inf (1 - exp(-l/tau_r)) with
    # oscillatory correction; robust: fit on the a column with clear saturation
    from scipy.optimize import curve_fit
    print("\nring-up fits (P = c0 + P_inf*(1-exp(-l/tau_r))*cos-shift-free):")
    for j in range(0, len(amps), 8):
        y = P[:, j] * 1e12
        if np.all(np.isnan(y)):
            continue
        fin = np.isfinite(y)
        try:
            popt, _ = curve_fit(
                lambda x, pinf, tau, c0: c0 + pinf * (1 - np.exp(-x / np.abs(tau))),
                lns[fin], y[fin], p0=[max(y[fin]), 2.0, min(y[fin])], maxfev=20000)
            pinf, tau, c0 = popt
            resid = y[fin] - (c0 + pinf * (1 - np.exp(-lns[fin] / abs(tau))))
            print(f"a={amps[j]:4d} mV: P_inf={pinf:8.2f}  tau_r={abs(tau):6.3f} ns  c0={c0:6.2f}  rms={np.std(resid):5.2f}")
        except Exception as e:
            print(f"a={amps[j]:4d} mV: fit fail {e}")
