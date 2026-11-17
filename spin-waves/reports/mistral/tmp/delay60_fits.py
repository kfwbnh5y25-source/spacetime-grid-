"""Fits for the 60-phase delay recheck: tau windows, half-contrast time, phi0 slew.

Uses reports/mistral/tmp/delay60_{name}_harmonics.npz produced by delay60_extract.py.
Fits use clean+repaired blocks (excluded_step and low_n dropped).
"""
import numpy as np
from scipy.optimize import curve_fit

BASE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp"
PRED_SLEW = {"delay": 360 * 61 / 1000, "delay2": 360 * 74 / 1000}
TAU0S = [0.0, 1.5, 2.0, 3.0, 5.0, 6.0]  # 5, 6 ns = supplementary (decay-onset window, delay.h5)
BND_LO = [0.0, 1e-4, -50e-6]           # a >= 0, tau > 0, c >= -50 uW  (W units)
BND_HI = [1e-2, 1e7, 1e-2]


def expdec(t, a, tau, c):
    return a * np.exp(-t / tau) + c


def gaussdec(t, a, sg, c):
    return a * np.exp(-t**2 / (2 * sg**2)) + c


def algdec(t, a, ta, c):
    return a / (1 + (t / ta) ** 2) + c


def linear(t, m, b):
    return m * t + b


def wfit(model, t, y, p0, w=None, bounds=None):
    """Weighted fit; sigma = 1/w. Returns (popt, perr, rms, wchi2/dof) or (None, err, ...)."""
    try:
        kw = dict(maxfev=40000)
        if bounds is not None:
            kw["bounds"] = bounds
        popt, pcov = curve_fit(model, t, y, p0=p0, sigma=1 / w if w is not None else None,
                               absolute_sigma=True, **kw)
    except Exception as e:
        return None, str(e), None, None
    r = y - model(t, *popt)
    rms = float(np.sqrt(np.mean(r**2)))
    dof = max(len(t) - len(popt), 1)
    wchi2 = float(np.sum((r * w) ** 2)) if w is not None else np.nan
    perr = np.sqrt(np.diag(pcov))
    return (popt, perr, rms, wchi2 / dof), None, None, None


def main():
    lines_all = []
    for name in ("delay", "delay2"):
        d = np.load(f"{BASE}/delay60_{name}_harmonics.npz", allow_pickle=True)
        t, A0, A1, A2, phi0, A3, A4 = d["t"], d["A0"], d["A1"], d["A2"], d["phi0"], d["A3"], d["A4"]
        status = d["status"]
        sigma_S = float(d["sigma_S"])
        delta = float(d["delta_MHz"])
        ok = (status == "clean") | (status == "repaired")
        sigma_A1 = sigma_S * np.sqrt(2.0 / d["nph"])   # DFT bin error from per-sweep noise
        use = ok & np.isfinite(A1)
        print(f"\n########## {name} (Delta = {delta:.0f} MHz) ##########")
        print(f"blocks used: {use.sum()} / {len(t)}  (sigma_A1 ~ {sigma_A1[use].mean()*1e9:.2f} nW)")

        pos = use & (t >= 0)
        tp, y = t[pos], A1[pos]
        w = 1.0 / sigma_A1[pos]

        # --- tau fits in windows
        print("\n-- A1 = a*exp(-t/tau)+c on positive delays (bounds: a>=0, tau in [1e-4,1e7], c>=-50uW) --\r")
        fits = {}
        for t0 in TAU0S:
            m = use & (t >= max(t0, 0.0))
            tt, yy, ww = t[m], A1[m], 1.0 / sigma_A1[m]
            cbase = max(float(np.mean(A1[use & (t >= 12)])), 1e-9)
            p0 = [max(float(np.mean(A1[use & (t >= t0) & (t < t0 + 1)]) - cbase), 1e-9), 3.0, cbase]
            res, err, _, _ = wfit(expdec, tt, yy, p0, w=ww, bounds=(BND_LO, BND_HI))
            if res is None:
                print(f"{t0:5.1f} FIT FAILED: {err}")
                continue
            popt, perr, rms, wred = res
            fits[t0] = (popt, perr)
            tag = "" if t0 <= 3 else "  [supplementary: decay-onset window]"
            print(f"{t0:5.1f} {popt[0]*1e6:8.2f} {popt[1]:9.3f} {perr[1]:8.3f} {popt[2]*1e6:8.2f} {rms*1e6:8.2f} {wred:9.2f} {len(tt):4d}{tag}")

        # --- weighted-fit cross-checks: model comparison on t>=0
        print("\n-- fit-quality cross-checks (t>=0 window, A1-weighted, bounded) --")
        m = use & (t >= 0)
        tt, yy, ww = t[m], A1[m], 1.0 / sigma_A1[m]
        cbase = max(float(np.mean(A1[use & (t >= 12)])), 1e-9)
        p00 = [max(float(np.mean(A1[use & (t >= 0) & (t <= 1)]) - cbase), 1e-9), 3.0, cbase]
        for nm, mdl in [("exponential a*exp(-t/tau)+c", expdec),
                        ("gaussian   a*exp(-t^2/2sg^2)+c", gaussdec),
                        ("algebraic  a/(1+(t/ta)^2)+c", algdec)]:
            res, err, _, _ = wfit(mdl, tt, yy, p00, w=ww, bounds=(BND_LO, BND_HI))
            if res:
                popt, perr, rms, wred = res
                print(f"   {nm:32s}: tau/sig/ta={popt[1]:9.3f}+-{perr[1]:7.3f} ns  c={popt[2]*1e6:6.2f} uW  wchi2/dof={wred:8.2f}")
            else:
                print(f"   {nm:32s}: FAILED {err}")

        # --- plateau + half-contrast time
        plat = float(np.mean(A1[use & (t >= 1) & (t <= 2)]))
        half = plat / 2
        print(f"\nplateau A1(1-2 ns) = {plat*1e6:.2f} uW ; half = {half*1e6:.2f} uW")
        # data-based crossing: bin-average A1 in 1 ns bins after 1 ns
        edges = np.arange(1, 16.0, 1.0)
        centers, avgs = [], []
        for lo in edges[:-1]:
            m = use & (t >= lo) & (t < lo + 1)
            if m.sum() >= 5:
                centers.append(lo + 0.5)
                avgs.append(np.mean(A1[m]))
        centers = np.array(centers); avgs = np.array(avgs)
        thalf_data = np.nan
        for i in range(1, len(avgs)):
            if avgs[i - 1] >= half >= avgs[i]:
                # linear interpolation
                thalf_data = centers[i - 1] + (centers[i] - centers[i - 1]) * (avgs[i - 1] - half) / (avgs[i - 1] - avgs[i])
                break
        if not np.isfinite(thalf_data):
            if avgs[-1] > half:
                thalf_data = np.inf
        th_str = f"{thalf_data:.2f} ns" if np.isfinite(thalf_data) else \
            ("NOT reached in range (A1 stays above half)" if thalf_data == np.inf else f"{thalf_data}")
        print(f"   binned-data half-contrast time: {th_str}")
        if 0.0 in fits:
            (a0, tau0, c0), _ = fits[0.0]
            if a0 > 0 and (plat / 2) > c0:
                t_half_fit = tau0 * np.log(a0 / (plat / 2 - c0))
                print(f"   fit-curve (t0=0 params) half-contrast time: {t_half_fit:.2f} ns")

        # also report A1 means at key ranges
        for lo, hi in [(0, 1.5), (4.5, 7), (7, 10)]:
            m = use & (t >= lo) & (t < hi)
            print(f"   A1[{lo}-{hi} ns] = {np.mean(A1[m])*1e6:7.2f} uW")

        # --- phi0 slew rate
        plat_w = float(np.mean(A1[use & (t >= -2) & (t <= 2)]))
        thr = 0.2 * max(plat_w, plat)
        msl = use & (t >= -2.0) & np.isfinite(phi0) & (A1 > thr)
        ts, ps, a1s = t[msl], phi0[msl], A1[msl]
        o = np.argsort(ts)
        ts, ps, a1s = ts[o], ps[o], a1s[o]
        pu = np.unwrap(ps)
        # weight by A1 (phi0 noise grows as A1 -> 0)
        sig = sigma_S * np.sqrt(2.0 / 60.0) / a1s  # ~ sigma_phi0
        r, err, _, _ = wfit(linear, ts, pu, [2 * np.pi * delta / 1e3, pu[0]], w=1 / sig)
        print(f"\n-- phi0 slew (t_gap >= -2 ns, A1 > {thr*1e6:.1f} uW = 20% plateau, {len(ts)} points) --")
        if r:
            popt, perr, rms, wred = r
            slope_deg = np.degrees(popt[0])
            slope_err = np.degrees(perr[0])
            pred = PRED_SLEW[name]
            print(f"   slope = {slope_deg:6.2f} +- {slope_err:.2f} deg/ns ; prediction 360*Delta/1000 = {pred:.2f} deg/ns ; ratio = {slope_deg/pred*100:.1f} %")
            print(f"   intercept = {np.degrees(popt[1]):.1f} deg")
            Delta_eff = popt[0] / (2 * np.pi) * 1e3
            print(f"   implied Delta_eff = {Delta_eff:.1f} MHz ({Delta_eff/delta*100:.1f}% of nominal)")
        else:
            print(f"   fit failed: {err}")
        # unweighted too
        r2, _, _, _ = wfit(linear, ts, pu, [2 * np.pi * delta / 1e3, pu[0]])
        if r2:
            print(f"   unweighted slope = {np.degrees(r2[0][0]):.2f} deg/ns")
        # linearity sub-windows
        for tlo, thi in [(-2.0, 2.0), (-2.0, 4.0), (0.0, 3.0)]:
            mm = (ts >= tlo) & (ts <= thi)
            if mm.sum() > 10:
                rw, _, _, _ = wfit(linear, ts[mm], pu[mm], [2 * np.pi * delta / 1e3, pu[mm][0]], w=1 / sig[mm])
                if rw:
                    print(f"   sub-window {tlo:+.0f}..{thi:+.0f} ns: slope = {np.degrees(rw[0][0]):6.2f} +- {np.degrees(rw[1][0]):.2f} deg/ns (n={mm.sum()})")

        # --- two-exponential model (is the decay multi-component?)
        def twoexp(t, a1, t1, a2, t2, c):
            return a1 * np.exp(-t / t1) + a2 * np.exp(-t / t2) + c
        m2 = use & (t >= 0)
        tt2, yy2 = t[m2], A1[m2]
        ww2 = 1.0 / sigma_A1[m2]
        res2, err2, _, _ = wfit(twoexp, tt2, yy2,
                                [130e-6, 1.5, 20e-6, 8.0, 1e-6], w=ww2,
                                bounds=([0, 1e-3, 0, 1e-3, -50e-6], [1e-2, 30.0, 1e-2, 30.0, 1e-2]))
        if res2:
            popt2, perr2, rms2, wred2 = res2
            print(f"\n-- two-exponential A1 = a1*e^(-t/t1)+a2*e^(-t/t2)+c (t>=0) --")
            print(f"   a1={popt2[0]*1e6:7.2f}+-{perr2[0]*1e6:.2f} uW  t1={popt2[1]:6.3f}+-{perr2[1]:.3f} ns")
            print(f"   a2={popt2[2]*1e6:7.2f}+-{perr2[2]*1e6:.2f} uW  t2={popt2[3]:7.3f}+-{perr2[3]:.3f} ns")
            print(f"   c ={popt2[4]*1e6:7.2f}+-{perr2[4]*1e6:.2f} uW  rms={rms2*1e6:.2f} uW  wchi2/dof={wred2:.2f}")

        # print A1 and A2 full rows summary
        print("\n-- A1(t), A2(t) [uW] every 1 ns (mean of available blocks, clean+repaired) --")
        for lo in np.arange(-8, 15, 1.0):
            m = use & (t >= lo) & (t < lo + 1)
            if m.sum():
                print(f"   t={lo:+5.1f}..{lo+1:+5.1f}: A1={np.mean(A1[m])*1e6:7.2f}  A2={np.mean(A2[m])*1e6:7.2f}  (n={m.sum()})")


if __name__ == "__main__":
    main()
