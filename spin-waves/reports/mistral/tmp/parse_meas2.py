"""
parse_meas2.py - ANALYSIS-ONLY: threshold vs phase from meas2_threshold_vs_phase h5.

Keys: TVP_a2_{a2:.1f}_p{ph:.4f}_t{t} (SIG/BL) + REFT1_a2_0_p{..} (tone1-only refs).
Per phase: line power vs a2 -> logistic fit -> threshold a2,th(+) (inflection of the
rising edge, in mV). Then c12 fit of a2,th(phi) on the 16-phase uniform grid -> mean
threshold h0, modulation amplitude (W/V in mV units) and phi0. Writes table + figure.
Selftest injects known thresholds.
"""
import argparse
import os
import sys

import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def line_power(fr, sig, bl, center_hz, win_mhz=0.75):
    d = np.asarray(sig, float) - np.asarray(bl, float)
    win = np.abs(fr / 1e6 - center_hz / 1e6) <= win_mhz
    noise_bins = ~win & (d < np.quantile(d, 0.9))
    return float(d[win].sum()), float(np.std(d[noise_bins]))


def c12(points):
    ph = np.array([p[0] for p in points]); pw = np.array([p[1] for p in points])
    c1 = np.mean(pw * np.exp(-1j * ph))
    return float(2*np.abs(c1)), float(np.mean(pw)), float(np.angle(c1))


def threshold_of(a2s, ys, ref=None, nsig=4.0):
    """Rising edge: anchor on the steepest point (max slope) of the line-power
    curve, then refine with a logistic fit in a window around it. Returns
    (a2th_mV [logistic midpoint], quality [fit baseline])."""
    a2s = np.asarray(a2s, float); ys = np.asarray(ys, float)
    if len(a2s) < 5:
        return None
    ysm = np.convolve(ys, np.ones(3) / 3.0, mode="same")
    dy = np.gradient(ysm, a2s)
    c0 = int(np.argmax(dy))
    lo, hi = max(0, c0 - 8), min(len(a2s), c0 + 9)
    try:
        from scipy.optimize import curve_fit
        def sig_(x, L, x0, k, b):
            return b + (L - b) / (1 + np.exp(-k * (x - x0)))
        p, _ = curve_fit(sig_, a2s[lo:hi], ysm[lo:hi],
                         p0=[ysm.max(), a2s[c0], 0.5, ysm.min()], maxfev=5000)
        return float(p[1]), float(p[3])
    except Exception:
        return float(a2s[c0]), float(ysm.min())


def selftest():
    path = "/tmp/selftest_meas2.h5"
    A2g = np.arange(268.0, 314.0, 2.0)
    TH = {round(ph, 4): 285.0 + 8.0 * np.cos(round(ph, 4) + 1.0)
          for ph in np.arange(16) * 2*np.pi/16}
    center = 2055.5e9
    rng_shift = np.random.default_rng(5)
    with h5py.File(path, "w") as g:
        fr = center + np.linspace(-2.5e6, 2.5e6, 401)
        for ph, th in TH.items():
            for t in range(2):
                for a2 in A2g:
                    base = np.abs(rng_shift.normal(5e-8, 8e-9, len(fr)))
                    win = np.abs(fr - center) <= 0.4e6
                    sig = base.copy()
                    p_line = 30e-6 / (1 + np.exp(-(a2 - th) / 2.0))
                    sig[win] += p_line / win.sum()
                    key = f"TVP_a2_{a2:.1f}_p{ph:.4f}_t{t}"
                    g.create_dataset(key + "__SIG/trace", data=np.stack([fr, sig])[None])
                    g.create_dataset(key + "__BL/trace", data=np.stack([fr, base])[None])
        res = main_parse(path, center)
        rows = res["rows"]
        okm = np.mean([abs(r[1] - TH[r[0]]) < 3.0 for r in rows])
        A1, A0, phi0 = c12([(r[0], r[1]) for r in rows])
        want1, want0 = 8.0, 285.0
        ok = len(rows) == 16 and okm >= 0.8 and abs(A1 - want1) < 3 and abs(A0 - want0) < 4
        print(f"selftest: thresholds recovered {okm*16:.0f}/16, "
              f"A0={A0:.1f} (want {want0}), A1={A1:.1f} (want {want1}) "
              f"-> {'PASS' if ok else 'FAIL'}")
        sys.exit(0 if ok else 1)


def main_parse(path, center):
    data = {}
    with h5py.File(path, "r") as f:
        for k in f:
            if not k.endswith("__SIG"):
                continue
            name = k[:-5]
            fr = f[k + "/trace"][0, 0]
            sig = f[k + "/trace"][0, 1]
            bl = f[name + "__BL/trace"][0, 1]
            p, n = line_power(fr, sig, bl, center)
            if name.startswith("REFT1"):
                continue
            try:
                parts = name.split("_")
                # TVP_a2_285.0_p0.3927_t0 -> parts[1]='a2', parts[2]=value
                a2 = float(parts[2])
                ph = float(parts[3][1:])
            except (ValueError, IndexError):
                continue
            data.setdefault(ph, {}).setdefault(a2, []).append((p, n))
    rows = []
    for ph in sorted(data):
        a2s = np.array(sorted(data[ph]), float)
        pw = np.array([np.mean([v[0] for v in data[ph][a]]) for a in a2s])
        ns = np.array([np.mean([v[1] for v in data[ph][a]]) for a in a2s])
        th = threshold_of(a2s, pw, ref=float(np.min(ns)) if len(ns) else 0.0)
        if th is not None:
            rows.append((ph, th[0]))
    return {"rows": rows, "data": data}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--out", default="parse_meas2_out.md")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    md = ["# parse_meas2: a2,th(phi)\n",
          "| file | phi deg | a2th mV |", "|---|---|---|"]
    for path in args.files:
        with h5py.File(next(iter([path])), "r") as f0:
            fr = f0[list(f0)[0] + "/trace"][0, 0]
            name = list(f0)[0][:-5]
            center = float(np.median(fr))
        res = main_parse(path, center)
        rows = res["rows"]
        md.append(f"\n## {os.path.basename(path)}\n")
        for ph, th in rows:
            md.append(f"| {np.degrees(ph):.1f} | {th:.1f} |")
        if len(rows) >= 8:
            A1, A0, phi0 = c12(rows)
            md.append(f"\nA0 = {A0:.1f} mV, modulation = {A1:.1f} mV, "
                      f"phi0 = {np.degrees(phi0):+.1f} deg, A1/A0 = {A1/A0:.3f} "
                      f"(threshold ANISOTROPY; expect W/V pattern)")
    with open(args.out, "w") as g:
        g.write("\n".join(md))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
