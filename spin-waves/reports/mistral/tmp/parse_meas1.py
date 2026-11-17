"""
parse_meas1.py - ANALYSIS-ONLY parser for meas1_longgap_memory h5 outputs.

Usage:
  python parse_meas1.py <meas1_*.h5> [<more.h5> ...] [--out table.md] [--figdir figdir]
  python parse_meas1.py --selftest

Auto-detects the line center (median of per-trace argmax of signal-minus-baseline),
then per (config, delay): cos-fit of the 16-phase line powers (c12 convention) ->
A0, A1, A1/A0, phi0; exponential tau_fit A1 = m*exp(-max(t,0)/tau)+c (grid scan),
slew_fit dphi0/dt_gap (unwrapped linear fit). Writes an .md table + figure grid.
Robustness: REF2 keys excluded; noise median per delay reported so flat/dead blocks
are visible rather than silently dropped.
"""
import argparse
import os
import sys

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("h5py required")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def line_power(fr, sig, bl, center_hz, win_mhz=1.5):
    d = np.asarray(sig, float) - np.asarray(bl, float)
    win = np.abs(fr / 1e6 - center_hz / 1e6) <= win_mhz
    noise_bins = ~win & (d < np.quantile(d, 0.9))
    return float(d[win].sum()), float(np.std(d[noise_bins]))


def c12(points):
    ph = np.array([p[0] for p in points])
    pw = np.array([p[1] for p in points])
    c1 = np.mean(pw * np.exp(-1j * ph))
    return float(2 * np.abs(c1)), float(np.mean(pw)), float(np.angle(c1))


def detect_center(fr, sig, bl):
    d = np.asarray(sig, float) - np.asarray(bl, float)
    k = int(np.argmax(np.convolve(np.abs(d), np.ones(9), "same")))
    fr0 = fr[k]
    sel = np.abs(fr - fr0) <= 1.5e6
    w = np.abs(d[sel])
    return float(np.sum(fr[sel] * w) / np.sum(w))


def parse_file(path):
    out = {"delays": {}}
    with h5py.File(path, "r") as f:
        keys = [k for k in f if k.endswith("__SIG") and k.startswith("LGM_")]
        if not keys:
            return None
        centers = []
        for k in keys[:50]:
            name = k[:-5]
            centers.append(detect_center(f[k + "/trace"][0, 0],
                                         f[k + "/trace"][0, 1],
                                         f[name + "__BL/trace"][0, 1]))
        out["center_mhz"] = float(np.median(centers)) / 1e6
        for k in keys:
            name = k[:-5]
            fr = f[k + "/trace"][0, 0]
            sig = f[k + "/trace"][0, 1]
            bl = f[name + "__BL/trace"][0, 1]
            p, n = line_power(fr, sig, bl, out["center_mhz"] * 1e6)
            try:
                parts = name.split("_")
                ci = int(parts[1][1:])
                dly = float(parts[2][1:])
                ph = float(parts[3][1:])
            except (ValueError, IndexError):
                print(f"  [skip unparsed] {name}")
                continue
            out["delays"].setdefault(ci, {}).setdefault(dly, []).append((ph, p, n))
    return out


def fit_config(cfg):
    rows = []
    for d in sorted(cfg):
        pts = [(ph, p) for ph, p, _ in cfg[d]]
        if len(pts) < 8:
            continue
        A1, A0, phi0 = c12(sorted(pts))
        nz = float(np.median([n for _, _, n in cfg[d]]))
        rows.append((d, A0, A1, A1 / max(A0, 1e-30), phi0, nz))
    return np.array(rows)


def tau_fit(rows):
    t = np.array([r[0] for r in rows])
    y = np.array([r[2] for r in rows])
    best = None
    for tau in np.linspace(0.3, 150, 400):
        X = np.c_[np.exp(-np.maximum(t, 0) / tau), np.ones_like(t)]
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = float(np.sqrt(np.mean((X @ c - y) ** 2)))
        if best is None or r < best[0]:
            best = (r, float(tau), float(c[0]), float(c[1]))
    return best


def slew_fit(rows):
    """Alias-safe winding fit: maximize vector strength R(m)=|<e^{i(phi-m t)}>| over a
    slope grid, then report residual circular std. Robust to sparse/irregular delays."""
    t = np.array([r[0] for r in rows])
    ph = np.array([r[4] for r in rows])
    if len(t) < 5:
        return None
    best = None
    for m in np.linspace(-2.0, 2.0, 8001):
        R = np.abs(np.mean(np.exp(1j * (ph - m * t))))
        if best is None or R > best[0]:
            best = (R, m)
    R, m = best
    res = np.angle(np.exp(1j * (ph - m * t)))   # circular residual, (-pi, pi]
    circ_std = float(np.sqrt(np.mean(res ** 2)))
    return float(m), circ_std, float(R)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*")
    ap.add_argument("--out", default="parse_meas1_out.md")
    ap.add_argument("--figdir", default="parse_meas1_figs")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    if not args.files:
        ap.error("give at least one h5 file (or --selftest)")
    os.makedirs(args.figdir, exist_ok=True)
    md = ["# parse_meas1 tables\n"]
    for path in args.files:
        res = parse_file(path)
        if res is None:
            md.append(f"\n## {os.path.basename(path)}: NO LGM_ keys\n")
            continue
        md.append(f"\n## {os.path.basename(path)} (center {res['center_mhz']:.2f} MHz)\n")
        md.append("| cfg | t_gap ns | A0 uW | A1 uW | A1/A0 | phi0 deg | noise uW |")
        md.append("|---|---|---|---|---|---|---|")
        fig, axs = plt.subplots(1, 2, figsize=(11, 4))
        for ci, cfg in sorted(res["delays"].items()):
            rows = fit_config(cfg)
            for r in rows:
                md.append(f"| {ci} | {r[0]:+.2f} | {r[1]*1e6:.1f} | {r[2]*1e6:.1f} |"
                          f" {r[3]:.3f} | {np.degrees(r[4]):+.1f} | {r[5]*1e6:.2f} |")
            axs[0].plot(rows[:, 0], rows[:, 3], "o-", label=f"c{ci}")
            axs[1].plot(rows[:, 0], np.degrees(rows[:, 4]), "o-", label=f"cfg {ci}")
            tf, sl = tau_fit(rows), slew_fit(rows)
            md.append(f"\ncfg {ci}: tau_fit = {tf[1]:.2f} ns (rms {tf[0]*1e6:.1f} uW)")
            if sl:
                md.append(f"slew = {np.degrees(sl[0]):+.2f} deg/ns "
                          f"(residual {sl[1]:.2f} deg)")
        axs[0].set_xlabel("t_gap (ns)"); axs[0].set_ylabel("A1/A0"); axs[0].legend()
        axs[1].set_xlabel("t_gap (ns)"); axs[1].set_ylabel("phi0 (deg)"); axs[1].legend()
        fig.tight_layout()
        fig.savefig(os.path.join(args.figdir,
                    os.path.basename(path).replace(".h5", ".png")), dpi=130)
    with open(args.out, "w") as g:
        g.write("\n".join(md))
    print(f"wrote {args.out} + figs in {args.figdir}")


def selftest():
    path = "/tmp/selftest_meas1.h5"
    A0t, A1t, tau_t, slew_t = 2.0e-4, 6.0e-5, 25.0, 0.30
    center = 2.082e9
    with h5py.File(path, "w") as g:
        fr = center + np.linspace(-15e6, 15e6, 1001)
        for d in (-8.0, 0.0, 5.0, 20.0, 60.0, 120.0, 200.0):
            for j in range(16):
                ph = j * 2 * np.pi / 16
                base = np.abs(np.random.normal(5e-8, 1e-8, len(fr)))
                win = np.abs(fr - center) <= 0.75e6
                sig = base.copy()
                sig[win] += (A0t + A1t * np.exp(-max(d, 0) / tau_t)
                             * np.cos(ph + slew_t * d)) / win.sum()
                key = f"LGM_c0_d{d:+.3f}_p{ph:.4f}_t0"
                g.create_dataset(key + "__SIG/trace", data=np.stack([fr, sig])[None])
                g.create_dataset(key + "__BL/trace", data=np.stack([fr, base])[None])
    res = parse_file(path)
    rows = fit_config(res["delays"][0])
    tf, sl = tau_fit(rows), slew_fit(rows)
    ok = (abs(rows[0, 1] - A0t) < 0.3e-4 and abs(tf[1] - tau_t) < 3.0
          and abs(np.degrees(sl[0]) - np.degrees(slew_t)) < 3.0)
    print(f"selftest A0 {rows[0,1]*1e6:.1f} uW (want {A0t*1e6:.0f}), "
          f"tau {tf[1]:.1f} ns (want {tau_t}), slew {np.degrees(sl[0]):.2f} deg/ns "
          f"(want {np.degrees(slew_t):.2f}) -> {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
