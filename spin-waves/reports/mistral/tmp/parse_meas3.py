"""
parse_meas3.py - ANALYSIS-ONLY: repeat-averaged phase-gating depth from meas3 h5.

Keys: DEP_{OP}_p{ph:.4f}_t{t} (SIG/BL). Per OP: for each trial t, cos-fit (c12) of the
line power vs phi -> per-trial depth 1 - Pmin/Pmax from (A0 +/- A1) Pmax/Pmin fit;
then mean +/- std/CI across trials (fresh baseline per block is the artifact kill:
flag trials whose A0 deviates > 20% from the median A0, and recompute excluding).
Also cos residual per block (unstructured = OK). Writes table + figure.
Selftest injects 5 trials with known A1/A0 and one drift-corrupted block.
"""
import argparse
import os
import sys

import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def line_power(fr, sig, bl, center_hz, win_mhz=1.5):
    d = np.asarray(sig, float) - np.asarray(bl, float)
    win = np.abs(fr / 1e6 - center_hz / 1e6) <= win_mhz
    noise_bins = ~win & (d < np.quantile(d, 0.9))
    return float(d[win].sum()), float(np.std(d[noise_bins]))


def c12(points):
    ph = np.array([p[0] for p in points]); pw = np.array([p[1] for p in points])
    c1 = np.mean(pw * np.exp(-1j * ph))
    A = np.c_[np.ones_like(ph), np.cos(ph), np.sin(ph)]
    coef, *_ = np.linalg.lstsq(A, pw, rcond=None)
    resid = float(np.std(pw - A @ coef))
    return float(2*np.abs(c1)), float(np.mean(pw)), float(np.angle(c1)), resid


def parse(path, center):
    data = {}
    with h5py.File(path, "r") as f:
        ks = [k for k in f if k.endswith("__SIG") and k.startswith("DEP_")]
        for k in ks:
            name = k[:-5]
            fr = f[k + "/trace"][0, 0]
            sig = f[k + "/trace"][0, 1]
            bl = f[name + "__BL/trace"][0, 1]
            p, n = line_power(fr, sig, bl, center)
            try:
                parts = name.split("_")
                op = parts[1]
                ph = float(parts[2][1:])
                t = int(parts[3][1:])
            except (ValueError, IndexError):
                continue
            data.setdefault(op, {}).setdefault(t, []).append((ph, p))
    out = {}
    for op, trials in data.items():
        per = []
        for t in sorted(trials):
            pts = sorted(trials[t])
            if len(pts) < 8:
                continue
            A1, A0, phi0, resid = c12(pts)
            pmax, pmin = A0 + abs(A1), A0 - abs(A1)
            depth = 1 - pmin / pmax if pmax > 0 else float("nan")
            per.append((t, A0, A1, A1 / max(A0, 1e-30), depth, resid, phi0))
        out[op] = np.array(per) if per else None
    return out


def summarize(op, per):
    A0s = per[:, 1]
    good = np.abs(A0s - np.median(A0s)) < 0.2 * np.median(A0s)
    depths = per[:, 4]
    rows = [f"| {int(t)} | {a0*1e6:.1f} | {a1*1e6:.1f} | {r:.3f} | {d:.3f} | {res*1e6:.1f}"
            f" | {'OK' if g else 'DRIFT?' } |"
            for (t, a0, a1, r, d, res, _) in per
            for g in [bool(good[int(t)])]]
    dd = depths[good] if good.any() else depths
    return rows, float(np.mean(dd)), float(np.std(dd)), int((~good).sum())


def selftest():
    path = "/tmp/selftest_meas3.h5"
    center = 2058.5e9
    A0t, A1t = 2.0e-4, 0.9e-4          # depth = 1-(A0-A1)/(A0+A1) = 0.45
    rng = np.random.default_rng(9)
    with h5py.File(path, "w") as g:
        fr = center + np.linspace(-15e6, 15e6, 1001)
        for t in range(5):
            drift = 1.4 if t == 3 else 1.0     # one drift-corrupted trial
            for j in range(16):
                ph = j * 2 * np.pi / 16
                base = np.abs(rng.normal(5e-8, 1e-8, len(fr))) * drift
                win = np.abs(fr - center) <= 0.75e6
                sig = base.copy()
                sig[win] += (A0t * drift + A1t * drift * np.cos(ph + 0.7)) / win.sum()
                key = f"DEP_OPX_p{ph:.4f}_t{t}"
                g.create_dataset(key + "__SIG/trace", data=np.stack([fr, sig])[None])
                g.create_dataset(key + "__BL/trace", data=np.stack([fr, base])[None])
    res = parse(path, center)
    per = res["OPX"]
    rows, dm, ds, ndr = summarize("OPX", per)
    want = 2 * A1t / (A0t + A1t)      # depth = 1 - Pmin/Pmax = 2A1/(A0+A1)
    ok = abs(dm - want) < 0.02 and ndr == 1
    print(f"selftest: depth {dm:.3f} +/- {ds:.3f} (want {want:.3f}), "
          f"flagged drifts {ndr} (want 1) -> {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--out", default="parse_meas3_out.md")
    ap.add_argument("--figdir", default="parse_meas3_figs")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    os.makedirs(args.figdir, exist_ok=True)
    md = ["# parse_meas3: repeat-averaged depth\n"]
    for path in args.files:
        with h5py.File(path, "r") as f0:
            center = float(np.median(f0[list(f0)[0] + "/trace"][0, 0]))
        res = parse(path, center)
        md.append(f"\n## {os.path.basename(path)}\n")
        for op, per in res.items():
            if per is None:
                continue
            rows, dm, ds, ndr = summarize(op, per)
            md.append(f"\n### {op}: depth = {dm:.3f} +/- {ds:.3f} "
                      f"({len(per)} trials, {ndr} drift-flagged)\n")
            md.append("| t | A0 uW | A1 uW | A1/A0 | depth | resid uW | flag |")
            md.append("|---|---|---|---|---|---|---|")
            md += rows
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(per[:, 0], per[:, 4], "o-")
            ax.axhline(dm, color="r", ls="--", label=f"mean {dm:.3f}")
            ax.set_xlabel("trial"); ax.set_ylabel("depth 1-Pmin/Pmax"); ax.legend()
            fig.tight_layout()
            fig.savefig(os.path.join(args.figdir,
                        f"{os.path.basename(path).replace('.h5','')}_{op}.png"), dpi=130)
    with open(args.out, "w") as g:
        g.write("\n".join(md))
    print(f"wrote {args.out} + figs")


if __name__ == "__main__":
    main()
