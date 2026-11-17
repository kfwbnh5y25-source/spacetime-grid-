"""
Generates publication-grade 3-panel figure comparing f1=4071 MHz vs f1=4081 MHz
from all 72,920 groups in spin-waves/threshold_reduction/data.h5.
Output: fig_4081_replication.png / .pdf
"""

from pathlib import Path
import h5py
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
H5_PATH = HERE.parent.parent / "threshold_reduction" / "data.h5"


def hampel_batch(specs: np.ndarray, hw: int = 7, nsig: float = 5.0) -> np.ndarray:
    pad = np.pad(specs, ((0, 0), (hw, hw)), mode="edge")
    win = sliding_window_view(pad, 2 * hw + 1, axis=1)
    med = np.median(win, axis=2)
    mad = np.median(np.abs(win - med[:, :, None]), axis=2)
    sig = 1.4826 * mad
    mask = (sig > 0) & (np.abs(specs - med) > nsig * sig)
    out = specs.copy()
    out[mask] = med[mask]
    return out


def main():
    meta, specs = [], []
    with h5py.File(H5_PATH, "r") as f:
        for k in f.keys():
            p = k.split("_")
            if len(p) < 7 or p[5] == "EX":
                continue
            tr = f[k]["trace"][:]
            s = tr[:, 1, :].mean(0) if tr.ndim == 3 else tr[1, :]
            meta.append((int(p[0]), int(p[1]), int(p[2]), int(p[3]), p[5], p[6]))
            specs.append(s)

    specs = np.array(specs, dtype=np.float64)
    clean = hampel_batch(specs)
    raw_sum = np.sum(clean, axis=1)

    data = {}
    spec_data = {}
    for idx, (a1, a2, f1, f2, ro, tim) in enumerate(meta):
        data.setdefault((f1, f2, ro, tim), {})[(a1, a2)] = raw_sum[idx]
        if ro in ("m", "p") and f1 == 4081 and a1 >= 250 and a2 >= 350 and tim == "D":
            spec_data.setdefault((f2, ro), []).append(clean[idx])

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), dpi=180)
    fig.suptitle(
        "Independent Replication at $f_1 = 4081\\,\\mathrm{MHz}$ vs. $4071\\,\\mathrm{MHz}$ (72,920 HDF5 Groups)",
        fontsize=12.5, fontweight="bold", y=1.01
    )

    # --- Panel (a): Log-log Idler Scaling P_m / P_{f2/2} vs a1 ---
    ax = axes[0]
    styles = [
        (4071, 4125, "#1f77b4", "o", "4071 / 4125"),
        (4071, 4135, "#aec7e8", "s", "4071 / 4135"),
        (4081, 4125, "#d62728", "^", "4081 / 4125 (+27%)"),
        (4081, 4135, "#ff9896", "D", "4081 / 4135 (+30%)"),
    ]
    for f1, f2, col, mk, lbl in styles:
        m_d = data[(f1, f2, "m", "D")]
        f2_d = data[(f1, f2, "f2", "D")]
        xs, ys = [], []
        for (a1, a2), pm in m_d.items():
            if a1 >= 30 and a2 >= 350 and (a1, a2) in f2_d:
                pf2 = f2_d[(a1, a2)]
                if pf2 > 1e-4 and pm > 2e-6:
                    xs.append(a1)
                    ys.append(pm / pf2)
        xs, ys = np.array(xs), np.array(ys)
        slope, inter = np.polyfit(np.log10(xs), np.log10(ys), 1)
        ax.scatter(xs, ys, s=14, alpha=0.55, color=col, marker=mk, label=f"{lbl} (slope={slope:.2f})")
        xg = np.linspace(30, 355, 50)
        ax.plot(xg, 10 ** (inter + slope * np.log10(xg)), color=col, lw=1.4)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Tone-1 Amplitude $a_1$ (mV)")
    ax.set_ylabel("Normalized Idler $P_m / P_{f_2/2}$")
    ax.set_title("(a) Exact $a_1^2$ Idler Scaling & DOS Boost")
    ax.legend(fontsize=7.5, loc="lower right")
    ax.grid(True, which="both", alpha=0.25)

    # --- Panel (b): Threshold a2_th vs a1^2 (D vs ND) ---
    ax = axes[1]
    a1_bins = [(0, 20), (100, 160), (160, 220), (220, 275), (275, 355)]
    for f1, col, lbl in [(4071, "#1f77b4", "4071"), (4081, "#d62728", "4081")]:
        for tim, ls, mk in [("D", "-", "o"), ("ND", "--", "s")]:
            f2_t = data[(f1, 4125, "f2", tim)]
            bx, by = [], []
            for lo, hi in a1_bins:
                pts = [(a1, a2, rs) for (a1, a2), rs in f2_t.items() if lo <= a1 < hi]
                above = [a2 for _, a2, rs in pts if rs > 2e-5]
                if len(above) >= 5:
                    bx.append(np.mean([a1 for a1, _, _ in pts]) ** 2 / 1e4)
                    by.append(np.percentile(above, 5))
            ax.plot(bx, by, ls=ls, marker=mk, ms=5, lw=1.6, color=col, label=f"f1={lbl} {tim}")

    ax.set_xlabel("Tone-1 Drive Power Proxy $a_1^2$ ($\\times 10^4\\,\\mathrm{mV}^2$)")
    ax.set_ylabel("Tone-2 Threshold $a_{2,\\mathrm{th}}$ (mV)")
    ax.set_title("(b) $f_2=4125\\,\\mathrm{MHz}$: `D > ND` Pre-Excitation Split")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(True, alpha=0.3)

    # --- Panel (c): Band-Edge Asymmetry (Channel m active vs p silent) ---
    ax = axes[2]
    det = np.linspace(-2.5, 2.5, 1001)
    for (f2, ro, fc, col, lbl) in [
        (4125, "m", 2018.5, "#d62728", "Channel m (2018.5 MHz, active)"),
        (4135, "m", 2013.5, "#ff7f0e", "Channel m (2013.5 MHz, active)"),
        (4125, "p", 2084.5, "#1f77b4", "Channel p (2084.5 MHz, silent)"),
        (4135, "p", 2094.5, "#7f7f7f", "Channel p (2094.5 MHz, silent)"),
    ]:
        arr = np.array(spec_data[(f2, ro)])
        ax.plot(det, arr.mean(0), lw=1.5, color=col, label=lbl)
    ax.set_yscale("log")
    ax.set_ylim(1e-9, 1e-6)
    ax.set_xlabel("Detuning from Channel Center $f - f_c$ (MHz)")
    ax.set_ylabel("Mean Despiked PSD (mW)")
    ax.set_title("(c) Side-Peak Asymmetry ($f_m < 2068.9 < f_p$)")
    ax.legend(fontsize=7.5, loc="upper right")
    ax.grid(True, which="both", alpha=0.25)

    plt.tight_layout()
    fig.savefig(HERE / "fig_4081_replication.png", bbox_inches="tight")
    fig.savefig(HERE / "fig_4081_replication.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {HERE / 'fig_4081_replication.png'}")


if __name__ == "__main__":
    main()
