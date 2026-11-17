"""
Generates the 4-figure publication suite for MANUSCRIPT_DRAFT.md in:
  spin-waves/reports/neuromorphic/figures/fig{1,2,3,4}_*.{png,pdf}
Uses pre-computed JSON audits + raw HDF5s where needed for full reproducibility.
"""

import json
from pathlib import Path
import h5py
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
NM_DIR = HERE.parent.parent / "neuromorphic"
IDLER_DIR = HERE.parent / "idler_math"
H5_TR = HERE.parent.parent / "threshold_reduction" / "data.h5"
FIG_DIR = HERE / "figures"
FIG_DIR.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.labelsize": 9.5,
    "legend.fontsize": 8.0,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
})


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


def make_fig1():
    params = json.loads((IDLER_DIR / "idler_extracted_params.json").read_text())
    # Load raw (a1, a2) idler scaling & channel p spectra from threshold_reduction/data.h5
    sums, p_sup, p_sub = {}, [], []
    with h5py.File(H5_TR, "r") as f:
        for k in f.keys():
            p = k.split("_")
            if len(p) < 7 or p[5] == "EX" or p[6] != "D":
                continue
            a1, a2, f1, f2, ro = int(p[0]), int(p[1]), int(p[2]), int(p[3]), p[5]
            if ro in ("m", "f2"):
                tr = f[k]["trace"][:]
                s = tr[:, 1, :].mean(0) if tr.ndim == 3 else tr[1, :]
                sums.setdefault((f1, f2, ro), {})[(a1, a2)] = float(np.sum(hampel_batch(s[None, :])[0]))
            elif ro == "p" and f1 == 4071 and f2 == 4125 and a2 >= 350:
                tr = f[k]["trace"][:]
                s = tr[:, 1, :].mean(0) if tr.ndim == 3 else tr[1, :]
                if a1 >= 275:
                    p_sup.append(s)
                elif 50 <= a1 < 220:
                    p_sub.append(s)

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.4), dpi=200)

    # (a) Log-log idler scaling
    ax = axes[0]
    for f1, f2, col, mk, lbl in [
        (4071, 4125, "#1f77b4", "o", "4071 / 4125"),
        (4071, 4135, "#aec7e8", "s", "4071 / 4135"),
        (4081, 4125, "#d62728", "^", "4081 / 4125 (+29%)"),
        (4081, 4135, "#ff9896", "D", "4081 / 4135 (+28%)"),
    ]:
        m_d, f2_d = sums[(f1, f2, "m")], sums[(f1, f2, "f2")]
        xs, ys = [], []
        for (a1, a2), pm in m_d.items():
            if a1 >= 35 and a2 >= 350 and (a1, a2) in f2_d and f2_d[(a1, a2)] > 1e-4 and pm > 2e-6:
                xs.append(a1); ys.append(pm / f2_d[(a1, a2)])
        xs, ys = np.array(xs), np.array(ys)
        sl, ic = np.polyfit(np.log10(xs), np.log10(ys), 1)
        ax.scatter(xs, ys, s=12, alpha=0.5, color=col, marker=mk, label=f"{lbl} (s={sl:.2f})")
        xg = np.linspace(35, 355, 40)
        ax.plot(xg, 10 ** (ic + sl * np.log10(xg)), color=col, lw=1.4)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Tone-1 Amplitude $a_1$ (mV)"); ax.set_ylabel("Normalized Idler $P_m / P_{f_2/2}$")
    ax.set_title("(a) Sub-Threshold Idler Law ($P_m \\propto a_1^2 P_{f_2/2}$)")
    ax.legend(loc="lower right"); ax.grid(True, which="both", alpha=0.25)

    # (b) Threshold reduction ND vs D
    ax = axes[1]
    for f1, col in [(4071, "#1f77b4"), (4081, "#d62728")]:
        pr = params["pairs"][f"{f1}_4125"]
        tc_nd, tc_d = pr["th_curves"]["ND"], pr["th_curves"]["D"]
        bx_nd = np.array(tc_nd["a1sq"]) / 1e4
        ax.errorbar(bx_nd, tc_nd["a2th"], yerr=tc_nd["se"], fmt="s", ms=4.5, color=col, capsize=3,
                    label=f"f1={f1} ND ($R^2={pr['nd_fit_r2']:.2f}$)")
        xg = np.linspace(0, bx_nd.max(), 40)
        ax.plot(xg, pr["a2_th_0_mv"] + pr["a2_th_slope_mv_per_mv2"] * xg * 1e4, "--", lw=1.5, color=col)
        bx_d = np.array(tc_d["a1sq"]) / 1e4
        ax.errorbar(bx_d, tc_d["a2th"], yerr=tc_d["se"], fmt="o-", ms=4.5, lw=1.1, color=col, mfc="white",
                    label=f"f1={f1} D (4B memory)")
    ax.set_xlabel("Tone-1 Power Proxy $a_1^2$ ($\\times 10^4\\,\\mathrm{mV}^2$)")
    ax.set_ylabel("Tone-2 Threshold $a_{2,\\mathrm{th}}$ (mV)")
    ax.set_title("(b) Two-Regime Threshold Reduction (`ND` vs. `D`)")
    ax.legend(loc="lower left"); ax.grid(True, alpha=0.3)

    # (c) Weak upper side peak p
    ax = axes[2]
    det = np.linspace(-2.5, 2.5, 1001)
    k_sm = np.ones(40) / 40.0
    mu_sup = np.convolve(hampel_batch(np.array(p_sup)).mean(0), k_sm, mode="same") * 1e9
    mu_sub = np.convolve(hampel_batch(np.array(p_sub)).mean(0), k_sm, mode="same") * 1e9
    ax.plot(det[30:-30], mu_sup[30:-30], lw=1.8, color="#1f77b4", label="Channel p, $a_1 > a_{1,\\mathrm{th}}$ ($t=+11.0$)")
    ax.plot(det[30:-30], mu_sub[30:-30], lw=1.2, ls="--", color="#7f7f7f", label="Channel p, $a_1 < a_{1,\\mathrm{th}}$ ($t=+0.8$)")
    ax.set_xlabel("Detuning $f - f_p$ (MHz)"); ax.set_ylabel("Smoothed Excess PSD (pW)")
    ax.set_title("(c) Weak Upper Side Peak (Turns On Above $a_{1,\\mathrm{th}}$)")
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig1_apparatus_and_physics.png", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig1_apparatus_and_physics.pdf", bbox_inches="tight")
    plt.close(fig)


def make_fig2():
    m_e04 = json.loads((NM_DIR / "e0_e4_metrics.json").read_text())
    # Re-run quick load of E0 pilots and E1FIX
    with h5py.File(NM_DIR / "data_e0_e4.h5", "r") as f:
        pks = sorted([k for k in f.keys() if k.startswith("E0_PILOT_") and k.endswith("__SIG")])
        t0 = None
        ts, ps = [], []
        for sk in pks:
            t_raw = f[sk]["timestamp"][0, 0]
            t_val = np.datetime64(t_raw.decode() if isinstance(t_raw, bytes) else str(t_raw)).astype("datetime64[s]").astype(float)
            if t0 is None:
                t0 = t_val
            s = f[sk]["trace"][:, 1, :].mean(0)
            b = f[sk[:-5] + "__BL"]["trace"][:, 1, :].mean(0)
            ts.append((t_val - t0) / 3600.0)
            ps.append(float(np.sum(np.maximum(hampel_batch((s - b)[None, :])[0], 0))))
    order = np.argsort(ts)
    ts, ps = np.array(ts)[order], np.array(ps)[order]
    ps /= np.median(ps)

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.4), dpi=200)

    # (a) Pilot stability
    ax = axes[0]
    ax.plot(ts, ps, "o-", ms=2.2, lw=0.9, color="#1f77b4")
    ax.axhline(1.0, color="k", ls="--", lw=0.8)
    ax.set_ylim(0.8, 1.25)
    ax.set_xlabel("Session Time (hours)"); ax.set_ylabel("Normalized Pilot Power")
    ax.set_title(f"(a) E0: 3.1-Hour Pilot Stability (CV={m_e04['E0']['pilot_cv_raw']*100:.1f}%)")
    ax.grid(True, alpha=0.3)

    # (b) E1FIX cyclic collapse
    ax = axes[1]
    gs = [0, 500, 2000]
    snrs = [0.490, 22.305, 13.733]
    accs = [23.3, 100.0, 100.0]
    ax.bar(["G=0 ns\n(Cyclic)", "G=500 ns", "G=2000 ns"], snrs, color=["#7f7f7f", "#d62728", "#d62728"], edgecolor="k")
    ax.axhline(1.0, color="k", ls="--", lw=0.9, label="Noise Floor (1.0)")
    for i, (s, a) in enumerate(zip(snrs, accs)):
        ax.text(i, s + 0.8, f"Acc={a:.0f}%", ha="center", fontsize=8.5, fontweight="bold")
    ax.set_ylim(0, 26)
    ax.set_ylabel("Cyclic Distance / Noise Floor")
    ax.set_title("(b) E1FIX: Integer-Cycle $G\\to 0$ Collapse Control")
    ax.legend(); ax.grid(True, axis="y", alpha=0.3)

    # (c) E2 vs E4 controls
    ax = axes[2]
    labels = ["E2: ovl", "E2: seq_d0", "E4: randph", "Total-E\n(Static)", "E4: offres", "E4: sub03x"]
    vals = [
        m_e04["E2"]["ovl"]["acc"] * 100,
        m_e04["E2"]["seq_d0"]["acc"] * 100,
        m_e04["E4"]["randph"]["acc"] * 100,
        m_e04["E2"]["ovl"]["total_energy_only_acc"] * 100,
        m_e04["E4"]["offres"]["acc"] * 100,
        m_e04["E4"]["sub03x"]["acc"] * 100,
    ]
    cols = ["#2ca02c", "#2ca02c", "#1f77b4", "#7f7f7f", "#7f7f7f", "#7f7f7f"]
    ax.bar(range(len(labels)), vals, color=cols, edgecolor="k")
    ax.axhline(50.0, color="k", ls="--", lw=0.9, label="Chance (50%)")
    ax.axhline(80.0, color="#2ca02c", ls=":", lw=1.1, label="Gate G1 (80%)")
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    ax.set_ylim(0, 112)
    ax.set_ylabel("5-Fold CV Accuracy (%)")
    ax.set_title("(c) Gate G1 (`AB` vs `BA`) vs. Artifact Controls")
    ax.legend(loc="upper right"); ax.grid(True, axis="y", alpha=0.3)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig2_disproof_and_gate_g1.png", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig2_disproof_and_gate_g1.pdf", bbox_inches="tight")
    plt.close(fig)


def make_fig3():
    params = json.loads((IDLER_DIR / "idler_extracted_params.json").read_text())
    me = params["memory_e3fix"]

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.4), dpi=200)

    # (a) Accuracy vs delay
    ax = axes[0]
    ds = [0, 10, 20, 40, 80, 150, 300]
    acc350 = [100.0, 100.0, 100.0, 43.3, 66.7, 53.3, 40.0]
    acc380 = [100.0, 100.0, 100.0, 63.3, 70.0, 56.7, 46.7]
    ax.plot(ds, acc350, "o-", color="#1f77b4", lw=1.5, label="350 mV")
    ax.plot(ds, acc380, "s-", color="#2ca02c", lw=1.5, label="380 mV")
    ax.axhline(50.0, color="k", ls="--", lw=0.9, label="Chance (50%)")
    ax.set_ylim(30, 108)
    ax.set_xlabel("Inter-Symbol Delay $d$ (ns)"); ax.set_ylabel("5-Fold CV Accuracy (%)")
    ax.set_title("(a) E3FIX: Inter-Symbol Memory Curve (`[A][d][B]`)")
    ax.legend(); ax.grid(True, alpha=0.3)

    # (b) Fisher J decay
    ax = axes[1]
    for amp, col in [("350", "#1f77b4"), ("380", "#2ca02c")]:
        ax.semilogy(me[amp]["delays_ns"], me[amp]["fisher_j"], "o-", lw=1.5, color=col,
                    label=f"{amp} mV ($\\tau_J={me[amp]['tau_J_ns']:.1f}$ ns)")
    ax.set_xlabel("Inter-Symbol Delay $d$ (ns)"); ax.set_ylabel("Order Fisher Separability $J$")
    ax.set_title("(b) Spin-Wave Relaxation Decay of Fisher $J$")
    ax.legend(); ax.grid(True, which="both", alpha=0.3)

    # (c) E5 Frozen Night-1 -> Session 2
    ax = axes[2]
    ax.bar(["Night 1 CV\n(n=320)", "Session 2 Frozen\n(n=80)"], [100.0, 100.0], color=["#1f77b4", "#2ca02c"], edgecolor="k")
    ax.axhline(50.0, color="k", ls="--", lw=0.9)
    ax.axhline(75.0, color="r", ls=":", lw=1.1, label="Gate E5 (75%)")
    ax.set_ylim(0, 112)
    ax.set_ylabel("Classification Accuracy (%)")
    ax.set_title("(c) E5: Frozen Readout Cross-Session Test (100.0%)")
    ax.legend(); ax.grid(True, axis="y", alpha=0.3)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig3_memory_curve_and_session2.png", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig3_memory_curve_and_session2.pdf", bbox_inches="tight")
    plt.close(fig)


def make_fig4():
    res_v2 = json.loads((NM_DIR / "multisymbol_v2_metrics.json").read_text())
    res_act = json.loads((NM_DIR / "multisymbol_active_metrics.json").read_text())["active"]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.5), dpi=200)

    # (a) S3 confusion matrix (v2: 92.2%)
    ax = axes[0]
    s3_cls = ["ABC", "ACB", "BAC", "BCA", "CAB", "CBA"]
    cm = np.array(res_v2["S3"]["cm"], dtype=float)
    cm_n = cm / cm.sum(axis=1, keepdims=True)
    ax.imshow(cm_n, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(6)); ax.set_yticks(range(6))
    ax.set_xticklabels(s3_cls, fontsize=8); ax.set_yticklabels(s3_cls, fontsize=8)
    for i in range(6):
        for j in range(6):
            ax.text(j, i, f"{int(cm[i,j])}", ha="center", va="center",
                    color="white" if cm_n[i,j] > 0.5 else "black", fontsize=8)
    ax.set_xlabel("Predicted"); ax.set_ylabel("True")
    ax.set_title(f"(a) 6-Class $S_3$ Permutations (Acc={res_v2['S3']['acc']*100:.1f}%)")

    # (b) Multi-symbol vs 0.3x controls (v2)
    ax = axes[1]
    b_names = ["S3 (6-cls)", "BIT3 (8-way)", "Matched-1hot", "Matched-2hot"]
    b_keys = ["S3", "BIT3_8way", "BIT3_1hot3", "BIT3_2hot3"]
    b_ctrl = ["S3_03x", "BIT3_03x_8way", "BIT3_03x_1hot3", "BIT3_03x_2hot3"]
    chances = [16.7, 12.5, 33.3, 33.3]
    x = np.arange(len(b_names))
    ax.bar(x - 0.18, [res_v2[k]["acc"]*100 for k in b_keys], width=0.35, color="#2ca02c", edgecolor="k", label="350-365 mV")
    ax.bar(x + 0.18, [res_v2[k]["acc"]*100 for k in b_ctrl], width=0.35, color="#7f7f7f", edgecolor="k", label="0.3x Ctrl")
    for xi, ch in zip(x, chances):
        ax.plot([xi - 0.38, xi + 0.38], [ch, ch], "r--", lw=1.2)
    ax.set_xticks(x); ax.set_xticklabels(b_names, rotation=15, ha="right", fontsize=8)
    ax.set_ylim(0, 108)
    ax.set_ylabel("5-Fold CV Accuracy (%)")
    ax.set_title("(b) Multi-Symbol & Matched-Count Order (v2)")
    ax.legend(); ax.grid(True, axis="y", alpha=0.3)

    # (c) Active learning trajectory
    ax = axes[2]
    traj = res_act["traj"]
    js = [t["j"] for t in traj]
    ni = sum(1 for t in traj if t["phase"] == "INIT")
    ax.plot(range(1, ni + 1), js[:ni], "o-", color="#1f77b4", label=f"Sobol Init (mean={res_act['init_mean_j']:.2f})")
    ax.plot(range(ni + 1, len(js) + 1), js[ni:], "s-", color="#d62728", label=f"BO UCB (mean={res_act['bo_mean_j']:.2f})")
    ax.plot(range(1, len(js) + 1), np.maximum.accumulate(js), "k--", lw=1.1, label=f"Best J={res_act['best_j']:.2f}")
    ax.set_xlabel("Evaluation Step"); ax.set_ylabel("Order Fisher Separability $J$")
    ax.set_title("(c) Bayesian Active-Learning Trajectory")
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig4_multisymbol_and_active.png", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig4_multisymbol_and_active.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    make_fig1()
    make_fig2()
    make_fig3()
    make_fig4()
    print("All 4 manuscript figures saved to", FIG_DIR)
