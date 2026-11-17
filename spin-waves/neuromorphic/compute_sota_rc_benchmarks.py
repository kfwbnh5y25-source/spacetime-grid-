"""
SOTA Neuromorphic Physical Reservoir Computing Benchmarks on Empirical Data:
  1. Kernel Rank (KR) — separation property across distinct temporal input classes.
  2. Generalization Rank (GR) — grouping property / noise-subspace rank across trial repetitions.
  3. Separation-to-Generalization Quality (Q = KR / GR & Fisher trace ratio).
  4. Linear Memory Capacity (MC) — lag-dependent recall profile (lag 0, 1, 2) & delay sweep tau_mem.
  5. Information Processing Capacity (IPC) — orthogonal Walsh decomposition (Linear, Quadratic, Parity)
     comparing physical reservoir vs. 0.3x sub-threshold control.

Outputs:
  - spin-waves/neuromorphic/sota_benchmarks.json
  - spin-waves/reports/neuromorphic/figures/fig_sota_benchmarks.png / .pdf
"""

import json
import re
from pathlib import Path
import h5py
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT_FIG = ROOT / "spin-waves" / "reports" / "neuromorphic" / "figures"
OUT_FIG.mkdir(parents=True, exist_ok=True)

H5_MULTI = HERE / "data_multisymbol_v2.h5"
H5_FOLLOW = HERE / "data_followup.h5"
H5_E0E4 = HERE / "data_e0_e4.h5"


def get_freq_mask(tr: np.ndarray, lo: float = 2030.0, hi: float = 2075.0) -> np.ndarray:
    freq = tr[0, 0, :] / 1e6
    return (freq >= lo) & (freq <= hi)


def compute_kr_gr(X_by_class: dict[str, np.ndarray], sv_thresh_pct: float = 1.0):
    """
    Computes Kernel Rank (KR) and Generalization Rank (GR) following Legenstein & Maass (2007).
    - KR: effective rank of class-mean state matrix (separation property).
    - GR: effective rank of pooled within-class trial variations (generalization/noise property).
    """
    classes = sorted(X_by_class.keys())
    means = np.array([X_by_class[c].mean(axis=0) for c in classes])  # (C, F)
    
    # KR from class centers
    u_m, s_m, vh_m = np.linalg.svd(means - means.mean(axis=0), full_matrices=False)
    kr_thresh = (sv_thresh_pct / 100.0) * s_m[0]
    kr = int(np.sum(s_m > kr_thresh))
    pr_kr = float((np.sum(s_m**2))**2 / np.sum(s_m**4)) if np.sum(s_m**4) > 0 else 0.0

    # GR from within-class deviations
    within_devs = []
    for c in classes:
        xc = X_by_class[c]
        within_devs.append(xc - means[classes.index(c)])
    W = np.vstack(within_devs)  # (Total_N, F)
    
    u_w, s_w, vh_w = np.linalg.svd(W, full_matrices=False)
    gr_thresh = (sv_thresh_pct / 100.0) * s_w[0]
    gr = int(np.sum(s_w > gr_thresh))
    pr_gr = float((np.sum(s_w**2))**2 / np.sum(s_w**4)) if np.sum(s_w**4) > 0 else 0.0

    # Fisher trace ratio: Tr(S_B) / Tr(S_W)
    sb_trace = float(np.sum(s_m**2) * (len(W) / len(classes)))
    sw_trace = float(np.sum(s_w**2))
    fisher_ratio = float(sb_trace / (sw_trace + 1e-12))
    
    # Q factor: KR / GR (dimension ratio)
    q_ratio = float(kr / max(1, gr))

    return {
        "num_classes": len(classes),
        "total_trials": len(W),
        "kr": kr,
        "kr_participation_ratio": pr_kr,
        "kr_singular_values": [float(x) for x in s_m[:min(10, len(s_m))]],
        "gr": gr,
        "gr_participation_ratio": pr_gr,
        "gr_singular_values": [float(x) for x in s_w[:min(10, len(s_w))]],
        "q_ratio": q_ratio,
        "fisher_trace_ratio": fisher_ratio,
    }


def compute_walsh_ipc(X: np.ndarray, U_bits: np.ndarray, n_splits: int = 5, seed: int = 42):
    """
    Computes Information Processing Capacity (IPC) using complete orthogonal Walsh basis.
    U_bits shape: (N, 3) in {-1, +1}.
    """
    u1, u2, u3 = U_bits[:, 0], U_bits[:, 1], U_bits[:, 2]
    walsh = {
        "deg1_u1_lag2": u1,
        "deg1_u2_lag1": u2,
        "deg1_u3_lag0": u3,
        "deg2_u1_u2": u1 * u2,
        "deg2_u2_u3": u2 * u3,
        "deg2_u1_u3": u1 * u3,
        "deg3_u1_u2_u3_parity": u1 * u2 * u3,
    }
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    caps = {}
    for name, target in walsh.items():
        preds = np.zeros(len(target))
        for tr, te in kf.split(X, target):
            sc = StandardScaler().fit(X[tr])
            clf = Ridge(alpha=10.0).fit(sc.transform(X[tr]), target[tr])
            preds[te] = clf.predict(sc.transform(X[te]))
        cov = np.cov(preds, target)[0, 1]
        var_p = np.var(preds)
        var_t = np.var(target)
        c = (cov**2) / (var_p * var_t + 1e-12) if (var_p > 1e-12 and cov > 0) else 0.0
        caps[name] = float(np.clip(c, 0.0, 1.0))

    c_deg1 = caps["deg1_u1_lag2"] + caps["deg1_u2_lag1"] + caps["deg1_u3_lag0"]
    c_deg2 = caps["deg2_u1_u2"] + caps["deg2_u2_u3"] + caps["deg2_u1_u3"]
    c_deg3 = caps["deg3_u1_u2_u3_parity"]
    c_tot = c_deg1 + c_deg2 + c_deg3
    return {
        "components": caps,
        "capacity_deg1_linear": float(c_deg1),
        "capacity_deg2_quadratic": float(c_deg2),
        "capacity_deg3_parity": float(c_deg3),
        "total_ipc": float(c_tot),
    }


def analyze_e3fix_memory(h5_follow_path: Path):
    """Analyzes continuous-time inter-symbol delay sweeps in E3FIX."""
    with h5py.File(h5_follow_path, "r") as f:
        e3_keys = [k for k in f.keys() if k.startswith("E3FIX_") and k.endswith("__SIG")]
        if not e3_keys:
            return None
        mask = get_freq_mask(f[e3_keys[0]]["trace"])

        data = {}  # amp -> delay -> order -> list of spectra
        for k in e3_keys:
            m = re.search(r"E3FIX_a(\d+)_([AB]{2})_d(\d+)_r(\d+)", k)
            if not m:
                continue
            amp, order, delay, rep = m.groups()
            delay = int(delay)
            bl_k = k.replace("__SIG", "__BL")
            if bl_k not in f:
                continue
            s = f[k]["trace"][:, 1, :].mean(0) - f[bl_k]["trace"][:, 1, :].mean(0)
            data.setdefault(amp, {}).setdefault(delay, {}).setdefault(order, []).append(s[mask])

    res = {}
    for amp, d_dict in data.items():
        delays = sorted(d_dict.keys())
        j_vals, d_list = [], []
        for d in delays:
            ab = np.array(d_dict[d]["AB"])
            ba = np.array(d_dict[d]["BA"])
            mu_diff = np.mean(ab, axis=0) - np.mean(ba, axis=0)
            var_pool = 0.5 * (np.var(ab, axis=0) + np.var(ba, axis=0)) + 1e-15
            j = float(np.sum(mu_diff**2 / var_pool))
            j_vals.append(j)
            d_list.append(d)

        # Fit exponential decay: J(d) = J_inf + (J_0 - J_inf) * exp(-d / tau)
        d_arr, j_arr = np.array(d_list, dtype=float), np.array(j_vals, dtype=float)
        try:
            popt, _ = curve_fit(
                lambda t, j0, jinf, tau: jinf + (j0 - jinf) * np.exp(-t / tau),
                d_arr, j_arr,
                p0=[j_arr[0], j_arr[-1], 30.0],
                bounds=([0, 0, 1.0], [np.inf, np.inf, 500.0])
            )
            tau_mem = float(popt[2])
        except Exception:
            tau_mem = float("nan")

        res[f"amp_{amp}mV"] = {
            "delays_ns": d_list,
            "fisher_j": j_vals,
            "tau_mem_ns": tau_mem,
            "j0": float(j_arr[0]),
            "j300ns": float(j_arr[-1]),
        }
    return res


def main():
    print("Computing SOTA Neuromorphic Physical Reservoir Benchmarks...")
    metrics = {}

    # 1. Multi-symbol v2 data: BIT3 and S3
    with h5py.File(H5_MULTI, "r") as f:
        sample_k = [k for k in f.keys() if "__SIG" in k][0]
        mask = get_freq_mask(f[sample_k]["trace"])

        # Load BIT3 classes
        bit3_by_cls = {}
        bit3_03x_by_cls = {}
        X_bit3, U_bit3 = [], []
        X_bit3_03x, U_bit3_03x = [], []

        for k in f.keys():
            if not k.endswith("__SIG"):
                continue
            bl_k = k.replace("__SIG", "__BL")
            if bl_k not in f:
                continue

            # BIT3 standard
            m_b = re.search(r"^BIT3_([01]{3})_p\d+_r\d+__SIG$", k)
            if m_b:
                w = m_b.group(1)
                s = f[k]["trace"][:, 1, :].mean(0) - f[bl_k]["trace"][:, 1, :].mean(0)
                bit3_by_cls.setdefault(w, []).append(s[mask])
                X_bit3.append(s[mask])
                U_bit3.append([1 if c == "1" else -1 for c in w])

            # BIT3 0.3x control
            m_b03 = re.search(r"^BIT3_03x_([01]{3})_p\d+_r\d+__SIG$", k)
            if m_b03:
                w = m_b03.group(1)
                s = f[k]["trace"][:, 1, :].mean(0) - f[bl_k]["trace"][:, 1, :].mean(0)
                bit3_03x_by_cls.setdefault(w, []).append(s[mask])
                X_bit3_03x.append(s[mask])
                U_bit3_03x.append([1 if c == "1" else -1 for c in w])

        # Load S3 classes
        s3_by_cls = {}
        s3_03x_by_cls = {}
        for k in f.keys():
            if not k.endswith("__SIG"):
                continue
            bl_k = k.replace("__SIG", "__BL")
            if bl_k not in f:
                continue

            m_s = re.search(r"^S3_([A-Z]{3})_p\d+_r\d+__SIG$", k)
            if m_s:
                w = m_s.group(1)
                s = f[k]["trace"][:, 1, :].mean(0) - f[bl_k]["trace"][:, 1, :].mean(0)
                s3_by_cls.setdefault(w, []).append(s[mask])

            m_s03 = re.search(r"^S3_03x_([A-Z]{3})_p\d+_r\d+__SIG$", k)
            if m_s03:
                w = m_s03.group(1)
                s = f[k]["trace"][:, 1, :].mean(0) - f[bl_k]["trace"][:, 1, :].mean(0)
                s3_03x_by_cls.setdefault(w, []).append(s[mask])

    bit3_by_cls = {k: np.array(v) for k, v in bit3_by_cls.items()}
    bit3_03x_by_cls = {k: np.array(v) for k, v in bit3_03x_by_cls.items()}
    s3_by_cls = {k: np.array(v) for k, v in s3_by_cls.items()}
    s3_03x_by_cls = {k: np.array(v) for k, v in s3_03x_by_cls.items()}

    # Compute KR and GR for BIT3 and S3
    metrics["KR_GR_BIT3"] = compute_kr_gr(bit3_by_cls)
    metrics["KR_GR_BIT3_03x"] = compute_kr_gr(bit3_03x_by_cls)
    metrics["KR_GR_S3"] = compute_kr_gr(s3_by_cls)
    metrics["KR_GR_S3_03x"] = compute_kr_gr(s3_03x_by_cls)

    # Compute IPC and Linear Memory Capacity
    metrics["IPC_BIT3_Reservoir"] = compute_walsh_ipc(np.array(X_bit3), np.array(U_bit3))
    metrics["IPC_BIT3_03x_Control"] = compute_walsh_ipc(np.array(X_bit3_03x), np.array(U_bit3_03x))

    # Compute continuous memory decay from E3FIX
    e3_res = analyze_e3fix_memory(H5_FOLLOW)
    if e3_res:
        metrics["E3FIX_Memory_Decay"] = e3_res

    # Save metrics JSON
    json_path = HERE / "sota_benchmarks.json"
    json_path.write_text(json.dumps(metrics, indent=2))
    print(f"Metrics saved to {json_path}")

    # --- Plot SOTA Benchmarks Summary Figure ---
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), dpi=180)

    # Panel a: Information Processing Capacity decomposition
    ax = axes[0, 0]
    bar_labels = ["Linear\n(MC)", "Quadratic\n(2-body)", "Cubic\n(Parity)", "Total\nIPC"]
    v2_ipc = metrics["IPC_BIT3_Reservoir"]
    c03_ipc = metrics["IPC_BIT3_03x_Control"]
    v2_vals = [v2_ipc["capacity_deg1_linear"], v2_ipc["capacity_deg2_quadratic"], v2_ipc["capacity_deg3_parity"], v2_ipc["total_ipc"]]
    c03_vals = [c03_ipc["capacity_deg1_linear"], c03_ipc["capacity_deg2_quadratic"], c03_ipc["capacity_deg3_parity"], c03_ipc["total_ipc"]]
    
    x = np.arange(len(bar_labels))
    ax.bar(x - 0.18, v2_vals, width=0.35, color="#2ca02c", edgecolor="k", label="Physical Reservoir (365 mV)")
    ax.bar(x + 0.18, c03_vals, width=0.35, color="#7f7f7f", edgecolor="k", label="0.3x Sub-th Control")
    ax.set_xticks(x)
    ax.set_xticklabels(bar_labels, fontsize=10)
    ax.set_ylabel("Capacity $C$ (Dimensionless)", fontsize=11)
    ax.set_title("(a) Information Processing Capacity (IPC Decomposition)", fontsize=11.5, fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(fontsize=9)
    for i, (v, c) in enumerate(zip(v2_vals, c03_vals)):
        ax.text(i - 0.18, v + 0.03, f"{v:.2f}", ha="center", fontsize=8.5, fontweight="bold")
        ax.text(i + 0.18, c + 0.03, f"{c:.2f}", ha="center", fontsize=8.5, color="#444")
    ax.set_ylim(0, 1.9)

    # Panel b: Linear Memory Decay across symbol lags (k=0, 1, 2)
    ax = axes[0, 1]
    lags = [0, 1, 2]
    c_lags_v2 = [v2_ipc["components"]["deg1_u3_lag0"], v2_ipc["components"]["deg1_u2_lag1"], v2_ipc["components"]["deg1_u1_lag2"]]
    c_lags_03 = [c03_ipc["components"]["deg1_u3_lag0"], c03_ipc["components"]["deg1_u2_lag1"], c03_ipc["components"]["deg1_u1_lag2"]]
    ax.plot(lags, c_lags_v2, "o-", color="#1f77b4", lw=2, markersize=8, label=f"Physical Reservoir (Total MC = {sum(c_lags_v2):.2f})")
    ax.plot(lags, c_lags_03, "s--", color="#7f7f7f", lw=1.5, markersize=6, label=f"0.3x Control (Total MC = {sum(c_lags_03):.2f})")
    ax.set_xticks(lags)
    ax.set_xticklabels(["Lag 0\n(symbol $u_3$)", "Lag 1\n(symbol $u_2$)", "Lag 2\n(symbol $u_1$)"], fontsize=9.5)
    ax.set_ylabel("Recall Capacity $R^2(u_{t-k})$", fontsize=11)
    ax.set_title("(b) Discrete Symbol Memory Recall Horizon", fontsize=11.5, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9.5)
    ax.set_ylim(-0.05, 0.8)

    # Panel c: Continuous Time Relaxation Memory Curve (E3FIX Delay Sweeps)
    ax = axes[1, 0]
    if e3_res:
        c_350 = e3_res.get("amp_350mV")
        c_380 = e3_res.get("amp_380mV")
        if c_350:
            ax.plot(c_350["delays_ns"], c_350["fisher_j"], "o-", color="#d62728", lw=1.8, label=f"350 mV ($\\tau={c_350['tau_mem_ns']:.1f}$ ns)")
        if c_380:
            ax.plot(c_380["delays_ns"], c_380["fisher_j"], "s-", color="#9467bd", lw=1.8, label=f"380 mV ($\\tau={c_380['tau_mem_ns']:.1f}$ ns)")
        ax.axhline(0, color="k", lw=0.6, linestyle=":")
        ax.set_xlabel("Inter-Symbol Gap $d$ (ns)", fontsize=11)
        ax.set_ylabel("Linear Fisher Separability $J(AB, BA)$", fontsize=11)
        ax.set_title("(c) Continuous Physical Relaxation Memory Curve", fontsize=11.5, fontweight="bold")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9.5)

    # Panel d: Kernel Rank vs. Generalization Rank (Separation Quality)
    ax = axes[1, 1]
    datasets_plot = ["BIT3 (8 words)", "S3 (6 perms)"]
    kr_vals = [metrics["KR_GR_BIT3"]["kr"], metrics["KR_GR_S3"]["kr"]]
    gr_vals = [metrics["KR_GR_BIT3"]["gr"], metrics["KR_GR_S3"]["gr"]]
    kr_03 = [metrics["KR_GR_BIT3_03x"]["kr"], metrics["KR_GR_S3_03x"]["kr"]]
    gr_03 = [metrics["KR_GR_BIT3_03x"]["gr"], metrics["KR_GR_S3_03x"]["gr"]]

    x_d = np.arange(len(datasets_plot))
    w = 0.2
    ax.bar(x_d - 1.5*w, kr_vals, width=w, color="#1f77b4", edgecolor="k", label="KR (Separation, Reservoir)")
    ax.bar(x_d - 0.5*w, gr_vals, width=w, color="#aec7e8", edgecolor="k", label="GR (Noise Subspace, Reservoir)")
    ax.bar(x_d + 0.5*w, kr_03, width=w, color="#ff7f0e", edgecolor="k", label="KR (0.3x Control)")
    ax.bar(x_d + 1.5*w, gr_03, width=w, color="#ffbb78", edgecolor="k", label="GR (0.3x Control)")
    ax.set_xticks(x_d)
    ax.set_xticklabels(datasets_plot, fontsize=10)
    ax.set_ylabel("Effective Matrix Rank (SVD > 1%)", fontsize=11)
    ax.set_title("(d) Kernel Rank vs. Generalization Rank", fontsize=11.5, fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")

    plt.tight_layout()
    fig.savefig(OUT_FIG / "fig_sota_benchmarks.png", bbox_inches="tight")
    fig.savefig(OUT_FIG / "fig_sota_benchmarks.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure to {OUT_FIG / 'fig_sota_benchmarks.png'}")

    # Print summary
    print("\n" + "="*70)
    print("SOTA NEUROMORPHIC BENCHMARK SUMMARY (NOVEMBER 2026 STANDARDS)")
    print("="*70)
    print(f"1. Linear Memory Capacity (MC):")
    print(f"   - Physical Reservoir: Total MC = {v2_ipc['capacity_deg1_linear']:.3f} symbols (Lag0={v2_ipc['components']['deg1_u3_lag0']:.2f}, Lag1={v2_ipc['components']['deg1_u2_lag1']:.2f}, Lag2={v2_ipc['components']['deg1_u1_lag2']:.2f})")
    print(f"   - 0.3x Sub-th Control: Total MC = {c03_ipc['capacity_deg1_linear']:.3f} symbols")
    print(f"   - Physical Relaxation Time (E3FIX): tau_mem = {e3_res['amp_350mV']['tau_mem_ns']:.1f} ns (350 mV) / {e3_res['amp_380mV']['tau_mem_ns']:.1f} ns (380 mV)")
    print(f"\n2. Information Processing Capacity (IPC Decomposition):")
    print(f"   - Degree 1 (Linear):     {v2_ipc['capacity_deg1_linear']:.3f} (Reservoir) vs. {c03_ipc['capacity_deg1_linear']:.3f} (0.3x)")
    print(f"   - Degree 2 (Quadratic):  {v2_ipc['capacity_deg2_quadratic']:.3f} (Reservoir) vs. {c03_ipc['capacity_deg2_quadratic']:.3f} (0.3x)")
    print(f"   - Degree 3 (Parity):     {v2_ipc['capacity_deg3_parity']:.3f} (Reservoir) vs. {c03_ipc['capacity_deg3_parity']:.3f} (0.3x)")
    print(f"   - Total IPC:             {v2_ipc['total_ipc']:.3f} (Reservoir) vs. {c03_ipc['total_ipc']:.3f} (0.3x)")
    print(f"\n3. Kernel Rank (KR) & Generalization Rank (GR):")
    print(f"   - BIT3 (8 words): KR = {metrics['KR_GR_BIT3']['kr']}, GR = {metrics['KR_GR_BIT3']['gr']}, Fisher ratio = {metrics['KR_GR_BIT3']['fisher_trace_ratio']:.2f}")
    print(f"   - S3 (6 perms):   KR = {metrics['KR_GR_S3']['kr']}, GR = {metrics['KR_GR_S3']['gr']}, Fisher ratio = {metrics['KR_GR_S3']['fisher_trace_ratio']:.2f}")
    print("="*70)


if __name__ == "__main__":
    main()
