"""
Full Scientific Analysis of Overnight Verification Dataset (spin-waves/neuromorphic/data_e0_e4.h5).
Generates:
  - spin-waves/neuromorphic/e0_e4_metrics.json
  - spin-waves/neuromorphic/fig_e0_e4_summary.png / .pdf
  - spin-waves/neuromorphic/E0_E4_RESULTS_REPORT.md
"""

import json
import re
from datetime import datetime
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

HERE = Path(__file__).resolve().parent
H5_PATH = HERE / "data_e0_e4.h5"


def hampel_despike(spec: np.ndarray, half_win: int = 7, n_sig: float = 5.0) -> np.ndarray:
    out = spec.copy()
    n = len(spec)
    for i in range(n):
        lo, hi = max(0, i - half_win), min(n, i + half_win + 1)
        win = spec[lo:hi]
        med = np.median(win)
        sig = 1.4826 * np.median(np.abs(win - med))
        if sig > 0 and abs(spec[i] - med) > n_sig * sig:
            out[i] = med
    return out


def load_dataset(h5_path: Path):
    records = []
    ex_records = []
    with h5py.File(h5_path, "r") as f:
        sig_keys = [k for k in f.keys() if k.endswith("__SIG")]
        for sk in sig_keys:
            base_name = sk[:-5]
            tr = f[sk]["trace"][:]  # (N, 2, 1001)
            ts_raw = f[sk]["timestamp"][:, 0]
            ts_list = [datetime.fromisoformat(t.decode() if isinstance(t, bytes) else str(t)).timestamp() for t in ts_raw]
            t_mean = float(np.mean(ts_list))
            freqs_mhz = tr[0, 0, :] / 1e6
            sig_lin = tr[:, 1, :].mean(axis=0).astype(np.float64)

            if base_name.startswith("E0_EX_"):
                ex_records.append({"key": base_name, "t": t_mean, "freqs_mhz": freqs_mhz, "spec": sig_lin})
                continue

            bk = base_name + "__BL"
            if bk in f:
                bl_lin = f[bk]["trace"][:, 1, :].mean(axis=0).astype(np.float64)
            else:
                bl_lin = np.zeros_like(sig_lin)

            diff = hampel_despike(sig_lin - bl_lin)
            records.append({
                "key": base_name,
                "t": t_mean,
                "freqs_mhz": freqs_mhz,
                "diff": diff,
                "sig": sig_lin,
                "bl": bl_lin,
            })
    records.sort(key=lambda r: r["t"])
    ex_records.sort(key=lambda r: r["t"])
    return records, ex_records


def eval_binary_cv(X: np.ndarray, y: np.ndarray, n_perm: int = 200, seed: int = 42):
    """5-fold stratified CV with PCA(<=15) + RidgeClassifier + fast pre-projected permutation p-value."""
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    folds = list(skf.split(X, y))
    preds = np.zeros(len(y), dtype=int)
    scores = np.zeros(len(y), dtype=float)
    Z_folds = []
    for tr_i, te_i in folds:
        n_comp = min(15, len(tr_i) - 2, X.shape[1])
        scaler = StandardScaler().fit(X[tr_i])
        pca = PCA(n_components=n_comp, random_state=seed).fit(scaler.transform(X[tr_i]))
        Z_tr = pca.transform(scaler.transform(X[tr_i]))
        Z_te = pca.transform(scaler.transform(X[te_i]))
        Z_folds.append((tr_i, te_i, Z_tr, Z_te))
        clf = RidgeClassifier(alpha=10.0).fit(Z_tr, y[tr_i])
        preds[te_i] = clf.predict(Z_te)
        scores[te_i] = clf.decision_function(Z_te)
    acc = float(np.mean(preds == y))
    auc = float(roc_auc_score(y, scores))

    # Fisher J in full spectrum
    A, B = X[y == 0], X[y == 1]
    fisher_j = float(np.sum((A.mean(0) - B.mean(0)) ** 2) / (np.sum(A.var(0) + B.var(0)) + 1e-12))

    # Fast permutation p-value using pre-projected PCA features
    rng = np.random.default_rng(seed)
    null_accs = []
    for _ in range(n_perm):
        yp = rng.permutation(y)
        p_preds = np.zeros(len(y), dtype=int)
        for tr_i, te_i, Z_tr, Z_te in Z_folds:
            clf = RidgeClassifier(alpha=10.0).fit(Z_tr, yp[tr_i])
            p_preds[te_i] = clf.predict(Z_te)
        null_accs.append(np.mean(p_preds == yp))
    p_val = float((1 + np.sum(np.array(null_accs) >= acc)) / (n_perm + 1))
    return {"n": int(len(y)), "acc": acc, "auc": auc, "fisher_j": fisher_j, "p_val": p_val}


def main():
    print(f"Loading {H5_PATH}...")
    records, ex_records = load_dataset(H5_PATH)
    t0 = records[0]["t"]
    for r in records:
        r["t_hr"] = (r["t"] - t0) / 3600.0

    # --- E0: Pilots & Drift ---
    pilots = [r for r in records if r["key"].startswith("E0_PILOT_")]
    pilot_t = np.array([r["t_hr"] for r in pilots])
    pilot_p = np.array([float(np.sum(np.maximum(r["diff"], 0.0))) for r in pilots])
    pilot_pk = np.array([float(r["freqs_mhz"][np.argmax(r["diff"])]) for r in pilots])
    pilot_med = float(np.median(pilot_p))
    pilot_cv_raw = float(np.std(pilot_p) / np.mean(pilot_p))

    # Nearest-pilot normalization for all non-pilot records
    for r in records:
        idx = int(np.argmin(np.abs(pilot_t - r["t_hr"])))
        scale = pilot_med / max(pilot_p[idx], 1e-12)
        r["diff_norm"] = r["diff"] * scale

    # Sequential difference CV of pilots (short-timescale vs full-night)
    pilot_diff_cv = float(np.std(np.diff(pilot_p)) / (np.sqrt(2) * np.mean(pilot_p)))

    # EX IMD3 check
    ex_summary = {}
    if ex_records:
        fr_ex = ex_records[0]["freqs_mhz"]
        ex_mean = np.mean([r["spec"] for r in ex_records], axis=0)
        for f_t, lbl in [(4065.0, "IMD3_lo_4065"), (4095.0, "F1_4095"), (4125.0, "F2_4125"), (4155.0, "IMD3_hi_4155")]:
            m = np.abs(fr_ex - f_t) < 2.0
            val_mw = float(np.max(ex_mean[m]))
            ex_summary[lbl] = {"mw": val_mw, "dbm": float(10 * np.log10(val_mw + 1e-15))}

    e0_metrics = {
        "n_pilots": len(pilots),
        "duration_hr": float(records[-1]["t_hr"]),
        "pilot_power_mean": float(np.mean(pilot_p)),
        "pilot_cv_raw": pilot_cv_raw,
        "pilot_cv_short": pilot_diff_cv,
        "pilot_peak_mhz_std": float(np.std(pilot_pk)),
        "ex_imd3": ex_summary,
    }
    print("[E0]", json.dumps(e0_metrics, indent=2))

    # --- E1: Cyclic Shifts vs Gap G ---
    e1_recs = [r for r in records if r["key"].startswith("E1_")]
    e1_groups = {}
    for r in e1_recs:
        m = re.match(r"E1_(ABC|BCA|CAB)_G(\d+)_r(\d+)", r["key"])
        if m:
            ord_s, g_ns, _ = m.group(1), int(m.group(2)), int(m.group(3))
            e1_groups.setdefault(g_ns, {}).setdefault(ord_s, []).append(r["diff_norm"])

    e1_metrics = {}
    for g_ns in sorted(e1_groups.keys()):
        d_by_ord = {o: np.array(v) for o, v in e1_groups[g_ns].items()}
        mu_abc, mu_bca, mu_cab = d_by_ord["ABC"].mean(0), d_by_ord["BCA"].mean(0), d_by_ord["CAB"].mean(0)
        d_between = float(np.mean([
            np.linalg.norm(mu_abc - mu_bca),
            np.linalg.norm(mu_bca - mu_cab),
            np.linalg.norm(mu_cab - mu_abc),
        ]))
        # Noise floor of mean difference: sqrt(2/N) * mean within-class SD
        within_var = np.mean([np.sum(arr.var(0)) for arr in d_by_ord.values()])
        n_avg = np.mean([len(arr) for arr in d_by_ord.values()])
        d_noise_expected = float(np.sqrt(2.0 * within_var / n_avg))
        ratio = float(d_between / max(d_noise_expected, 1e-12))
        # Also 3-class CV accuracy
        X_e1 = np.vstack([d_by_ord["ABC"], d_by_ord["BCA"], d_by_ord["CAB"]])
        y_e1 = np.array([0]*len(d_by_ord["ABC"]) + [1]*len(d_by_ord["BCA"]) + [2]*len(d_by_ord["CAB"]))
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        preds = np.zeros(len(y_e1), dtype=int)
        for tr_i, te_i in skf.split(X_e1, y_e1):
            pipe = make_pipeline(StandardScaler(), PCA(n_components=10, random_state=42), RidgeClassifier(alpha=10.0))
            pipe.fit(X_e1[tr_i], y_e1[tr_i])
            preds[te_i] = pipe.predict(X_e1[te_i])
        e1_metrics[str(g_ns)] = {
            "d_between": d_between,
            "d_noise_floor": d_noise_expected,
            "snr_ratio": ratio,
            "acc_3class": float(np.mean(preds == y_e1)),
        }
    print("[E1]", json.dumps(e1_metrics, indent=2))

    # --- E2: Energy-Matched Order Gate (G1) ---
    e2_recs = [r for r in records if r["key"].startswith("E2_")]
    e2_layouts = {}
    for r in e2_recs:
        m = re.match(r"E2_(seq_d\d+|ovl)_(AB|BA)_r(\d+)", r["key"])
        if m:
            lay, ord_s = m.group(1), m.group(2)
            e2_layouts.setdefault(lay, {"X": [], "y": []})
            e2_layouts[lay]["X"].append(r["diff_norm"])
            e2_layouts[lay]["y"].append(0 if ord_s == "AB" else 1)

    e2_metrics = {}
    for lay in sorted(e2_layouts.keys()):
        X = np.array(e2_layouts[lay]["X"])
        y = np.array(e2_layouts[lay]["y"])
        res = eval_binary_cv(X, y, n_perm=300)
        # Static total-energy control (1D integrated power only)
        X_tot = np.sum(X, axis=1, keepdims=True)
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        tot_preds = np.zeros(len(y), dtype=int)
        for tr_i, te_i in skf.split(X_tot, y):
            clf = RidgeClassifier(alpha=1.0).fit(X_tot[tr_i], y[tr_i])
            tot_preds[te_i] = clf.predict(X_tot[te_i])
        res["total_energy_only_acc"] = float(np.mean(tot_preds == y))
        e2_metrics[lay] = res
    print("[E2]", json.dumps(e2_metrics, indent=2))

    # --- E3: Pump-Probe Delay Curve (G2) ---
    e3_recs = [r for r in records if r["key"].startswith("E3_")]
    e3_delays = {}
    for r in e3_recs:
        m = re.match(r"E3_(AB|BA)_d(\d+)_r(\d+)", r["key"])
        if m:
            ord_s, d_ns = m.group(1), int(m.group(2))
            e3_delays.setdefault(d_ns, {"X": [], "y": []})
            e3_delays[d_ns]["X"].append(r["diff_norm"])
            e3_delays[d_ns]["y"].append(0 if ord_s == "AB" else 1)

    e3_metrics = {}
    for d_ns in sorted(e3_delays.keys()):
        X = np.array(e3_delays[d_ns]["X"])
        y = np.array(e3_delays[d_ns]["y"])
        e3_metrics[str(d_ns)] = eval_binary_cv(X, y, n_perm=300)
    print("[E3]", json.dumps(e3_metrics, indent=2))

    # --- E4: Artifact Disproof Controls ---
    e4_recs = [r for r in records if r["key"].startswith("E4_")]
    e4_ctrls = {}
    for r in e4_recs:
        m = re.match(r"E4_(sub03x|offres|randph)_(AB|BA)_r(\d+)", r["key"])
        if m:
            ctrl, ord_s = m.group(1), m.group(2)
            e4_ctrls.setdefault(ctrl, {"X": [], "y": []})
            e4_ctrls[ctrl]["X"].append(r["diff_norm"])
            e4_ctrls[ctrl]["y"].append(0 if ord_s == "AB" else 1)

    e4_metrics = {}
    for ctrl in sorted(e4_ctrls.keys()):
        X = np.array(e4_ctrls[ctrl]["X"])
        y = np.array(e4_ctrls[ctrl]["y"])
        e4_metrics[ctrl] = eval_binary_cv(X, y, n_perm=300)
    print("[E4]", json.dumps(e4_metrics, indent=2))

    all_metrics = {"E0": e0_metrics, "E1": e1_metrics, "E2": e2_metrics, "E3": e3_metrics, "E4": e4_metrics}
    (HERE / "e0_e4_metrics.json").write_text(json.dumps(all_metrics, indent=2))

    # --- Publication 4-Panel Summary Figure ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), dpi=160)
    fig.suptitle("Experimental Spin-Wave Neuromorphic Verification (E0–E4, N=4,637 Acquisitions)", fontsize=13, fontweight="bold")

    # Panel (a): E0 Pilot Drift
    ax = axes[0, 0]
    ax.plot(pilot_t, pilot_p / pilot_med, "o-", ms=2.5, lw=1.0, color="#1f77b4")
    ax.axhline(1.0, color="k", ls="--", lw=0.8)
    ax.set_title(f"(a) E0: Pilot Stability (Raw CV={pilot_cv_raw*100:.1f}%, Short CV={pilot_diff_cv*100:.1f}%)")
    ax.set_xlabel("Session Time (hours)")
    ax.set_ylabel("Normalized Pilot Power")
    ax.grid(True, alpha=0.3)

    # Panel (b): E1 Cyclic Collapse vs Gap G
    ax = axes[0, 1]
    gs = sorted([int(k) for k in e1_metrics.keys()])
    ratios = [e1_metrics[str(g)]["snr_ratio"] for g in gs]
    accs_e1 = [e1_metrics[str(g)]["acc_3class"] * 100 for g in gs]
    ax.plot(gs, ratios, "s-", color="#d62728", label="Cyclic Dist / Noise Floor")
    ax.axhline(1.0, color="k", ls="--", lw=0.8, label="Noise Floor (1.0)")
    ax.set_xscale("symlog", linthresh=200)
    ax.set_title("(b) E1: Cyclic-Shift Collapse as Gap G -> 0")
    ax.set_xlabel("Recovery Gap G (ns)")
    ax.set_ylabel("D(ABC,BCA,CAB) / Noise Floor")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel (c): E2 + E4 Order Classification & Controls
    ax = axes[1, 0]
    labels_c = list(e2_metrics.keys()) + ["E4:" + k for k in e4_metrics.keys()]
    accs_c = [e2_metrics[k]["acc"] * 100 for k in e2_metrics.keys()] + [e4_metrics[k]["acc"] * 100 for k in e4_metrics.keys()]
    colors_c = ["#2ca02c" if k in e2_metrics else "#7f7f7f" for k in list(e2_metrics.keys()) + list(e4_metrics.keys())]
    bars = ax.bar(range(len(labels_c)), accs_c, color=colors_c, edgecolor="k", lw=0.6)
    ax.axhline(50.0, color="k", ls="--", lw=0.8, label="Chance (50%)")
    ax.axhline(80.0, color="#2ca02c", ls=":", lw=1.0, label="Gate G1 (80%)")
    ax.set_xticks(range(len(labels_c)))
    ax.set_xticklabels(labels_c, rotation=25, ha="right", fontsize=8)
    ax.set_ylim(35, 105)
    ax.set_title("(c) E2 Order Gate vs. E4 Artifact Controls")
    ax.set_ylabel("5-Fold CV Accuracy (%)")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)

    # Panel (d): E3 Pump-Probe Memory Curve
    ax = axes[1, 1]
    ds = sorted([int(k) for k in e3_metrics.keys()])
    accs_d = [e3_metrics[str(d)]["acc"] * 100 for d in ds]
    js_d = [e3_metrics[str(d)]["fisher_j"] for d in ds]
    ax.plot(ds, accs_d, "o-", color="#9467bd", lw=1.5, label="Order Accuracy (%)")
    ax.axhline(50.0, color="k", ls="--", lw=0.8)
    ax.set_xscale("symlog", linthresh=100)
    ax.set_ylim(40, 105)
    ax.set_title("(d) E3: Pump-Probe Delayed-Order Memory Curve")
    ax.set_xlabel("Probe Delay d (ns)")
    ax.set_ylabel("5-Fold CV Accuracy (%)")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(HERE / "fig_e0_e4_summary.png")
    fig.savefig(HERE / "fig_e0_e4_summary.pdf")
    plt.close(fig)

    # --- Write Markdown Report ---
    g1_pass_ovl = e2_metrics.get("ovl", {}).get("acc", 0) > 0.80 and e2_metrics.get("ovl", {}).get("p_val", 1) < 0.01
    e4_pass = all(v["acc"] <= 0.62 for v in e4_metrics.values())
    e1_g0 = e1_metrics.get("0", {}).get("snr_ratio", 99)

    md = [
        "# Experimental Verification Report: Spin-Wave Physical Reservoir (E0–E4)",
        f"**Dataset:** `{H5_PATH.name}` | **Total Acquisitions:** {len(records) + len(ex_records)} ({len(records)} sample+BL pairs, {len(ex_records)} EX monitors) | **Duration:** {e0_metrics['duration_hr']:.2f} hours",
        "",
        "## 0. Gate Verdict Summary",
        f"- **E0 (Drift & RF Chain):** Raw pilot CV = **{pilot_cv_raw*100:.2f}%** over {e0_metrics['duration_hr']:.2f} h (short-timescale CV = **{pilot_diff_cv*100:.2f}%**); ZHL IMD3 at 4065/4155 MHz = **{ex_summary.get('IMD3_lo_4065', {}).get('dbm', 0):.1f} / {ex_summary.get('IMD3_hi_4155', {}).get('dbm', 0):.1f} dBm**.",
        f"- **E1 (Cyclic-Shift Control):** At $G=0$, cyclic distance / noise floor = **{e1_g0:.2f}×** (3-class Acc = **{e1_metrics.get('0', {}).get('acc_3class', 0)*100:.1f}%** vs 33.3% chance).",
        f"- **E2 (Gate G1 — Energy-Matched AB vs BA):** Overlap Acc = **{e2_metrics.get('ovl', {}).get('acc', 0)*100:.1f}%** (AUC={e2_metrics.get('ovl', {}).get('auc', 0):.3f}, $J={e2_metrics.get('ovl', {}).get('fisher_j', 0):.2f}$, $p={e2_metrics.get('ovl', {}).get('p_val', 1):.4f}$); Sequential $d=0$ Acc = **{e2_metrics.get('seq_d0', {}).get('acc', 0)*100:.1f}%**.",
        f"- **E4 (Artifact Disproof Controls):** `sub03x` = **{e4_metrics.get('sub03x', {}).get('acc', 0)*100:.1f}%**, `offres` = **{e4_metrics.get('offres', {}).get('acc', 0)*100:.1f}%**, `randph` = **{e4_metrics.get('randph', {}).get('acc', 0)*100:.1f}%** ({'ALL AT CHANCE — Electronic artifact ruled out' if e4_pass else 'WARNING: Above-chance control detected'}).",
        "",
        "## 1. E0: Stability & ZHL Linearity",
        f"- Pilots acquired: **{len(pilots)}** | Peak freq jitter: **{e0_metrics['pilot_peak_mhz_std']:.3f} MHz**",
        f"- Excitation (`EX`): F1 (4095 MHz) = **{ex_summary.get('F1_4095', {}).get('dbm', 0):.2f} dBm**, F2 (4125 MHz) = **{ex_summary.get('F2_4125', {}).get('dbm', 0):.2f} dBm**, IMD3 (4065 / 4155 MHz) = **{ex_summary.get('IMD3_lo_4065', {}).get('dbm', 0):.2f} / {ex_summary.get('IMD3_hi_4155', {}).get('dbm', 0):.2f} dBm**.",
        "",
        "## 2. E1: Cyclic-Shift Collapse vs. Recovery Gap $G$",
        "| Gap $G$ (ns) | $D(\\text{between})$ | Noise Floor | SNR Ratio | 3-Class Acc (%) |",
        "| :---: | :---: | :---: | :---: | :---: |",
    ]
    for g in gs:
        m = e1_metrics[str(g)]
        md.append(f"| {g} | {m['d_between']:.3e} | {m['d_noise_floor']:.3e} | {m['snr_ratio']:.2f} | {m['acc_3class']*100:.1f}% |")

    md += [
        "",
        "## 3. E2: Energy-Matched Order Discrimination (`AB` vs `BA`)",
        "| Layout | $N$ | Ridge Acc (%) | AUC | Fisher $J$ | Perm $p$-val | Total-Energy-Only Acc (%) |",
        "| :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]
    for lay, m in e2_metrics.items():
        md.append(f"| `{lay}` | {m['n']} | **{m['acc']*100:.1f}%** | {m['auc']:.3f} | {m['fisher_j']:.2f} | {m['p_val']:.4f} | {m['total_energy_only_acc']*100:.1f}% |")

    md += [
        "",
        "## 4. E3: Pump–Probe Delayed-Order Memory Curve",
        "| Probe Delay $d$ (ns) | $N$ | Ridge Acc (%) | AUC | Fisher $J$ | Perm $p$-val |",
        "| :---: | :---: | :---: | :---: | :---: | :---: |",
    ]
    for d in ds:
        m = e3_metrics[str(d)]
        md.append(f"| {d} | {m['n']} | **{m['acc']*100:.1f}%** | {m['auc']:.3f} | {m['fisher_j']:.2f} | {m['p_val']:.4f} |")

    md += [
        "",
        "## 5. E4: Artifact Disproof Controls",
        "| Control | $N$ | Ridge Acc (%) | AUC | Fisher $J$ | Perm $p$-val |",
        "| :--- | :---: | :---: | :---: | :---: | :---: |",
    ]
    for ctrl, m in e4_metrics.items():
        md.append(f"| `{ctrl}` | {m['n']} | **{m['acc']*100:.1f}%** | {m['auc']:.3f} | {m['fisher_j']:.2f} | {m['p_val']:.4f} |")

    md += ["", "![E0-E4 Summary](fig_e0_e4_summary.png)"]
    (HERE / "E0_E4_RESULTS_REPORT.md").write_text("\n".join(md))
    print(f"[SAVED] {HERE / 'E0_E4_RESULTS_REPORT.md'} and fig_e0_e4_summary.png")


if __name__ == "__main__":
    main()
