"""
Offline Evaluation & Figure Generator for Unified High-Accuracy Benchmarks:
  - Task 1: Tuned Nonlinear Channel Equalization (EQ-Plus)
  - Task 2: Single-Speaker Spoken Digit Recognition (Jackson 4-Class)

Outputs:
  - spin-waves/neuromorphic/high_accuracy_metrics2.json
  - spin-waves/reports/neuromorphic/figures/fig_high_accuracy.png / .pdf
"""

import matplotlib.pyplot as plt
import argparse
import json
import re
from pathlib import Path
import h5py
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
import matplotlib
matplotlib.use("Agg")

HERE = Path(__file__).resolve().parent
DEFAULT_H5 = HERE / "data_high_accuracy.h5"
REPORT_FIGS = HERE.parent / "reports" / "neuromorphic" / "figures"
REPORT_FIGS.mkdir(parents=True, exist_ok=True)


def hampel_filter(specs: np.ndarray, hw: int = 7, nsig: float = 5.0) -> np.ndarray:
    """Removes single-bin spur noise from spectrum analyzer traces."""
    pad = np.pad(specs, ((0, 0), (hw, hw)), mode="edge")
    win = sliding_window_view(pad, 2 * hw + 1, axis=1)
    med = np.median(win, axis=2)
    mad = np.median(np.abs(win - med[:, :, None]), axis=2)
    sig = 1.4826 * mad
    mask = (sig > 0) & (np.abs(specs - med) > nsig * sig)
    out = specs.copy()
    out[mask] = med[mask]
    return out


def evaluate_cv(X: np.ndarray, y: np.ndarray, n_splits: int = 5, alpha: float = 1.0, n_perm: int = 200, seed: int = 42):
    """Evaluates cross-validation accuracy and permutation p-value using pure NumPy Ridge."""
    classes = np.unique(y)
    K = len(classes)
    N = len(y)
    actual_splits = min(n_splits, min(np.sum(y == c) for c in classes))

    rng = np.random.default_rng(seed)
    folds = [[] for _ in range(actual_splits)]
    for c in classes:
        c_idx = rng.permutation(np.where(y == c)[0])
        for idx, item in enumerate(c_idx):
            folds[idx % actual_splits].append(item)

    Y_onehot = np.zeros((N, K))
    for i, c in enumerate(y):
        Y_onehot[i, np.where(classes == c)[0][0]] = 1.0

    preds = np.zeros(N, dtype=int)
    for k in range(actual_splits):
        te_idx = np.array(folds[k])
        tr_idx = np.array([i for i in range(N) if i not in te_idx])
        X_tr, Y_tr = X[tr_idx], Y_onehot[tr_idx]
        X_te = X[te_idx]

        mu = X_tr.mean(axis=0)
        std = X_tr.std(axis=0) + 1e-8
        Z_tr = (X_tr - mu) / std
        Z_te = (X_te - mu) / std

        d = Z_tr.shape[1]
        try:
            W = np.linalg.solve(Z_tr.T @ Z_tr + alpha *
                                np.eye(d), Z_tr.T @ Y_tr)
        except np.linalg.LinAlgError:
            W = np.linalg.pinv(Z_tr.T @ Z_tr + alpha *
                               np.eye(d)) @ (Z_tr.T @ Y_tr)

        logits = Z_te @ W
        preds[te_idx] = classes[np.argmax(logits, axis=1)]

    acc = float(np.mean(preds == y))

    # Permutation test
    nulls = []
    for _ in range(n_perm):
        yp = rng.permutation(y)
        Yp_onehot = np.zeros((N, K))
        for i, c in enumerate(yp):
            Yp_onehot[i, np.where(classes == c)[0][0]] = 1.0

        p_preds = np.zeros(N, dtype=int)
        for k in range(actual_splits):
            te_idx = np.array(folds[k])
            tr_idx = np.array([i for i in range(N) if i not in te_idx])
            X_tr, Y_tr = X[tr_idx], Yp_onehot[tr_idx]
            X_te = X[te_idx]
            mu = X_tr.mean(axis=0)
            std = X_tr.std(axis=0) + 1e-8
            Z_tr = (X_tr - mu) / std
            Z_te = (X_te - mu) / std
            d = Z_tr.shape[1]
            try:
                W = np.linalg.solve(Z_tr.T @ Z_tr + alpha *
                                    np.eye(d), Z_tr.T @ Y_tr)
            except np.linalg.LinAlgError:
                W = np.linalg.pinv(Z_tr.T @ Z_tr + alpha *
                                   np.eye(d)) @ (Z_tr.T @ Y_tr)
            p_preds[te_idx] = classes[np.argmax(Z_te @ W, axis=1)]
        nulls.append(np.mean(p_preds == yp))

    p_val = float((1 + np.sum(np.array(nulls) >= acc)) / (n_perm + 1))

    # Confusion matrix
    cm = np.zeros((K, K), dtype=int)
    for true_c, pred_c in zip(y, preds):
        i = np.where(classes == true_c)[0][0]
        j = np.where(classes == pred_c)[0][0]
        cm[i, j] += 1

    return {"acc": acc, "p_val": p_val, "cm": cm.tolist(), "classes": [int(c) for c in classes], "n": N}


def analyze_dataset(h5_path: Path):
    dryrun_json = Path(str(h5_path) + ".dryrun.json")
    if not h5_path.exists() and not dryrun_json.exists():
        raise FileNotFoundError(f"Dataset not found: {h5_path}")

    keys, specs, meta_dict = [], [], {}
    if h5_path.exists():
        print(f"Loading HDF5 dataset: {h5_path.name}")
        with h5py.File(h5_path, "r") as f:
            freqs = None
            for k in sorted(f.keys()):
                if not k.endswith("__SIG"):
                    continue
                base_k = k[:-5]
                tr = f[k]["trace"][:]
                if freqs is None:
                    freqs = tr[0, 0, :] / 1e6
                s = tr[:, 1, :].mean(axis=0).astype(np.float64)
                bl_k = base_k + "__BL"
                b = f[bl_k]["trace"][:, 1, :].mean(axis=0).astype(
                    np.float64) if bl_k in f else np.zeros_like(s)
                keys.append(base_k)
                specs.append(s - b)
                meta_dict[base_k] = json.loads(
                    f[k].attrs.get("meta_json", "{}"))
    else:
        print(f"Loading dryrun JSON dataset: {dryrun_json.name}")
        content = json.loads(dryrun_json.read_text())
        manifest = content.get("sample", []) if isinstance(
            content, dict) else content
        freqs = np.linspace(2000.0, 2080.0, 1001)
        rng = np.random.default_rng(42)
        for entry in manifest:
            base_k = entry["key"]
            keys.append(base_k)
            specs.append(rng.normal(0, 1, size=1001))
            meta_dict[base_k] = entry.get("meta", entry)

    clean_specs = hampel_filter(np.array(specs))
    band_mask = (freqs >= 2030.0) & (freqs <= 2075.0)
    data = {k: v[band_mask] for k, v in zip(keys, clean_specs)}

    results = {}

    # 1. Channel Equalization
    eq_X, eq_y = [], []
    eq_ctrl_X, eq_ctrl_y = [], []
    for k, v in data.items():
        if k.startswith("EQ_TUNED"):
            eq_X.append(v)
            eq_y.append(meta_dict[k]["target_lag1"])
        elif k.startswith("EQ_CTRL03X"):
            eq_ctrl_X.append(v)
            eq_ctrl_y.append(meta_dict[k]["target_lag1"])

    if eq_X:
        results["eq_lag1"] = evaluate_cv(np.array(eq_X), np.array(eq_y))
        print(
            f"[EQ-Plus Lag 1] Accuracy: {results['eq_lag1']['acc']*100:.1f}% (p={results['eq_lag1']['p_val']:.4f}, N={results['eq_lag1']['n']})")
    if eq_ctrl_X:
        results["eq_ctrl03x"] = evaluate_cv(
            np.array(eq_ctrl_X), np.array(eq_ctrl_y))
        print(
            f"[EQ 0.3x Control] Accuracy: {results['eq_ctrl03x']['acc']*100:.1f}% (Chance: 50.0%)")

    # 2. Spoken Digits (Jackson 4-Class)
    d_X, d_y = [], []
    d_ctrl_X, d_ctrl_y = [], []
    for k, v in data.items():
        if k.startswith("DIGIT_JACKSON_d"):
            d_X.append(v)
            d_y.append(meta_dict[k]["digit"])
        elif k.startswith("DIGIT_JACKSON_CTRL03X"):
            d_ctrl_X.append(v)
            d_ctrl_y.append(meta_dict[k]["digit"])

    if d_X:
        results["digits_jackson"] = evaluate_cv(np.array(d_X), np.array(d_y))
        print(
            f"[Jackson Spoken Digits] 4-Class Accuracy: {results['digits_jackson']['acc']*100:.1f}% (p={results['digits_jackson']['p_val']:.4f}, N={results['digits_jackson']['n']})")
    if d_ctrl_X:
        results["digits_ctrl03x"] = evaluate_cv(
            np.array(d_ctrl_X), np.array(d_ctrl_y))
        print(
            f"[Digits 0.3x Control] Accuracy: {results['digits_ctrl03x']['acc']*100:.1f}% (Chance: 25.0%)")

    metrics_out = HERE / "high_accuracy_metrics.json"
    metrics_out.write_text(json.dumps(results, indent=2))
    print(f"Saved metrics to: {metrics_out}")

    # Plot 2-Panel Summary Figure
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=200)
    fig.suptitle("High-Accuracy Physical Reservoir Benchmarks (BiYIG Microdisk)",
                 fontsize=13, fontweight="bold")

    # Panel A: Channel Equalization
    ax = axes[0]
    labels = ["Lag 1 (s(t-1))", "0.3x Control"]
    accs = [results.get("eq_lag1", {}).get("acc", 0.0) * 100,
            results.get("eq_ctrl03x", {}).get("acc", 0.5) * 100]
    colors = ["#1f77b4", "#7f7f7f"]
    bars = ax.bar(labels, accs, color=colors, edgecolor="k", width=0.45)
    ax.axhline(50.0, color="r", linestyle="--",
               lw=1.2, label="Chance Level (50%)")
    ax.set_ylim(0, 105)
    ax.set_ylabel("5-Fold Cross-Validation Accuracy (%)", fontsize=10)
    ax.set_title(
        "(a) Tuned Nonlinear Channel Equalization (EQ-Plus)", fontsize=10.5)
    ax.legend(loc="lower right")
    ax.grid(True, axis="y", alpha=0.3)
    for b in bars:
        h = b.get_height()
        ax.text(b.get_x() + b.get_width() / 2.0, h + 2.0,
                f"{h:.1f}%", ha="center", va="bottom", fontweight="bold")

    # Panel B: Spoken Digits Confusion Matrix or Bar
    ax = axes[1]
    if "digits_jackson" in results and "cm" in results["digits_jackson"]:
        cm = np.array(results["digits_jackson"]["cm"], dtype=float)
        cm_norm = cm / np.maximum(cm.sum(axis=1, keepdims=True), 1e-12)
        im = ax.imshow(cm_norm, cmap="Blues", vmin=0, vmax=1)
        ax.set_xticks(range(4))
        ax.set_yticks(range(4))
        ax.set_xticklabels(["'0'", "'1'", "'2'", "'3'"])
        ax.set_yticklabels(["'0'", "'1'", "'2'", "'3'"])
        ax.set_xlabel("Predicted Digit", fontsize=10)
        ax.set_ylabel("True Digit", fontsize=10)
        for i in range(4):
            for j in range(4):
                val = int(cm[i, j])
                pct = cm_norm[i, j] * 100
                ax.text(j, i, f"{val}\n({pct:.0f}%)", ha="center", va="center",
                        color="white" if pct > 50 else "black", fontsize=9, fontweight="bold")
        ax.set_title(
            f"(b) Jackson Spoken Digits (Acc={results['digits_jackson']['acc']*100:.1f}%, 0.3x={results.get('digits_ctrl03x', {}).get('acc', 0.25)*100:.1f}%)", fontsize=10.5)
        fig.colorbar(im, ax=ax, fraction=0.046,
                     pad=0.04, label="Normalized Recall")

    plt.tight_layout()
    fig_png = REPORT_FIGS / "fig_high_accuracy2.png"
    fig_pdf = REPORT_FIGS / "fig_high_accuracy2.pdf"
    fig.savefig(fig_png, bbox_inches="tight")
    fig.savefig(fig_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Generated publication figures: {fig_png}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5", type=str, default=str(DEFAULT_H5),
                        help="Path to data_high_accuracy.h5")
    args = parser.parse_args()
    analyze_dataset(Path(args.h5))


if __name__ == "__main__":
    main()
