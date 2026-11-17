"""
V3 Analysis & Publication Figure Generator for Spoken Digit Recognition:
  - 4-Class Classification Accuracy (Digits 0, 1, 2, 3) via 5-Fold Stratified CV
  - Leave-One-Speaker-Out Cross-Validation (Held-out human speaker generalization)
  - Balanced 0.3x Sub-Threshold Control (True Chance: 25.0%)
  - 500-Iteration Permutation Significance Test
  - 4x4 Confusion Matrix
  - Pilot Drift Verification (CV %)
  - Publication-Quality Figure: fig_spoken_digits_v3.{png,pdf}
"""

import json
import re
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import RidgeClassifier
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT_FIG = ROOT / "spin-waves" / "reports" / "neuromorphic" / "figures"
OUT_FIG.mkdir(parents=True, exist_ok=True)
H5_PATH = HERE / "data_spoken_digits_v3.h5"

DIGIT_NAMES = ["0 (zero)", "1 (one)", "2 (two)", "3 (three)"]


def main():
    if not H5_PATH.exists():
        print(f"[FAIL] Data file {H5_PATH} does not exist yet. Run run_spoken_digits_benchmark_v3.py first.")
        return

    with h5py.File(H5_PATH, "r") as f:
        sample_k = [k for k in f.keys() if "__SIG" in k][0]
        freq = f[sample_k]["trace"][0, 0, :] / 1e6
        mask = (freq >= 2030) & (freq <= 2075)

        def load_block(prefix: str):
            X, y, speakers = [], [], []
            keys = sorted([k for k in f.keys() if k.startswith(prefix) and k.endswith("__SIG")])
            for k in keys:
                m = re.search(r"_d(\d+)_([a-z]+)_r(\d+)", k)
                if not m:
                    continue
                d, spk, rep = m.groups()
                bl_k = k.replace("__SIG", "__BL")
                if bl_k not in f:
                    continue
                s = f[k]["trace"][:, 1, :].mean(0) - f[bl_k]["trace"][:, 1, :].mean(0)
                X.append(s[mask])
                y.append(int(d))
                speakers.append(spk)
            return np.array(X), np.array(y), np.array(speakers)

        X_res, y_res, spk_res = load_block("DIGIT_V3_d")
        X_ctrl, y_ctrl, _ = load_block("DIGIT_V3_CTRL03X_d")

        # Pilot traces for drift tracking
        pilot_keys = sorted(
            [k for k in f.keys() if k.startswith("PILOT_") and k.endswith("__SIG")],
            key=lambda x: int(x.split("_")[1])
        )
        pilot_pows = []
        for k in pilot_keys:
            bl_k = k.replace("__SIG", "__BL")
            s = f[k]["trace"][:, 1, :].mean(0) - f[bl_k]["trace"][:, 1, :].mean(0)
            pilot_pows.append(float(np.sum(s)))

    pilot_cv = float(np.std(pilot_pows) / np.mean(pilot_pows) * 100) if pilot_pows else 0.0

    print("=" * 70)
    print("V3 SPOKEN DIGIT RECOGNITION BENCHMARK RESULTS")
    print("=" * 70)
    print(f"Loaded: {len(y_res)} Reservoir trials across {len(set(spk_res))} speakers")
    print(f"Loaded: {len(y_ctrl)} Balanced sub-threshold 0.3x control trials")

    # 1. 5-Fold Stratified Cross-Validation (Standard Benchmark)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    preds_res = np.zeros(len(y_res), dtype=int)
    for tr, te in skf.split(X_res, y_res):
        sc = StandardScaler().fit(X_res[tr])
        clf = RidgeClassifier(alpha=10.0).fit(sc.transform(X_res[tr]), y_res[tr])
        preds_res[te] = clf.predict(sc.transform(X_res[te]))
    acc_res = accuracy_score(y_res, preds_res)
    cm_res = confusion_matrix(y_res, preds_res)

    # 2. Balanced Control Evaluation (5 trials per digit, Chance = 25.0%)
    if len(y_ctrl) >= 4:
        skf_c = StratifiedKFold(n_splits=min(5, len(y_ctrl)//4), shuffle=True, random_state=42)
        preds_ctrl = np.zeros(len(y_ctrl), dtype=int)
        for tr, te in skf_c.split(X_ctrl, y_ctrl):
            sc = StandardScaler().fit(X_ctrl[tr])
            clf = RidgeClassifier(alpha=10.0).fit(sc.transform(X_ctrl[tr]), y_ctrl[tr])
            preds_ctrl[te] = clf.predict(sc.transform(X_ctrl[te]))
        acc_ctrl = accuracy_score(y_ctrl, preds_ctrl)
    else:
        acc_ctrl = 0.25

    # 3. Leave-One-Speaker-Out Cross-Validation (Cross-Speaker Generalization)
    speakers = sorted(set(spk_res))
    loso_accs = {}
    for held_out in speakers:
        tr_idx = np.where(spk_res != held_out)[0]
        te_idx = np.where(spk_res == held_out)[0]
        sc = StandardScaler().fit(X_res[tr_idx])
        clf = RidgeClassifier(alpha=10.0).fit(sc.transform(X_res[tr_idx]), y_res[tr_idx])
        pred_ho = clf.predict(sc.transform(X_res[te_idx]))
        loso_accs[held_out] = float(accuracy_score(y_res[te_idx], pred_ho))
    mean_loso = float(np.mean(list(loso_accs.values())))

    # 4. Permutation Test (500 Permutations)
    rng = np.random.default_rng(42)
    null_accs = []
    folds = list(skf.split(X_res, y_res))
    for _ in range(500):
        yp = rng.permutation(y_res)
        pp = np.zeros(len(yp), dtype=int)
        for tr, te in folds:
            sc = StandardScaler().fit(X_res[tr])
            clf = RidgeClassifier(alpha=10.0).fit(sc.transform(X_res[tr]), yp[tr])
            pp[te] = clf.predict(sc.transform(X_res[te]))
        null_accs.append(accuracy_score(yp, pp))
    pval = float(np.mean(np.array(null_accs) >= acc_res))

    print(f"\n[PHYSICAL RESERVOIR V3 ACCURACY]")
    print(f"  5-Fold CV Accuracy:            {acc_res * 100:.1f}% (Chance: 25.0%)")
    print(f"  Cross-Speaker Generalization:  {mean_loso * 100:.1f}% (Held-out human speaker)")
    print(f"  Permutation Significance:      p = {pval:.5f}")
    print(f"  Balanced 0.3x Control Acc:     {acc_ctrl * 100:.1f}% (Chance: 25.0%)")
    print(f"  Drift Stability (Pilots):      CV = {pilot_cv:.2f}%\n")
    print("Confusion Matrix:\n", cm_res)

    out_metrics = {
        "dataset": str(H5_PATH),
        "n_trials": len(y_res),
        "accuracy_5fold_cv": acc_res,
        "cross_speaker_generalization": mean_loso,
        "loso_per_speaker": loso_accs,
        "control_03x_acc": acc_ctrl,
        "permutation_pval": pval,
        "pilot_drift_cv": pilot_cv,
        "confusion_matrix": cm_res.tolist(),
    }
    (HERE / "spoken_digits_v3_metrics.json").write_text(json.dumps(out_metrics, indent=2))
    print(f"Metrics saved to {HERE / 'spoken_digits_v3_metrics.json'}")

    # --- Plot Publication-Quality Figure ---
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.5), dpi=180)

    # (a) Confusion Matrix
    ax = axes[0]
    cm_norm = cm_res.astype(float) / cm_res.sum(axis=1, keepdims=True)
    im = ax.imshow(cm_norm, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(4)); ax.set_yticks(range(4))
    ax.set_xticklabels(DIGIT_NAMES, fontsize=9)
    ax.set_yticklabels(DIGIT_NAMES, fontsize=9)
    ax.set_xlabel("Predicted Digit", fontsize=10.5)
    ax.set_ylabel("True Spoken Digit", fontsize=10.5)
    ax.set_title(f"(a) 4-Class Confusion Matrix (Acc={acc_res*100:.1f}%)", fontsize=11, fontweight="bold")
    for i in range(4):
        for j in range(4):
            ax.text(j, i, f"{cm_res[i, j]}\n({cm_norm[i, j]*100:.0f}%)", ha="center", va="center",
                    color="white" if cm_norm[i, j] > 0.5 else "black", fontsize=9, fontweight="bold")

    # (b) Reservoir vs Controls & Cross-Speaker Generalization
    ax = axes[1]
    b_labels = ["5-Fold CV\n(Reservoir)", "Cross-Speaker\n(Held-out)", "Balanced 0.3x\n(Control)"]
    b_vals = [acc_res * 100, mean_loso * 100, acc_ctrl * 100]
    colors = ["#2ca02c", "#1f77b4", "#7f7f7f"]
    ax.bar(b_labels, b_vals, color=colors, edgecolor="k", width=0.55)
    ax.axhline(25.0, color="r", linestyle="--", lw=1.2, label="Chance Level (25.0%)")
    ax.set_ylabel("Classification Accuracy (%)", fontsize=10.5)
    ax.set_ylim(0, 105)
    ax.set_title("(b) Generalization & Balanced Controls", fontsize=11, fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(fontsize=9, loc="upper right")
    for idx, v in enumerate(b_vals):
        ax.text(idx, v + 2.0, f"{v:.1f}%", ha="center", fontweight="bold", fontsize=9.5)

    # (c) Permutation Significance Null Distribution
    ax = axes[2]
    ax.hist(np.array(null_accs) * 100, bins=20, color="#aec7e8", edgecolor="k", alpha=0.8, density=True)
    ax.axvline(acc_res * 100, color="#d62728", lw=2.2, label=f"Physical Reservoir: {acc_res*100:.1f}%\n($p = {pval:.4f}$)")
    ax.axvline(25.0, color="k", linestyle="--", lw=1.0, label="Chance Level (25.0%)")
    ax.set_xlabel("Classification Accuracy (%)", fontsize=10.5)
    ax.set_ylabel("Probability Density", fontsize=10.5)
    ax.set_title("(c) Permutation Significance Test", fontsize=11, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.5, loc="upper right")

    plt.tight_layout()
    fig.savefig(OUT_FIG / "fig_spoken_digits_v3.png", bbox_inches="tight")
    fig.savefig(OUT_FIG / "fig_spoken_digits_v3.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure to {OUT_FIG / 'fig_spoken_digits_v3.png'}")


if __name__ == "__main__":
    main()
