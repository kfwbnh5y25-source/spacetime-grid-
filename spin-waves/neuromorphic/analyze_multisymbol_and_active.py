"""
Comprehensive Analysis of:
  1. spin-waves/neuromorphic/data_multisymbol.h5 (6-class S3 permutations + 8-class BIT3 + matched-count order + 3-bit parity + 0.3x controls)
  2. spin-waves/neuromorphic/data_active.h5 (Bayesian Active Learning trajectory + canonical AB/BA readout)

Outputs:
  - multisymbol_active_metrics.json
  - fig_multisymbol_and_active.png / .pdf
"""

import json
import re
from pathlib import Path
import h5py
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent


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


def load_h5(path: Path):
    keys, specs = [], []
    with h5py.File(path, "r") as f:
        for sk in sorted(f.keys()):
            if not sk.endswith("__SIG"):
                continue
            name = sk[:-5]
            s = f[sk]["trace"][:, 1, :].mean(0).astype(np.float64)
            bk = name + "__BL"
            b = f[bk]["trace"][:, 1, :].mean(0).astype(np.float64) if bk in f else np.zeros_like(s)
            keys.append(name)
            specs.append(s - b)
    clean = hampel_batch(np.array(specs)) if specs else np.zeros((0, 1001))
    return dict(zip(keys, clean))


def eval_multiclass_cv(X: np.ndarray, y: np.ndarray, n_comp: int = 20, n_perm: int = 200, seed: int = 42):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    folds = list(skf.split(X, y))
    preds = np.zeros(len(y), dtype=int)
    Z_folds = []
    for tr_i, te_i in folds:
        nc = min(n_comp, len(tr_i) - 2, X.shape[1])
        sc = StandardScaler().fit(X[tr_i])
        pca = PCA(n_components=nc, random_state=seed).fit(sc.transform(X[tr_i]))
        Z_tr, Z_te = pca.transform(sc.transform(X[tr_i])), pca.transform(sc.transform(X[te_i]))
        Z_folds.append((tr_i, te_i, Z_tr, Z_te))
        clf = RidgeClassifier(alpha=10.0).fit(Z_tr, y[tr_i])
        preds[te_i] = clf.predict(Z_te)
    acc = float(np.mean(preds == y))
    cm = confusion_matrix(y, preds)

    rng = np.random.default_rng(seed)
    nulls = []
    for _ in range(n_perm):
        yp = rng.permutation(y)
        pp = np.zeros(len(y), dtype=int)
        for tr_i, te_i, Z_tr, Z_te in Z_folds:
            pp[te_i] = RidgeClassifier(alpha=10.0).fit(Z_tr, yp[tr_i]).predict(Z_te)
        nulls.append(np.mean(pp == yp))
    p_val = float((1 + np.sum(np.array(nulls) >= acc)) / (n_perm + 1))
    return {"n": int(len(y)), "n_classes": int(len(np.unique(y))), "acc": acc, "p_val": p_val, "cm": cm.tolist()}


def main():
    results = {}

    # 1. Analyze data_multisymbol.h5
    ms = load_h5(HERE / "data_multisymbol.h5")
    print(f"Loaded {len(ms)} trials from data_multisymbol.h5")

    # Task A: S3 (6-class) and S3_03x
    s3_classes = ["ABC", "ACB", "BAC", "BCA", "CAB", "CBA"]
    for lbl in ["S3", "S3_03x"]:
        X, y = [], []
        for k, v in ms.items():
            m = re.match(rf"^{lbl}_([ABC]{{3}})_r\d+$", k)
            if m and m.group(1) in s3_classes:
                X.append(v)
                y.append(s3_classes.index(m.group(1)))
        results[lbl] = eval_multiclass_cv(np.array(X), np.array(y), n_comp=20)
        print(f"[{lbl}] 6-class Acc = {results[lbl]['acc']*100:.1f}% (chance 16.7%, p={results[lbl]['p_val']:.4f})")

    # Task B: BIT3 (8-class), matched-count 3-class ({001,010,100} and {011,101,110}), and 3-bit parity
    bit3_words = ["000", "001", "010", "011", "100", "101", "110", "111"]
    for lbl in ["BIT3", "BIT3_03x"]:
        X_all, y_all, y_par = [], [], []
        X_w1, y_w1, X_w2, y_w2 = [], [], [], []
        w1_set = ["001", "010", "100"]
        w2_set = ["011", "101", "110"]
        for k, v in ms.items():
            m = re.match(rf"^{lbl}_([01]{{3}})_r\d+$", k)
            if m:
                w = m.group(1)
                X_all.append(v)
                y_all.append(bit3_words.index(w))
                y_par.append(sum(int(c) for c in w) % 2)
                if w in w1_set:
                    X_w1.append(v); y_w1.append(w1_set.index(w))
                elif w in w2_set:
                    X_w2.append(v); y_w2.append(w2_set.index(w))
        results[f"{lbl}_8way"] = eval_multiclass_cv(np.array(X_all), np.array(y_all), n_comp=20)
        results[f"{lbl}_parity"] = eval_multiclass_cv(np.array(X_all), np.array(y_par), n_comp=20)
        results[f"{lbl}_1hot3"] = eval_multiclass_cv(np.array(X_w1), np.array(y_w1), n_comp=15)
        results[f"{lbl}_2hot3"] = eval_multiclass_cv(np.array(X_w2), np.array(y_w2), n_comp=15)
        print(f"[{lbl}] 8-way={results[f'{lbl}_8way']['acc']*100:.1f}% | Parity={results[f'{lbl}_parity']['acc']*100:.1f}% | "
              f"Matched-1hot={results[f'{lbl}_1hot3']['acc']*100:.1f}% | Matched-2hot={results[f'{lbl}_2hot3']['acc']*100:.1f}%")

    # 2. Analyze data_active.h5
    act = load_h5(HERE / "data_active.h5")
    pil_p = [np.sum(np.maximum(v, 0)) for k, v in act.items() if k.startswith("PILOT_")]
    pil_med = float(np.median(pil_p)) if pil_p else 1.0
    groups = {}
    for k, v in act.items():
        m = re.match(r"^(INIT|BO)_(\d+)_(AB|BA)_r\d+$", k)
        if m:
            groups.setdefault((m.group(1), int(m.group(2))), {"AB": [], "BA": []})[m.group(3)].append(v / pil_med)

    traj = []
    X_act, y_act = [], []
    for (ph, idx), d in sorted(groups.items(), key=lambda x: (0 if x[0][0] == "INIT" else 1, x[0][1])):
        A, B = np.array(d["AB"]), np.array(d["BA"])
        j_lin = float(np.sum((A.mean(0) - B.mean(0))**2) / (np.sum(A.var(0) + B.var(0)) + 1e-12))
        traj.append({"phase": ph, "idx": idx, "j": j_lin})
        X_act.extend(list(A) + list(B))
        y_act.extend([0]*len(A) + [1]*len(B))

    act_cv = eval_multiclass_cv(np.array(X_act), np.array(y_act), n_comp=15)
    bo_summary = json.loads((HERE / "data_active.h5.bo_summary.json").read_text()) if (HERE / "data_active.h5.bo_summary.json").exists() else {}
    results["active"] = {
        "n_configs": len(traj),
        "best_j": max(t["j"] for t in traj),
        "init_mean_j": float(np.mean([t["j"] for t in traj if t["phase"] == "INIT"])),
        "bo_mean_j": float(np.mean([t["j"] for t in traj if t["phase"] == "BO"])),
        "pooled_ab_ba_acc": act_cv["acc"],
        "pooled_p_val": act_cv["p_val"],
        "bo_summary": bo_summary,
        "traj": traj,
    }
    print(f"[ACTIVE] {len(traj)} configs | best J={results['active']['best_j']:.2f} | "
          f"INIT mean={results['active']['init_mean_j']:.2f} -> BO mean={results['active']['bo_mean_j']:.2f} | "
          f"Pooled AB/BA Acc={act_cv['acc']*100:.1f}% (p={act_cv['p_val']:.4f})")

    (HERE / "multisymbol_active_metrics.json").write_text(json.dumps(results, indent=2))

    # --- 4-Panel Figure ---
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.5), dpi=170)
    fig.suptitle("Multi-Symbol Reservoir Computing & Bayesian Active Learning Results", fontsize=13, fontweight="bold")

    # (a) S3 6-class Confusion Matrix
    ax = axes[0, 0]
    cm_s3 = np.array(results["S3"]["cm"], dtype=float)
    cm_s3_norm = cm_s3 / cm_s3.sum(axis=1, keepdims=True)
    im = ax.imshow(cm_s3_norm, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(6)); ax.set_yticks(range(6))
    ax.set_xticklabels(s3_classes); ax.set_yticklabels(s3_classes)
    for i in range(6):
        for j in range(6):
            ax.text(j, i, f"{int(cm_s3[i,j])}", ha="center", va="center",
                    color="white" if cm_s3_norm[i,j] > 0.5 else "black", fontsize=8)
    ax.set_xlabel("Predicted Permutation"); ax.set_ylabel("True Permutation")
    ax.set_title(f"(a) 6-Class $S_3$ Permutations (Acc={results['S3']['acc']*100:.1f}%, 0.3x={results['S3_03x']['acc']*100:.1f}%)")

    # (b) BIT3 Benchmarks vs 0.3x Controls
    ax = axes[0, 1]
    b_names = ["S3 (6-cls)", "BIT3 (8-way)", "Matched-1hot", "Matched-2hot", "3-Bit Parity"]
    b_keys = ["S3", "BIT3_8way", "BIT3_1hot3", "BIT3_2hot3", "BIT3_parity"]
    b_ctrl = ["S3_03x", "BIT3_03x_8way", "BIT3_03x_1hot3", "BIT3_03x_2hot3", "BIT3_03x_parity"]
    chances = [16.7, 12.5, 33.3, 33.3, 50.0]
    x = np.arange(len(b_names))
    ax.bar(x - 0.18, [results[k]["acc"]*100 for k in b_keys], width=0.35, color="#2ca02c", edgecolor="k", label="350 mV Reservoir")
    ax.bar(x + 0.18, [results[k]["acc"]*100 for k in b_ctrl], width=0.35, color="#7f7f7f", edgecolor="k", label="0.3x Control")
    for xi, ch in zip(x, chances):
        ax.plot([xi - 0.38, xi + 0.38], [ch, ch], "r--", lw=1.2)
    ax.set_xticks(x); ax.set_xticklabels(b_names, rotation=15, ha="right", fontsize=8.5)
    ax.set_ylim(0, 108)
    ax.set_ylabel("5-Fold CV Accuracy (%)")
    ax.set_title("(b) Multi-Symbol & Matched-Count Order vs. 0.3x Controls")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)

    # (c) BIT3 8-way Confusion Matrix
    ax = axes[1, 0]
    cm_b = np.array(results["BIT3_8way"]["cm"], dtype=float)
    cm_b_norm = cm_b / cm_b.sum(axis=1, keepdims=True)
    ax.imshow(cm_b_norm, cmap="Greens", vmin=0, vmax=1)
    ax.set_xticks(range(8)); ax.set_yticks(range(8))
    ax.set_xticklabels(bit3_words, fontsize=8); ax.set_yticklabels(bit3_words, fontsize=8)
    for i in range(8):
        for j in range(8):
            ax.text(j, i, f"{int(cm_b[i,j])}", ha="center", va="center",
                    color="white" if cm_b_norm[i,j] > 0.5 else "black", fontsize=7.5)
    ax.set_xlabel("Predicted 3-Bit Word"); ax.set_ylabel("True 3-Bit Word")
    ax.set_title(f"(c) 8-Way 3-Bit Word Classification (Acc={results['BIT3_8way']['acc']*100:.1f}%)")

    # (d) Active Learning Fisher J Trajectory
    ax = axes[1, 1]
    js = [t["j"] for t in traj]
    n_init = sum(1 for t in traj if t["phase"] == "INIT")
    ax.plot(range(1, n_init + 1), js[:n_init], "o-", color="#1f77b4", label=f"Sobol Init (mean J={results['active']['init_mean_j']:.2f})")
    ax.plot(range(n_init + 1, len(js) + 1), js[n_init:], "s-", color="#d62728", label=f"BO UCB (mean J={results['active']['bo_mean_j']:.2f})")
    ax.plot(range(1, len(js) + 1), np.maximum.accumulate(js), "k--", lw=1.2, label=f"Best J={results['active']['best_j']:.2f}")
    ax.set_xlabel("Evaluation Step")
    ax.set_ylabel("Order Fisher Separability $J$")
    ax.set_title(f"(d) Bayesian Active Learning (Pooled AB/BA={act_cv['acc']*100:.1f}%)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(HERE / "fig_multisymbol_and_active.png", bbox_inches="tight")
    fig.savefig(HERE / "fig_multisymbol_and_active.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {HERE / 'fig_multisymbol_and_active.png'}")


if __name__ == "__main__":
    main()
