"""
Analysis of spin-waves/neuromorphic/data_multisymbol_v2.h5:
  - S3 (6-class permutations across varied phases) vs S3_03x
  - BIT3 (8-way, Matched-1hot, Matched-2hot, 3-bit parity) vs BIT3_03x
  - PHASE_GATE (AB vs BA across dphi in {0, pi/4, pi/2, 3pi/4, pi})

Outputs:
  - multisymbol_v2_metrics.json
  - fig_multisymbol_v2.png / .pdf
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
H5_V2 = HERE / "data_multisymbol_v2.h5"


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


def eval_cv(X: np.ndarray, y: np.ndarray, n_comp: int = 20, n_perm: int = 200, seed: int = 42):
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
        preds[te_i] = RidgeClassifier(alpha=1.0).fit(Z_tr, y[tr_i]).predict(Z_te)
    acc = float(np.mean(preds == y))
    cm = confusion_matrix(y, preds)

    rng = np.random.default_rng(seed)
    nulls = []
    for _ in range(n_perm):
        yp = rng.permutation(y)
        pp = np.zeros(len(y), dtype=int)
        for tr_i, te_i, Z_tr, Z_te in Z_folds:
            pp[te_i] = RidgeClassifier(alpha=1.0).fit(Z_tr, yp[tr_i]).predict(Z_te)
        nulls.append(np.mean(pp == yp))
    p_val = float((1 + np.sum(np.array(nulls) >= acc)) / (n_perm + 1))
    return {"n": int(len(y)), "acc": acc, "p_val": p_val, "cm": cm.tolist()}


def main():
    keys, specs = [], []
    with h5py.File(H5_V2, "r") as f:
        fr = None
        for sk in sorted(f.keys()):
            if not sk.endswith("__SIG"):
                continue
            name = sk[:-5]
            tr = f[sk]["trace"][:]
            if fr is None:
                fr = tr[0, 0, :] / 1e6
            s = tr[:, 1, :].mean(0).astype(np.float64)
            bk = name + "__BL"
            b = f[bk]["trace"][:, 1, :].mean(0).astype(np.float64) if bk in f else np.zeros_like(s)
            keys.append(name)
            specs.append(s - b)

    clean = hampel_batch(np.array(specs))
    band_mask = (fr >= 2030.0) & (fr <= 2075.0)
    data = {k: v[band_mask] for k, v in zip(keys, clean)}
    print(f"Loaded {len(data)} trials from {H5_V2.name}")

    results = {}
    s3_cls = ["ABC", "ACB", "BAC", "BCA", "CAB", "CBA"]
    for lbl in ["S3", "S3_03x"]:
        X, y = [], []
        for k, v in data.items():
            m = re.match(rf"^{lbl}_([ABC]{{3}})_p\d+_r\d+$", k)
            if m and m.group(1) in s3_cls:
                X.append(v); y.append(s3_cls.index(m.group(1)))
        results[lbl] = eval_cv(np.array(X), np.array(y), n_comp=20)
        print(f"[{lbl}] 6-class Acc = {results[lbl]['acc']*100:.1f}% (chance 16.7%, p={results[lbl]['p_val']:.4f})")

    bit3_words = ["000", "001", "010", "011", "100", "101", "110", "111"]
    w1_set, w2_set = ["001", "010", "100"], ["011", "101", "110"]
    for lbl in ["BIT3", "BIT3_03x"]:
        X_all, y_all, y_par, X_w1, y_w1, X_w2, y_w2 = [], [], [], [], [], [], []
        for k, v in data.items():
            m = re.match(rf"^{lbl}_([01]{{3}})_p\d+_r\d+$", k)
            if m:
                w = m.group(1)
                X_all.append(v); y_all.append(bit3_words.index(w)); y_par.append(sum(int(c) for c in w) % 2)
                if w in w1_set: X_w1.append(v); y_w1.append(w1_set.index(w))
                elif w in w2_set: X_w2.append(v); y_w2.append(w2_set.index(w))
        results[f"{lbl}_8way"] = eval_cv(np.array(X_all), np.array(y_all), n_comp=20)
        results[f"{lbl}_parity"] = eval_cv(np.array(X_all), np.array(y_par), n_comp=20)
        results[f"{lbl}_1hot3"] = eval_cv(np.array(X_w1), np.array(y_w1), n_comp=15)
        results[f"{lbl}_2hot3"] = eval_cv(np.array(X_w2), np.array(y_w2), n_comp=15)
        print(f"[{lbl}] 8-way={results[f'{lbl}_8way']['acc']*100:.1f}% | Parity={results[f'{lbl}_parity']['acc']*100:.1f}% | "
              f"Matched-1hot={results[f'{lbl}_1hot3']['acc']*100:.1f}% | Matched-2hot={results[f'{lbl}_2hot3']['acc']*100:.1f}%")

    # Task C: Phase Gate
    ph_res = {}
    for pi_idx, ph_deg in enumerate([0, 45, 90, 135, 180]):
        X, y = [], []
        for k, v in data.items():
            m = re.match(rf"^PH_(AB|BA)_p{pi_idx}_r\d+$", k)
            if m:
                X.append(v); y.append(0 if m.group(1) == "AB" else 1)
        if X:
            X_a, y_a = np.array(X), np.array(y)
            A, B = X_a[y_a == 0], X_a[y_a == 1]
            j = float(np.sum((A.mean(0) - B.mean(0))**2) / (np.sum(A.var(0) + B.var(0)) + 1e-12))
            r_cv = eval_cv(X_a, y_a, n_comp=min(5, len(X_a)//2))
            ph_res[ph_deg] = {"acc": r_cv["acc"], "fisher_j": j, "ab_sum": float(A.sum(1).mean()), "ba_sum": float(B.sum(1).mean())}
    results["phase_gate"] = ph_res
    print("[PHASE_GATE]", json.dumps(ph_res, indent=2))

    (HERE / "multisymbol_v2_metrics.json").write_text(json.dumps(results, indent=2))

    # --- 3-Panel Summary Figure ---
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6), dpi=180)
    fig.suptitle("Multi-Symbol v2 (Phase-Varied, Equalized Drive): S3, 3-Bit Words & Phase Gating", fontsize=12.5, fontweight="bold")

    ax = axes[0]
    cm = np.array(results["S3"]["cm"], dtype=float)
    cm_n = cm / cm.sum(axis=1, keepdims=True)
    ax.imshow(cm_n, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(6)); ax.set_yticks(range(6))
    ax.set_xticklabels(s3_cls, fontsize=8); ax.set_yticklabels(s3_cls, fontsize=8)
    for i in range(6):
        for j in range(6):
            ax.text(j, i, f"{int(cm[i,j])}", ha="center", va="center",
                    color="white" if cm_n[i,j] > 0.5 else "black", fontsize=8)
    ax.set_xlabel("Predicted"); ax.set_ylabel("True")
    ax.set_title(f"(a) 6-Class $S_3$ (Acc={results['S3']['acc']*100:.1f}%, 0.3x={results['S3_03x']['acc']*100:.1f}%)")

    ax = axes[1]
    b_names = ["S3 (6-cls)", "BIT3 (8-way)", "Matched-1hot", "Matched-2hot", "3-Bit Parity"]
    b_keys = ["S3", "BIT3_8way", "BIT3_1hot3", "BIT3_2hot3", "BIT3_parity"]
    b_ctrl = ["S3_03x", "BIT3_03x_8way", "BIT3_03x_1hot3", "BIT3_03x_2hot3", "BIT3_03x_parity"]
    chances = [16.7, 12.5, 33.3, 33.3, 50.0]
    x = np.arange(len(b_names))
    ax.bar(x - 0.18, [results[k]["acc"]*100 for k in b_keys], width=0.35, color="#2ca02c", edgecolor="k", label="v2 Reservoir")
    ax.bar(x + 0.18, [results[k]["acc"]*100 for k in b_ctrl], width=0.35, color="#7f7f7f", edgecolor="k", label="0.3x Ctrl")
    for xi, ch in zip(x, chances):
        ax.plot([xi - 0.38, xi + 0.38], [ch, ch], "r--", lw=1.2)
    ax.set_xticks(x); ax.set_xticklabels(b_names, rotation=15, ha="right", fontsize=8)
    ax.set_ylim(0, 108); ax.set_ylabel("5-Fold CV Accuracy (%)")
    ax.set_title("(b) Benchmarks vs. 0.3x Controls")
    ax.legend(fontsize=8); ax.grid(True, axis="y", alpha=0.3)

    ax = axes[2]
    degs = sorted(ph_res.keys())
    ax.plot(degs, [ph_res[d]["fisher_j"] for d in degs], "o-", color="#d62728", lw=1.6, label="Fisher J(AB, BA)")
    ax2 = ax.twinx()
    ax2.plot(degs, [ph_res[d]["acc"]*100 for d in degs], "s--", color="#1f77b4", lw=1.3, label="AB/BA Acc (%)")
    ax2.set_ylim(40, 110)
    ax.set_xlabel("Relative Phase $\\Delta\\phi$ (deg)")
    ax.set_ylabel("Fisher Separability $J$", color="#d62728")
    ax2.set_ylabel("AB/BA Accuracy (%)", color="#1f77b4")
    ax.set_title("(c) Explicit 4B Phase-Gating Sweep")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(HERE / "fig_multisymbol_v2.png", bbox_inches="tight")
    fig.savefig(HERE / "fig_multisymbol_v2.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {HERE / 'fig_multisymbol_v2.png'}")


if __name__ == "__main__":
    main()
