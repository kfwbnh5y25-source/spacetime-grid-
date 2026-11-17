"""
Analysis of Follow-Up Dataset (spin-waves/neuromorphic/data_followup.h5):
  - E1FIX: Integer-cycle G=0 cyclic collapse check
  - E3FIX: True inter-symbol delay [A][d][B] vs [B][d][A] at fixed T_L = 3.5 us (350 & 380 mV)
  - E5: Night-2 held-out generalization of frozen Night-1 (data_e0_e4.h5 E2_ovl) classifier
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
from sklearn.metrics import roc_auc_score
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


def load_h5(path: Path, prefix_filter=None):
    keys, specs = [], []
    with h5py.File(path, "r") as f:
        for sk in sorted(f.keys()):
            if not sk.endswith("__SIG"):
                continue
            name = sk[:-5]
            if prefix_filter and not name.startswith(prefix_filter):
                continue
            s = f[sk]["trace"][:, 1, :].mean(0).astype(np.float64)
            bk = name + "__BL"
            b = f[bk]["trace"][:, 1, :].mean(0).astype(np.float64) if bk in f else np.zeros_like(s)
            keys.append(name)
            specs.append(s - b)
    clean = hampel_batch(np.array(specs)) if specs else np.zeros((0, 1001))
    return dict(zip(keys, clean))


def eval_cv(X: np.ndarray, y: np.ndarray, seed: int = 42):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    preds, scores = np.zeros(len(y), dtype=int), np.zeros(len(y), dtype=float)
    for tr_i, te_i in skf.split(X, y):
        nc = min(10, len(tr_i) - 2, X.shape[1])
        sc = StandardScaler().fit(X[tr_i])
        pca = PCA(n_components=nc, random_state=seed).fit(sc.transform(X[tr_i]))
        clf = RidgeClassifier(alpha=10.0).fit(pca.transform(sc.transform(X[tr_i])), y[tr_i])
        preds[te_i] = clf.predict(pca.transform(sc.transform(X[te_i])))
        scores[te_i] = clf.decision_function(pca.transform(sc.transform(X[te_i])))
    A, B = X[y == 0], X[y == 1]
    j = float(np.sum((A.mean(0) - B.mean(0)) ** 2) / (np.sum(A.var(0) + B.var(0)) + 1e-12))
    return float(np.mean(preds == y)), float(roc_auc_score(y, scores)), j


def main():
    d_fu = load_h5(HERE / "data_followup.h5")
    print(f"Loaded {len(d_fu)} trials from data_followup.h5")

    # 1. E1FIX
    e1_res = {}
    for g in [0, 500, 2000]:
        by_o = {o: [v for k, v in d_fu.items() if k.startswith(f"E1FIX_{o}_G{g}_")] for o in ["ABC", "BCA", "CAB"]}
        if all(by_o.values()):
            mus = [np.mean(v, axis=0) for v in by_o.values()]
            d_bet = float(np.mean([np.linalg.norm(mus[0]-mus[1]), np.linalg.norm(mus[1]-mus[2]), np.linalg.norm(mus[2]-mus[0])]))
            w_var = np.mean([np.sum(np.array(v).var(0)) for v in by_o.values()])
            d_nf = float(np.sqrt(2.0 * w_var / len(by_o["ABC"])))
            X = np.vstack([by_o["ABC"], by_o["BCA"], by_o["CAB"]])
            y = np.array([0]*len(by_o["ABC"]) + [1]*len(by_o["BCA"]) + [2]*len(by_o["CAB"]))
            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
            pr = np.zeros(len(y), dtype=int)
            for tr, te in skf.split(X, y):
                sc = StandardScaler().fit(X[tr])
                pca = PCA(n_components=10, random_state=42).fit(sc.transform(X[tr]))
                clf = RidgeClassifier(alpha=10.0).fit(pca.transform(sc.transform(X[tr])), y[tr])
                pr[te] = clf.predict(pca.transform(sc.transform(X[te])))
            e1_res[g] = {"snr_ratio": d_bet / max(d_nf, 1e-12), "acc_3class": float(np.mean(pr == y))}
    print("[E1FIX]", json.dumps(e1_res, indent=2))

    # 2. E3FIX
    e3_res = {}
    for amp in [350, 380]:
        e3_res[amp] = {}
        for d in [0, 10, 20, 40, 80, 150, 300]:
            ab = [v for k, v in d_fu.items() if k.startswith(f"E3FIX_a{amp}_AB_d{d}_")]
            ba = [v for k, v in d_fu.items() if k.startswith(f"E3FIX_a{amp}_BA_d{d}_")]
            if ab and ba:
                X = np.vstack([ab, ba])
                y = np.array([0]*len(ab) + [1]*len(ba))
                acc, auc, j = eval_cv(X, y)
                e3_res[amp][d] = {"n": len(y), "acc": acc, "auc": auc, "fisher_j": j}
    print("[E3FIX]", json.dumps(e3_res, indent=2))

    # 3. E5 Cross-session generalization (Train on Night 1 E2_ovl, test on Night 2 E5_ovl)
    d_n1_ovl = load_h5(HERE / "data_e0_e4.h5", prefix_filter="E2_ovl_")
    d_n1_pil = load_h5(HERE / "data_e0_e4.h5", prefix_filter="E0_PILOT_")
    pil1_med = float(np.median([np.sum(np.maximum(v, 0)) for v in d_n1_pil.values()]))
    pil2_vals = [np.sum(np.maximum(v, 0)) for k, v in d_fu.items() if k.startswith("E0_PILOT_")]
    pil2_med = float(np.median(pil2_vals)) if pil2_vals else pil1_med
    scale_2to1 = pil1_med / max(pil2_med, 1e-12)

    X1 = np.array(list(d_n1_ovl.values()))
    y1 = np.array([0 if "_AB_" in k else 1 for k in d_n1_ovl.keys()])
    e5_items = [(k, v * scale_2to1) for k, v in d_fu.items() if k.startswith("E5_ovl_")]
    X2 = np.array([v for _, v in e5_items])
    y2 = np.array([0 if "_AB_" in k else 1 for k, _ in e5_items])

    sc = StandardScaler().fit(X1)
    pca = PCA(n_components=15, random_state=42).fit(sc.transform(X1))
    clf = RidgeClassifier(alpha=10.0).fit(pca.transform(sc.transform(X1)), y1)
    acc_e5 = float(clf.score(pca.transform(sc.transform(X2)), y2))
    auc_e5 = float(roc_auc_score(y2, clf.decision_function(pca.transform(sc.transform(X2)))))
    print(f"[E5 Cross-Session] pilot_scale={scale_2to1:.3f} | Night1->Night2 Held-Out Acc={acc_e5*100:.1f}%, AUC={auc_e5:.3f} (n={len(y2)})")

    # Save figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), dpi=170)
    fig.suptitle("Follow-Up Verification: Cyclic Control, True Memory Curve & Night-2 Generalization", fontsize=12, fontweight="bold")

    ax = axes[0]
    gs = sorted(e1_res.keys())
    ax.plot(gs, [e1_res[g]["snr_ratio"] for g in gs], "s-", color="#d62728")
    ax.axhline(1.0, color="k", ls="--", lw=0.8)
    ax.set_title("(a) E1FIX: Integer-Cycle Cyclic Collapse")
    ax.set_xlabel("Gap G (ns)")
    ax.set_ylabel("Cyclic Dist / Noise Floor")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ds = [0, 10, 20, 40, 80, 150, 300]
    for amp, col in [(350, "#1f77b4"), (380, "#2ca02c")]:
        ax.plot(ds, [e3_res[amp][d]["acc"]*100 for d in ds], "o-", color=col, label=f"{amp} mV")
    ax.axhline(50.0, color="k", ls="--", lw=0.8)
    ax.set_ylim(30, 105)
    ax.set_title("(b) E3FIX: True Inter-Symbol Delay [A][d][B]")
    ax.set_xlabel("Inter-Symbol Delay d (ns)")
    ax.set_ylabel("5-Fold CV Accuracy (%)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.bar(["Night 1 CV\n(n=320)", "Night 2 Frozen\n(n=80)"], [100.0, acc_e5 * 100], color=["#1f77b4", "#2ca02c"], edgecolor="k")
    ax.axhline(50.0, color="k", ls="--", lw=0.8)
    ax.axhline(75.0, color="r", ls=":", lw=1.0, label="Gate E5 (75%)")
    ax.set_ylim(0, 110)
    ax.set_title(f"(c) E5: Frozen Readout Generalization ({acc_e5*100:.1f}%)")
    ax.set_ylabel("Classification Accuracy (%)")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)

    plt.tight_layout()
    fig.savefig(HERE / "fig_followup_summary.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
