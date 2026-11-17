"""
Rank-saturation curve: accuracy vs. PCA components for each task.
Saves fig_rank_saturation.{png,pdf} to spin-waves/reports/neuromorphic/figures/
"""
import h5py, numpy as np, re
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeClassifier
from sklearn.model_selection import StratifiedKFold
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "reports" / "neuromorphic" / "figures"
OUT.mkdir(exist_ok=True)


def get_freq_mask(f_h5, key, lo=2030, hi=2075):
    tr = f_h5[key]["trace"]
    freq = tr[0, 0, :] / 1e6
    return (freq >= lo) & (freq <= hi)


def cv_at_rank(X, y, rank, n_splits=5, seed=42):
    if len(y) == 0:
        return 0.0
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    preds = np.zeros(len(y), dtype=int)
    r = min(rank, len(y) // n_splits - 2, X.shape[1] - 1)
    if r < 1:
        return float(np.mean(np.zeros(len(y), dtype=int) == y))
    for tr_idx, te_idx in skf.split(X, y):
        sc = StandardScaler().fit(X[tr_idx])
        pca = PCA(n_components=r, random_state=seed).fit(sc.transform(X[tr_idx]))
        Z_tr = pca.transform(sc.transform(X[tr_idx]))
        Z_te = pca.transform(sc.transform(X[te_idx]))
        preds[te_idx] = RidgeClassifier(alpha=1.0).fit(Z_tr, y[tr_idx]).predict(Z_te)
    return float(np.mean(preds == y))


# ── Load datasets ──────────────────────────────────────────────────────────────

datasets = {}

# 1. E2: AB vs BA
X, y = [], []
with h5py.File(HERE / "data_e0_e4.h5", "r") as f:
    keys = [k for k in f.keys() if "__SIG" in k and "E2" in k and "ovl" in k]
    mask = None
    for k in keys:
        bl_k = k.replace("__SIG", "__BL")
        if bl_k not in f:
            continue
        if mask is None:
            mask = get_freq_mask(f, k)
        s = f[k]["trace"][:, 1, :].mean(0) - f[bl_k]["trace"][:, 1, :].mean(0)
        X.append(s[mask])
        y.append(1 if "BA" in k else 0)
datasets["E2: AB vs BA"] = (np.array(X), np.array(y), 0.5)
print(f"E2 loaded: {len(X)} trials")

# 2. S3: 6-class permutations  (key format: S3_ABC_p0_r0__SIG, excluding 03x)
s3_cls = ["ABC", "ACB", "BAC", "BCA", "CAB", "CBA"]
X, y = [], []
with h5py.File(HERE / "data_multisymbol_v2.h5", "r") as f:
    mask = None
    for k in f.keys():
        m = re.match(r"^S3_([A-Z]{3})_p\d+_r\d+__SIG$", k)
        if not m or m.group(1) not in s3_cls:
            continue
        bl_k = k.replace("__SIG", "__BL")
        if bl_k not in f:
            continue
        if mask is None:
            mask = get_freq_mask(f, k)
        s = f[k]["trace"][:, 1, :].mean(0) - f[bl_k]["trace"][:, 1, :].mean(0)
        X.append(s[mask])
        y.append(s3_cls.index(m.group(1)))
datasets["S3: 6-class perms"] = (np.array(X), np.array(y), 1 / 6)
print(f"S3 loaded: {len(X)} trials")

# 3. BIT3 8-way  (key format: BIT3_000_p0_r0__SIG)
bit3_cls = ["000", "001", "010", "011", "100", "101", "110", "111"]
X_b, y_b = [], []
with h5py.File(HERE / "data_multisymbol_v2.h5", "r") as f:
    mask = None
    for k in f.keys():
        m = re.match(r"^BIT3_([01]{3})_p\d+_r\d+__SIG$", k)
        if not m or m.group(1) not in bit3_cls:
            continue
        bl_k = k.replace("__SIG", "__BL")
        if bl_k not in f:
            continue
        if mask is None:
            mask = get_freq_mask(f, k)
        s = f[k]["trace"][:, 1, :].mean(0) - f[bl_k]["trace"][:, 1, :].mean(0)
        X_b.append(s[mask])
        y_b.append(bit3_cls.index(m.group(1)))
X_b, y_b = np.array(X_b), np.array(y_b)
datasets["BIT3: 8-way words"] = (X_b, y_b, 0.125)
print(f"BIT3 loaded: {len(y_b)} trials")

# 4. Matched-count subsets
idx_2hot = [i for i, lbl in enumerate(y_b) if bin(lbl).count("1") == 2]
idx_1hot = [i for i, lbl in enumerate(y_b) if bin(lbl).count("1") == 1]
lab_2hot = [["011", "101", "110"].index(bit3_cls[y_b[i]]) for i in idx_2hot]
lab_1hot = [["001", "010", "100"].index(bit3_cls[y_b[i]]) for i in idx_1hot]
datasets["BIT3-2hot {011,101,110}"] = (X_b[idx_2hot], np.array(lab_2hot), 1 / 3)
datasets["BIT3-1hot {001,010,100}"] = (X_b[idx_1hot], np.array(lab_1hot), 1 / 3)
print(f"2hot: {len(idx_2hot)}, 1hot: {len(idx_1hot)}")

# ── Sweep ranks ────────────────────────────────────────────────────────────────

ranks = list(range(1, 51))
colors = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd"]

fig, ax = plt.subplots(figsize=(11, 6.5), dpi=180)

sat_ranks = {}
for (label, (X, y, chance)), color in zip(datasets.items(), colors):
    accs = [cv_at_rank(X, y, r) for r in ranks]
    amax = max(accs)
    sat_idx = next((i for i, a in enumerate(accs) if a >= 0.95 * amax), len(accs) - 1)
    sat_ranks[label] = (ranks[sat_idx], amax)

    ax.plot(ranks, [a * 100 for a in accs], "-o", color=color,
            markersize=3, lw=1.8, label=f"{label}  [sat@{ranks[sat_idx]}, max={amax*100:.1f}%]")
    ax.axhline(chance * 100, color=color, lw=0.9, linestyle="--", alpha=0.45)
    ax.axvline(ranks[sat_idx], color=color, lw=0.8, linestyle=":", alpha=0.6)
    ax.text(ranks[sat_idx] + 0.4, chance * 100 + 1.5, str(ranks[sat_idx]),
            color=color, fontsize=8)

ax.set_xlabel("PCA Components (Reservoir Rank)", fontsize=13)
ax.set_ylabel("5-Fold CV Accuracy (%)", fontsize=13)
ax.set_title("Rank-Saturation Curve — Useful Dimensionality of the Physical Reservoir",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=8.5, loc="lower right", framealpha=0.9)
ax.grid(True, alpha=0.3)
ax.set_xlim(1, 50)
ax.set_ylim(0, 108)
ax.text(0.02, 0.97, "Dashed lines = chance level\nDotted verticals = 95%-saturation rank",
        transform=ax.transAxes, fontsize=8, va="top", alpha=0.6)

plt.tight_layout()
fig.savefig(OUT / "fig_rank_saturation.png", bbox_inches="tight")
fig.savefig(OUT / "fig_rank_saturation.pdf", bbox_inches="tight")
plt.close()
print(f"\nSaved to {OUT / 'fig_rank_saturation.png'}")

print("\nSaturation summary:")
for label, (sat_r, amax) in sat_ranks.items():
    print(f"  {label:35s}: 95%-sat rank = {sat_r:2d},  max acc = {amax*100:.1f}%")
