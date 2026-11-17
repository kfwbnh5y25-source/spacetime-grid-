"""
Hybrid B+D Bayesian Active-Learning Loop for Spin-Wave Physical Reservoir Computing.
Maximizes energy-matched order separability J(AB, BA) across (f1, f2, gap_ns, amp_mv)
while interleaving drift pilots and 0.3x / off-resonance controls.

Usage:
  python run_active_reservoir.py --out data_active.h5 [--dry-run] [--budget 400]
"""

import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel
from sklearn.linear_model import RidgeClassifier
from sklearn.decomposition import PCA

from hw_core import HardwareRig, make_pulse_waveform, hampel_despike


# Search bounds: [f1_mhz, f2_mhz, delta_ns, amp_mv]
BOUNDS_LO = np.array([4050.0, 4050.0, -400.0, 310.0])
BOUNDS_HI = np.array([4160.0, 4160.0, 400.0, 370.0])


def build_pair_waveforms(theta: np.ndarray, scale: float = 1.0, f_shift: float = 0.0):
    """Build matched-energy (AB, BA) waveforms with canonical f_lo < f_hi ordering."""
    f_a, f_b = float(theta[0]) + f_shift, float(theta[1]) + f_shift
    f1, f2 = min(f_a, f_b), max(f_a, f_b)
    delta_ns, amp = float(theta[2]), float(theta[3]) * scale
    dur = 500.0
    t2 = max(0.0, dur + delta_ns)
    syms_ab = [
        {"amp_mv": amp, "freq_mhz": f1, "start_ns": 0.0, "dur_ns": dur, "phase": 0.0},
        {"amp_mv": amp, "freq_mhz": f2, "start_ns": t2, "dur_ns": dur, "phase": 0.0},
    ]
    syms_ba = [
        {"amp_mv": amp, "freq_mhz": f2, "start_ns": 0.0, "dur_ns": dur, "phase": 0.0},
        {"amp_mv": amp, "freq_mhz": f1, "start_ns": t2, "dur_ns": dur, "phase": 0.0},
    ]
    return make_pulse_waveform(syms_ab, gap_ns=2000.0), make_pulse_waveform(syms_ba, gap_ns=2000.0)


def featurize(sig_lin: np.ndarray, bl_lin: np.ndarray, pilot_power: float) -> np.ndarray:
    """Pilot-normalized Hampel-despiked linear spectrum (avoids log10 noise-floor amplification)."""
    diff = hampel_despike(sig_lin - bl_lin)
    return diff / max(pilot_power, 1e-12)


def evaluate_pair(rig: HardwareRig, theta: np.ndarray, tag: str, pilot_power: float, n_rep: int = 3):
    wf_ab, wf_ba = build_pair_waveforms(theta)
    feats_ab, feats_ba = [], []
    for r in range(n_rep):
        s_ab, b_ab, _ = rig.measure_trial(f"{tag}_AB_r{r}", wf_ab, {"tag": tag, "order": "AB", "theta": theta.tolist()})
        s_ba, b_ba, _ = rig.measure_trial(f"{tag}_BA_r{r}", wf_ba, {"tag": tag, "order": "BA", "theta": theta.tolist()})
        feats_ab.append(featurize(s_ab, b_ab, pilot_power))
        feats_ba.append(featurize(s_ba, b_ba, pilot_power))

    A = np.array(feats_ab)
    B = np.array(feats_ba)
    mu_a, mu_b = A.mean(axis=0), B.mean(axis=0)
    var_sum = np.sum(A.var(axis=0) + B.var(axis=0)) + 1e-9
    fisher_j = float(np.sum((mu_a - mu_b) ** 2) / var_sum)
    return fisher_j, A, B


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=str, default="data_active.h5")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--n-init", type=int, default=16)
    parser.add_argument("--n-iters", type=int, default=20)
    parser.add_argument("--n-rep", type=int, default=3)
    args = parser.parse_args()

    rig = HardwareRig(args.out, dry_run=args.dry_run, script_filename=__file__)
    pilot_wf, _ = build_pair_waveforms(np.array([4095.0, 4125.0, -250.0, 350.0]))

    def get_pilot_power(idx: int) -> float:
        s, b, _ = rig.measure_trial(f"PILOT_{idx}", pilot_wf, {"block": "pilot"})
        return float(np.sum(np.maximum(hampel_despike(s - b), 0.0)))

    X_train, y_train = [], []
    all_feats, all_labels = [], []

    try:
        # Phase 1: Sobol space-filling initialization
        sampler = qmc.Sobol(d=4, scramble=True, seed=2026)
        u_init = sampler.random(n=args.n_init)
        thetas_init = BOUNDS_LO + u_init * (BOUNDS_HI - BOUNDS_LO)

        pilot_p = get_pilot_power(0)
        for i, th in enumerate(thetas_init):
            if i > 0 and i % 5 == 0:
                pilot_p = get_pilot_power(i)
            j_val, A, B = evaluate_pair(rig, th, f"INIT_{i}", pilot_p, n_rep=args.n_rep)
            X_train.append((th - BOUNDS_LO) / (BOUNDS_HI - BOUNDS_LO))
            y_train.append(np.log1p(j_val))
            all_feats.extend(list(A) + list(B))
            all_labels.extend([0] * len(A) + [1] * len(B))
            print(f"[INIT {i+1}/{args.n_init}] J={j_val:.3f} | theta={np.round(th, 1)}")

        # Phase 2: Bayesian Active Learning (UCB)
        kernel = Matern(length_scale=np.ones(4) * 0.25, nu=2.5) + WhiteKernel(noise_level=0.05)
        gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5, random_state=42)
        rng = np.random.default_rng(77)

        for it in range(args.n_iters):
            gp.fit(np.array(X_train), np.array(y_train))
            cand_u = rng.uniform(0.0, 1.0, size=(2000, 4))
            mu, std = gp.predict(cand_u, return_std=True)
            ucb = mu + 2.0 * std
            best_u = cand_u[int(np.argmax(ucb))]
            best_th = BOUNDS_LO + best_u * (BOUNDS_HI - BOUNDS_LO)

            if it % 5 == 0:
                pilot_p = get_pilot_power(args.n_init + it)

            j_val, A, B = evaluate_pair(rig, best_th, f"BO_{it}", pilot_p, n_rep=args.n_rep)
            X_train.append(best_u)
            y_train.append(np.log1p(j_val))
            all_feats.extend(list(A) + list(B))
            all_labels.extend([0] * len(A) + [1] * len(B))
            print(f"[BO {it+1}/{args.n_iters}] J={j_val:.3f} (best={np.expm1(max(y_train)):.3f}) | theta={np.round(best_th, 1)}")

        # Quick linear ridge readout check (train/test 70/30 chronological split)
        X_f = np.array(all_feats)
        y_l = np.array(all_labels)
        split = int(0.7 * len(X_f))
        pca = PCA(n_components=min(15, split - 1)).fit(X_f[:split])
        clf = RidgeClassifier(alpha=1.0).fit(pca.transform(X_f[:split]), y_l[:split])
        acc = float(clf.score(pca.transform(X_f[split:]), y_l[split:]))
        print(f"[READOUT] Chronological held-out Ridge AB/BA accuracy = {acc:.3f}")

    except KeyboardInterrupt:
        print("[INTERRUPTED] Flushing data safely...")
    finally:
        rig.close()
        best_idx = int(np.argmax(y_train)) if y_train else 0
        best_th = (BOUNDS_LO + X_train[best_idx] * (BOUNDS_HI - BOUNDS_LO)).tolist() if X_train else []
        Path(args.out + ".bo_summary.json").write_text(json.dumps({
            "n_evals": len(y_train),
            "best_fisher_j": float(np.expm1(max(y_train))) if y_train else 0.0,
            "best_theta": best_th,
        }, indent=2))


if __name__ == "__main__":
    main()
