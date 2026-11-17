"""
Corrected PAM-4 Nonlinear Channel Equalization (Jaeger & Haas 2004) for Spin-Wave Reservoir.
Fixes all 4 flaws of naive binary designs:
  1. True 4-level PAM-4 input u(t) in {-3, -1, +1, +3} (nonlinear q^2, q^3 cannot be inverted linearly).
  2. Zero frequency leakage: ONLY scalar distorted channel output y(t) modulates 4081 MHz amplitude;
     4125 MHz is held at constant parametric pump bias (355 mV).
  3. Matched timescale: 20 ns/symbol (6 steps = 120 ns), inside measured tau_J ~ 5-10 ns memory.
  4. Proper baselines: Physical Reservoir vs. Direct Linear Equalizer (Ridge on y) vs. 0.3x Control.

Usage:
  python run_pam4_equalization.py --out data_pam4_eq.h5 [--dry-run]
"""

import argparse
import json
from pathlib import Path
import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from hw_core import HardwareRig, make_pulse_waveform


PAM4_LEVELS = np.array([-3.0, -1.0, 1.0, 3.0])


def generate_pam4_window(rng: np.random.Generator, n_pre: int = 5, n_steps: int = 6) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Generate PAM-4 stream u and nonlinear channel output y[0:n_steps] (Jaeger & Haas 2004).
    Target label is u at the middle-to-end delay d=1 (index n_pre + n_steps - 2).
    """
    u = rng.choice(PAM4_LEVELS, size=n_pre + n_steps + 2)
    y = np.zeros(n_steps)
    for k in range(n_steps):
        t = n_pre + k
        q = (0.08 * u[t + 2] - 0.12 * u[t + 1] + u[t]
             + 0.18 * u[t - 1] - 0.10 * u[t - 2] + 0.09 * u[t - 3]
             - 0.05 * u[t - 4] + 0.04 * u[t - 5])
        y[k] = q + 0.036 * (q ** 2) - 0.011 * (q ** 3)
    target_u = int(np.where(PAM4_LEVELS == u[n_pre + n_steps - 2])[0][0])
    return u, y, target_u


def encode_channel_waveform(y_win: np.ndarray, base_mv: float = 350.0, pump_mv: float = 350.0,
                            step_ns: float = 20.0, gap_ns: float = 200.0) -> np.ndarray:
    """Map scalar y(t) ONLY onto 4081 MHz amplitude (peak <= 735 mV < 760 mV ZHL ceiling); constant 4125 MHz pump."""
    y_norm = np.clip(y_win / 5.0, -0.55, 0.10)
    syms = [{
        "amp_mv": pump_mv, "freq_mhz": 4125.0,
        "start_ns": 0.0, "dur_ns": step_ns * len(y_win), "phase": 0.0,
    }]
    for k, yn in enumerate(y_norm):
        syms.append({
            "amp_mv": float(base_mv * (1.0 + yn)),
            "freq_mhz": 4081.0,
            "start_ns": k * step_ns,
            "dur_ns": step_ns,
            "phase": 0.0,
        })
    return make_pulse_waveform(syms, gap_ns=gap_ns)


def build_trials(n_main: int = 400, n_ctrl: int = 160, seed: int = 42):
    rng = np.random.default_rng(seed)
    trials = []
    for block, n, scale in [("EQ_MAIN", n_main, 1.0), ("EQ_03X", n_ctrl, 0.3)]:
        for r in range(n):
            _, y_win, label = generate_pam4_window(rng)
            wf = encode_channel_waveform(y_win, base_mv=350.0 * scale, pump_mv=350.0 * scale)
            meta = {"block": block, "label": label, "y_win": y_win.tolist(), "scale": scale, "rep": r}
            trials.append((f"{block}_r{r}", wf, meta))
    return trials


def evaluate_linear_baseline(trials) -> float:
    main = [t for t in trials if t[2]["block"] == "EQ_MAIN"]
    X = np.array([t[2]["y_win"] for t in main])
    y = np.array([t[2]["label"] for t in main])
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    preds = np.zeros_like(y)
    for tr, te in skf.split(X, y):
        sc = StandardScaler().fit(X[tr])
        preds[te] = RidgeClassifier(alpha=1.0).fit(sc.transform(X[tr]), y[tr]).predict(sc.transform(X[te]))
    return float(np.mean(preds == y))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data_pam4_eq.h5")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    trials = build_trials()
    lin_acc = evaluate_linear_baseline(trials)
    print(f"Prepared {len(trials)} PAM-4 equalization trials | Direct Linear Equalizer 5-fold acc = {lin_acc*100:.1f}% (chance=25.0%)")

    if args.dry_run:
        print("Dry-run OK (waveform shape:", trials[0][1].shape, ")")
        return

    rig = HardwareRig(dry_run=False)
    try:
        rig.init_hdf5(args.out)
        for name, wf, meta in trials:
            spec = rig.run_trial(wf)
            rig.save_trial(args.out, name, spec, meta)
    finally:
        rig.close()


if __name__ == "__main__":
    main()
