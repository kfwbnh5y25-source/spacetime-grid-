"""
Unified High-Accuracy Neuromorphic Benchmark Runner:
  Task 1: Tuned Nonlinear Channel Equalization (EQ-Plus)
          - Inter-symbol overlap tuned to 80 ns (matches 4-magnon physical memory lifetime tau_mem)
          - Inverts nonlinear ISI channel: q(t) = s(t) + 0.18*s(t-1) - 0.22*s(t-2), y(t) = q(t) - 0.15*q(t)^3
          - Target: High-accuracy recovery of delayed symbol s(t-1)
  Task 2: Single-Speaker Isolated Spoken Digit Recognition (FSDD Jackson 4-Class)
          - Gold-standard benchmark speaker (Jackson) following Torrejon et al. (Nature 2017)
          - Eliminates cross-speaker acoustic variance
          - Sequential single-tone phoneme drive with 100 ns overlap, V_peak <= 375 mV (no amp compression)

Features:
  - Both tasks unified into ONE script and ONE HDF5 file (data_high_accuracy.h5)
  - Real-time online accuracy calculation after every trial with color-coded status flags:
      [SUCCESS] (Acc >= 85%) | [PARTIAL] (70% <= Acc < 85%) | [FAIL] (Acc < 70%) | [WARMUP]
  - Balanced 0.3x sub-threshold controls for both tasks
  - Zero dependencies beyond numpy/scipy/h5py (pure numpy Ridge CV fallback)
  - Offline --dry-run support

Usage:
  python run_high_accuracy_benchmarks.py --dry-run
  python run_high_accuracy_benchmarks.py --task all --out data_high_accuracy.h5
  python run_high_accuracy_benchmarks.py --task eq
  python run_high_accuracy_benchmarks.py --task digits
"""

import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np

from hw_core import HardwareRig, make_pulse_waveform

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = str(HERE / "data_high_accuracy.h5")
FEATURES_FILE = HERE / "spoken_digits_features_v2.json"

# Hardware RF Calibration Parameters
F_A_MHZ = 4081.0  # Dominant low formant / Symbol +1 carrier
F_B_MHZ = 4125.0  # Dominant high formant / Symbol -1 carrier
BASE_AMP_MV = 365.0
AMP_B_MV = 350.0

# Band for subharmonic readout integration
SUBH_LOW_HZ = 2030e6
SUBH_HIGH_HZ = 2075e6


# ==============================================================================
# Pure NumPy Online Cross-Validation Engine (Zero sklearn dependency required)
# ==============================================================================

def pure_numpy_ridge_cv(X_list: list[np.ndarray], y_list: list[int], n_splits: int = 5, alpha: float = 1.0, seed: int = 42) -> float:
    """Computes Stratified K-Fold CV accuracy with L2 Ridge regression in pure NumPy."""
    y = np.array(y_list)
    classes = np.unique(y)
    K = len(classes)
    N = len(y)

    if K < 2:
        return float(np.nan)
    min_count = min(np.sum(y == c) for c in classes)
    if min_count < 2:
        return float(np.nan)

    actual_splits = min(n_splits, min_count)
    if actual_splits < 2:
        return float(np.nan)

    X = np.array(X_list)
    rng = np.random.default_rng(seed)

    # Stratified fold generation
    folds = [[] for _ in range(actual_splits)]
    for c in classes:
        c_indices = rng.permutation(np.where(y == c)[0])
        for idx, item in enumerate(c_indices):
            folds[idx % actual_splits].append(item)

    # 1-hot target matrix
    Y_onehot = np.zeros((N, K))
    for i, c in enumerate(y):
        Y_onehot[i, np.where(classes == c)[0][0]] = 1.0

    preds = np.zeros(N, dtype=int)
    for k in range(actual_splits):
        te_idx = np.array(folds[k])
        tr_idx = np.array([i for i in range(N) if i not in te_idx])
        if len(tr_idx) == 0 or len(te_idx) == 0:
            continue

        X_tr, Y_tr = X[tr_idx], Y_onehot[tr_idx]
        X_te = X[te_idx]

        # Standardize using training distribution
        mu = X_tr.mean(axis=0)
        std = X_tr.std(axis=0) + 1e-8
        Z_tr = (X_tr - mu) / std
        Z_te = (X_te - mu) / std

        # Normal equations: W = (Z_tr^T Z_tr + alpha * I)^(-1) Z_tr^T Y_tr
        d = Z_tr.shape[1]
        try:
            W = np.linalg.solve(Z_tr.T @ Z_tr + alpha * np.eye(d), Z_tr.T @ Y_tr)
        except np.linalg.LinAlgError:
            W = np.linalg.pinv(Z_tr.T @ Z_tr + alpha * np.eye(d)) @ (Z_tr.T @ Y_tr)

        logits = Z_te @ W
        preds[te_idx] = classes[np.argmax(logits, axis=1)]

    return float(np.mean(preds == y))


# ==============================================================================
# Task 1: High-Accuracy Tuned Channel Equalization (EQ-Plus)
# ==============================================================================

def generate_channel_eq_trials(n_trials: int, seed: int = 42) -> list[dict]:
    """Generates random binary symbol sequences passing through a nonlinear ISI channel."""
    rng = np.random.default_rng(seed)
    trials = []
    for r in range(n_trials):
        s = rng.choice([-1, 1], size=3)
        # Channel ISI:
        q = np.zeros(3)
        q[0] = s[0]
        q[1] = s[1] + 0.18 * s[0]
        q[2] = s[2] + 0.18 * s[1] - 0.22 * s[0]
        # Cubic amplifier saturation:
        y = q - 0.15 * (q ** 3)
        y_norm = np.clip(y / 1.5, -1.0, 1.0)

        trials.append({
            "trial_id": r,
            "symbols": [int(x) for x in s],
            "channel_out": [float(x) for x in y_norm],
            "target_lag0": int(s[2] > 0),
            "target_lag1": int(s[1] > 0),  # Primary memory target
            "target_lag2": int(s[0] > 0),
        })
    return trials


def build_eq_waveform_tuned(trial: dict, scale: float = 1.0, dur_ns: float = 300.0, step_ns: float = 220.0, gap_ns: float = 1500.0) -> np.ndarray:
    """
    Constructs the RF pulse sequence with an 80 ns overlap between consecutive symbols:
    dur_ns=300, step_ns=220 -> 80 ns overlap perfectly couples consecutive symbols
    via the measured 4-magnon interaction lifetime (tau_mem = 12-60 ns).
    """
    syms = []
    t_cur = 0.0
    for i in range(3):
        s_val = trial["symbols"][i]
        freq = F_A_MHZ if s_val > 0 else F_B_MHZ
        # Modulate amplitude around supra-threshold operating point
        amp = (BASE_AMP_MV + 25.0 * trial["channel_out"][i]) * scale
        syms.append({
            "amp_mv": float(amp),
            "freq_mhz": float(freq),
            "start_ns": float(t_cur),
            "dur_ns": float(dur_ns),
            "phase": 0.0,
        })
        t_cur += step_ns

    return make_pulse_waveform(syms, gap_ns=gap_ns)


# ==============================================================================
# Task 2: Single-Speaker Spoken Digit Recognition (Jackson 4-Class)
# ==============================================================================

def load_jackson_dataset() -> list[dict]:
    """Loads acoustic formant features for single speaker Jackson (digits 0, 1, 2, 3)."""
    if not FEATURES_FILE.exists():
        raise FileNotFoundError(f"Features file not found: {FEATURES_FILE}")
    all_feats = json.loads(FEATURES_FILE.read_text())
    jackson_feats = [f for f in all_feats if f["speaker"] == "jackson"]
    return jackson_feats


def build_digit_waveform_jackson(trial: dict, scale: float = 1.0, dphi: float = 0.0, dur_ns: float = 450.0, step_ns: float = 350.0, gap_ns: float = 1500.0) -> np.ndarray:
    """
    5-step phoneme transition waveform:
    Dominant formant selects carrier (Tone A for vowels, Tone B for consonants)
    100 ns overlap between adjacent phonemes (dur_ns=450, step_ns=350).
    Peak amplitude <= 375 mV -> Zero amplifier compression.
    """
    low_env = trial["low_formant"]
    high_env = trial["high_formant"]
    n_steps = len(low_env)

    syms = []
    t_cur = 0.0
    for i in range(n_steps):
        lo_val = low_env[i]
        hi_val = high_env[i]

        if lo_val >= hi_val:
            freq = F_A_MHZ
            base_amp = BASE_AMP_MV
            mod = 15.0 * (lo_val - 0.5)
        else:
            freq = F_B_MHZ
            base_amp = AMP_B_MV
            mod = 15.0 * (hi_val - 0.5)

        amp = (base_amp + mod) * scale
        syms.append({
            "amp_mv": float(amp),
            "freq_mhz": float(freq),
            "start_ns": float(t_cur),
            "dur_ns": float(dur_ns),
            "phase": float(i * dphi),
        })
        t_cur += step_ns

    return make_pulse_waveform(syms, gap_ns=gap_ns)


def build_pilot_wf() -> np.ndarray:
    """Standard 2-tone reference pilot for drift tracking."""
    syms = [
        {"amp_mv": 365.0, "freq_mhz": F_A_MHZ, "start_ns": 0.0, "dur_ns": 400.0, "phase": 0.0},
        {"amp_mv": 350.0, "freq_mhz": F_B_MHZ, "start_ns": 200.0, "dur_ns": 400.0, "phase": 0.0},
    ]
    return make_pulse_waveform(syms, gap_ns=1500.0)


# ==============================================================================
# Unified Execution Engine
# ==============================================================================

def run_benchmarks(out_h5: str, task: str = "all", dry_run: bool = False, n_eq: int = 80, n_eq_ctrl: int = 20, digit_reps: int = 3, n_digit_ctrl: int = 16, seed: int = 42):
    print("=" * 78)
    print("UNIFIED HIGH-ACCURACY PHYSICAL RESERVOIR BENCHMARK RUNNER")
    print(f"HDF5 Target:      {out_h5}")
    print(f"Target Tasks:     {task.upper()} | Dry-run: {dry_run}")
    print("=" * 78)

    rig = HardwareRig(out_h5, dry_run=dry_run, script_filename=__file__)
    pilot_wf = build_pilot_wf()

    # Prepare Task 1 trials
    eq_trials = generate_channel_eq_trials(n_eq, seed=seed)
    eq_ctrl_trials = generate_channel_eq_trials(n_eq_ctrl, seed=seed + 100)

    # Prepare Task 2 trials (Single Speaker Jackson, 20 base utterances * reps)
    jackson_base = load_jackson_dataset()
    phases = [0.0, 0.5 * np.pi, np.pi]
    digit_trials = []
    rng = np.random.default_rng(seed)
    for r in range(digit_reps):
        p_idx = r % len(phases)
        dphi = phases[p_idx]
        for item in jackson_base:
            digit_trials.append({
                "item": item,
                "dphi": dphi,
                "rep": r,
                "digit": item["digit"],
                "scale": 1.0,
            })
    rng.shuffle(digit_trials)

    # Sub-threshold controls for Digits (balanced per digit)
    digit_ctrl_trials = []
    per_digit_ctrl = n_digit_ctrl // 4
    by_digit = {}
    for item in jackson_base:
        by_digit.setdefault(item["digit"], []).append(item)
    for d in sorted(by_digit.keys()):
        for item in by_digit[d][:per_digit_ctrl]:
            digit_ctrl_trials.append({
                "item": item,
                "dphi": 0.0,
                "rep": 0,
                "digit": d,
                "scale": 0.3,
            })
    rng.shuffle(digit_ctrl_trials)

    # Calculate totals
    n_eq_total = len(eq_trials) + len(eq_ctrl_trials) if task in ["all", "eq"] else 0
    n_digit_total = len(digit_trials) + len(digit_ctrl_trials) if task in ["all", "digits"] else 0
    n_pilots = (n_eq_total // 25 + 2 if n_eq_total > 0 else 0) + (n_digit_total // 25 + 2 if n_digit_total > 0 else 0)
    total_acquisitions = n_eq_total + n_digit_total + n_pilots

    sec_per_trial = 0.02 if dry_run else 4.8
    est_total_min = (total_acquisitions * sec_per_trial) / 60.0

    print(f"Total Planned Acquisitions:     {total_acquisitions}")
    if task in ["all", "eq"]:
        print(f"  • Task 1 (Channel Equalization): {len(eq_trials)} res + {len(eq_ctrl_trials)} ctrl")
    if task in ["all", "digits"]:
        print(f"  • Task 2 (Jackson Spoken Digits): {len(digit_trials)} res + {len(digit_ctrl_trials)} ctrl")
    print(f"Estimated Total Runtime:        {est_total_min:.1f} minutes")
    print("-" * 78)

    completed = 0
    t0 = time.time()

    def print_progress(tag_name: str, desc: str, acc_val: float, target_acc: float = 0.85):
        nonlocal completed
        elapsed = time.time() - t0
        pct = (completed / max(1, total_acquisitions)) * 100
        rem_sec = (elapsed / max(1, completed)) * (total_acquisitions - completed) if completed > 0 else 0

        if np.isnan(acc_val):
            status = "\033[94m[WARMUP]\033[0m"
            acc_str = "Acc = ... "
        elif acc_val >= target_acc:
            status = "\033[92m[SUCCESS]\033[0m"
            acc_str = f"Acc = {acc_val*100:4.1f}%"
        elif acc_val >= 0.70:
            status = "\033[93m[PARTIAL]\033[0m"
            acc_str = f"Acc = {acc_val*100:4.1f}%"
        else:
            status = "\033[91m[FAIL]\033[0m"
            acc_str = f"Acc = {acc_val*100:4.1f}%"

        print(f"[{completed:3d}/{total_acquisitions:3d} | {pct:4.1f}%] {status} ({acc_str}) {tag_name:<28s} | {desc:<20s} | Rem: {rem_sec/60.0:4.1f} min", flush=True)

    pilot_idx = 0
    try:
        # Initial Reference Pilot
        sig, bl, freqs = rig.measure_trial(f"PILOT_{pilot_idx}", pilot_wf, {"type": "pilot", "idx": pilot_idx})
        completed += 1
        mask = (freqs >= SUBH_LOW_HZ) & (freqs <= SUBH_HIGH_HZ)
        print_progress("PILOT_0", "System baseline check", np.nan)

        # ----------------------------------------------------------------------
        # RUN TASK 1: TUNED CHANNEL EQUALIZATION
        # ----------------------------------------------------------------------
        if task in ["all", "eq"]:
            print(f"\n>>> EXECUTING TASK 1: TUNED CHANNEL EQUALIZATION (80 ns Overlap, Target: s(t-1)) <<<")
            accum_X_eq, accum_y_eq = [], []
            accum_X_eq_ctrl, accum_y_eq_ctrl = [], []

            for i, tr in enumerate(eq_trials):
                wf = build_eq_waveform_tuned(tr, scale=1.0)
                key = f"EQ_TUNED_t{tr['trial_id']}_s{''.join(str(int(x>0)) for x in tr['symbols'])}"
                meta = {"block": "EQ_TUNED", "symbols": tr["symbols"], "target_lag1": tr["target_lag1"], "target_lag0": tr["target_lag0"]}
                sig, bl, freqs = rig.measure_trial(key, wf, meta)
                completed += 1

                spec = sig[mask] - bl[mask]
                accum_X_eq.append(spec)
                accum_y_eq.append(tr["target_lag1"])

                cur_acc = pure_numpy_ridge_cv(accum_X_eq, accum_y_eq)
                print_progress(f"EQ Trial {i+1}/{len(eq_trials)}", f"Target s(t-1)={tr['target_lag1']}", cur_acc, target_acc=0.85)

                if (i + 1) % 25 == 0:
                    pilot_idx += 1
                    rig.measure_trial(f"PILOT_{pilot_idx}", pilot_wf, {"type": "pilot", "idx": pilot_idx})
                    completed += 1
                    print_progress(f"PILOT_{pilot_idx}", f"Checkpoint EQ {i+1}", cur_acc)

            # Sub-threshold Controls for EQ
            for j, tr in enumerate(eq_ctrl_trials):
                wf = build_eq_waveform_tuned(tr, scale=0.3)
                key = f"EQ_CTRL03X_t{tr['trial_id']}_s{''.join(str(int(x>0)) for x in tr['symbols'])}"
                meta = {"block": "EQ_CTRL03X", "symbols": tr["symbols"], "target_lag1": tr["target_lag1"]}
                sig, bl, freqs = rig.measure_trial(key, wf, meta)
                completed += 1

                spec = sig[mask] - bl[mask]
                accum_X_eq_ctrl.append(spec)
                accum_y_eq_ctrl.append(tr["target_lag1"])

                cur_ctrl_acc = pure_numpy_ridge_cv(accum_X_eq_ctrl, accum_y_eq_ctrl)
                print_progress(f"EQ Ctrl {j+1}/{len(eq_ctrl_trials)}", "0.3x sub-threshold", cur_ctrl_acc, target_acc=0.85)

        # ----------------------------------------------------------------------
        # RUN TASK 2: SINGLE-SPEAKER SPOKEN DIGIT RECOGNITION (JACKSON 4-CLASS)
        # ----------------------------------------------------------------------
        if task in ["all", "digits"]:
            print(f"\n>>> EXECUTING TASK 2: SINGLE-SPEAKER SPOKEN DIGITS (Jackson 4-Class: '0','1','2','3') <<<")
            accum_X_d, accum_y_d = [], []
            accum_X_d_ctrl, accum_y_d_ctrl = [], []

            for i, dt in enumerate(digit_trials):
                item = dt["item"]
                d = dt["digit"]
                r = dt["rep"]
                dphi = dt["dphi"]
                key = f"DIGIT_JACKSON_d{d}_r{r}_p{int(dphi*180/np.pi)}"
                wf = build_digit_waveform_jackson(item, scale=1.0, dphi=dphi)
                meta = {"block": "DIGIT_JACKSON", "digit": d, "speaker": "jackson", "rep": r, "dphi": dphi}
                sig, bl, freqs = rig.measure_trial(key, wf, meta)
                completed += 1

                spec = sig[mask] - bl[mask]
                accum_X_d.append(spec)
                accum_y_d.append(d)

                cur_acc = pure_numpy_ridge_cv(accum_X_d, accum_y_d)
                print_progress(f"Digit Trial {i+1}/{len(digit_trials)}", f"Digit '{d}' (Jackson)", cur_acc, target_acc=0.85)

                if (i + 1) % 25 == 0:
                    pilot_idx += 1
                    rig.measure_trial(f"PILOT_{pilot_idx}", pilot_wf, {"type": "pilot", "idx": pilot_idx})
                    completed += 1
                    print_progress(f"PILOT_{pilot_idx}", f"Checkpoint Digit {i+1}", cur_acc)

            # Sub-threshold Controls for Digits
            for j, dt in enumerate(digit_ctrl_trials):
                item = dt["item"]
                d = dt["digit"]
                key = f"DIGIT_JACKSON_CTRL03X_d{d}_j{j}"
                wf = build_digit_waveform_jackson(item, scale=0.3)
                meta = {"block": "DIGIT_JACKSON_CTRL03X", "digit": d, "speaker": "jackson", "scale": 0.3}
                sig, bl, freqs = rig.measure_trial(key, wf, meta)
                completed += 1

                spec = sig[mask] - bl[mask]
                accum_X_d_ctrl.append(spec)
                accum_y_d_ctrl.append(d)

                cur_ctrl_acc = pure_numpy_ridge_cv(accum_X_d_ctrl, accum_y_d_ctrl)
                print_progress(f"Digit Ctrl {j+1}/{len(digit_ctrl_trials)}", f"0.3x Ctrl: Digit '{d}'", cur_ctrl_acc, target_acc=0.85)

        # Final Reference Pilot
        pilot_idx += 1
        rig.measure_trial(f"PILOT_{pilot_idx}", pilot_wf, {"type": "pilot", "idx": pilot_idx})
        completed += 1
        print_progress(f"PILOT_{pilot_idx}", "Final run checkpoint", np.nan)

    finally:
        rig.close()

    total_time = time.time() - t0
    print("\n" + "=" * 78)
    print("ALL REQUESTED BENCHMARKS COMPLETED")
    print(f"Total Acquisitions:             {completed}/{total_acquisitions}")
    print(f"Total Execution Time:           {total_time:.1f} s ({total_time / 60.0:.2f} min)")
    print(f"HDF5 File Saved:                {out_h5}")
    print("=" * 78)


def main():
    parser = argparse.ArgumentParser(description="Unified High-Accuracy Neuromorphic Benchmark Runner")
    parser.add_argument("--out", type=str, default=DEFAULT_OUT, help="Path to output HDF5 dataset")
    parser.add_argument("--task", type=str, default="all", choices=["all", "eq", "digits"], help="Which task to run")
    parser.add_argument("--dry-run", action="store_true", help="Simulate without hardware connection")
    parser.add_argument("--n-eq", type=int, default=80, help="Number of standard EQ trials")
    parser.add_argument("--n-eq-ctrl", type=int, default=20, help="Number of 0.3x EQ control trials")
    parser.add_argument("--digit-reps", type=int, default=3, help="Repetitions for 20 Jackson recordings (e.g. 3 => 60 trials)")
    parser.add_argument("--n-digit-ctrl", type=int, default=16, help="Number of 0.3x digit controls (4 per class)")
    args = parser.parse_args()

    run_benchmarks(
        out_h5=args.out,
        task=args.task,
        dry_run=args.dry_run,
        n_eq=args.n_eq,
        n_eq_ctrl=args.n_eq_ctrl,
        digit_reps=args.digit_reps,
        n_digit_ctrl=args.n_digit_ctrl,
    )


if __name__ == "__main__":
    main()
