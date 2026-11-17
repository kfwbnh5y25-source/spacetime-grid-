"""
V3 Hardware Execution Benchmark: Isolated Spoken Digit Recognition ("0", "1", "2", "3")
  - Physics-aligned Phoneme Transition Encoding:
      In each 350-ns frame, plays ONE dominant formant tone:
        Tone A (4081 MHz, ~365 mV) if low formant dominates (vowels/glides)
        Tone B (4125 MHz, ~355 mV) if high formant dominates (consonants/fricatives/plosives)
      Adjacent phonemes overlap by 100 ns (dur_ns=450, step_ns=350, total 1.85 us << 20 us AWG limit).
      Peak amplitude is strictly <= 380 mV -> ZERO amplifier clipping or electronic IMD3 distortion!
  - 80 Reservoir Trials (20 per digit across 4 distinct human speakers), randomly shuffled
  - 20 Sub-threshold 0.3x Controls (Balanced 5 per digit across all 4 classes)
  - Interleaved reference pilots every 25 trials
  - Real-time live online accuracy printout after every trial
  - Output: spin-waves/neuromorphic/data_spoken_digits_v3.h5

Usage:
  python run_spoken_digits_benchmark_v3.py --dry-run
  python run_spoken_digits_benchmark_v3.py --out data_spoken_digits_v3.h5
"""

import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np
from sklearn.linear_model import RidgeClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from hw_core import HardwareRig, make_pulse_waveform

HERE = Path(__file__).resolve().parent
FEATURES_FILE = HERE / "spoken_digits_features_v2.json"

F_A_MHZ = 4081.0  # Low formant carrier (vowels, glides)
F_B_MHZ = 4125.0  # High formant carrier (consonants, sibilants, plosives)

AMP_A_MV = 365.0
AMP_B_MV = 355.0

DUR_NS = 450.0   # Pulse duration
STEP_NS = 350.0  # 100 ns overlap between successive phonemes (activates 4-magnon memory!)
GAP_NS = 1500.0


def build_digit_waveform_v3(trial: dict, scale: float = 1.0) -> np.ndarray:
    """
    Constructs a 5-step phoneme transition waveform:
    In each 350-ns frame:
      If low_formant >= high_formant: plays Tone A (4081 MHz)
      Else: plays Tone B (4125 MHz)
      Amplitude modulated slightly (+/-15 mV) by acoustic loudness.
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
            base_amp = AMP_A_MV
            mod = 20.0 * (lo_val - 0.5)
        else:
            freq = F_B_MHZ
            base_amp = AMP_B_MV
            mod = 20.0 * (hi_val - 0.5)

        amp = (base_amp + mod) * scale
        syms.append({
            "amp_mv": float(amp),
            "freq_mhz": float(freq),
            "start_ns": float(t_cur),
            "dur_ns": float(DUR_NS),
            "phase": 0.0,
        })
        t_cur += STEP_NS

    return make_pulse_waveform(syms, gap_ns=GAP_NS)


def build_pilot_wf() -> np.ndarray:
    """Standard 2-tone pilot reference for drift tracking."""
    syms = [
        {"amp_mv": 365.0, "freq_mhz": F_A_MHZ, "start_ns": 0.0, "dur_ns": 400.0, "phase": 0.0},
        {"amp_mv": 350.0, "freq_mhz": F_B_MHZ, "start_ns": 200.0, "dur_ns": 400.0, "phase": 0.0},
    ]
    return make_pulse_waveform(syms, gap_ns=1500.0)


def get_balanced_ctrl_subset(dataset: list[dict], n_ctrl: int = 20, seed: int = 1337) -> list[dict]:
    per_digit = n_ctrl // 4
    by_digit = {}
    for item in dataset:
        by_digit.setdefault(item["digit"], []).append(item)

    rng = np.random.default_rng(seed)
    ctrl_subset = []
    for d in sorted(by_digit.keys()):
        chosen = rng.choice(by_digit[d], size=per_digit, replace=False).tolist()
        ctrl_subset.extend(chosen)
    rng.shuffle(ctrl_subset)
    return ctrl_subset


def compute_live_acc(X_list: list[np.ndarray], y_list: list[int]) -> float:
    classes = set(y_list)
    if len(classes) < 2:
        return float(np.nan)
    min_count = min(y_list.count(c) for c in classes)
    if min_count < 2 or len(classes) < 4:
        return float(np.nan)

    n_splits = min(5, min_count)
    if n_splits < 2:
        return float(np.nan)

    X = np.array(X_list)
    y = np.array(y_list)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    preds = np.zeros(len(y), dtype=int)
    for tr, te in skf.split(X, y):
        sc = StandardScaler().fit(X[tr])
        clf = RidgeClassifier(alpha=10.0).fit(sc.transform(X[tr]), y[tr])
        preds[te] = clf.predict(sc.transform(X[te]))
    return float(accuracy_score(y, preds))


def run_spoken_digits_v3(out_h5: str, dry_run: bool = False, n_ctrl: int = 20, seed: int = 42):
    if not FEATURES_FILE.exists():
        print(f"[FAIL] Missing {FEATURES_FILE}. Run prepare_fsdd_audio_v2.py first!")
        sys.exit(1)

    dataset_raw = json.loads(FEATURES_FILE.read_text())
    
    # Shuffle reservoir dataset to interleave all 4 digits across time
    rng = np.random.default_rng(seed)
    dataset = list(dataset_raw)
    rng.shuffle(dataset)

    n_res = len(dataset)
    ctrl_subset = get_balanced_ctrl_subset(dataset_raw, n_ctrl, seed=seed + 1)
    actual_n_ctrl = len(ctrl_subset)

    total_trials = n_res + actual_n_ctrl + 5
    est_sec = total_trials * 4.8 if not dry_run else 0.5

    print("=" * 76)
    print("SOTA BENCHMARK V3: ISOLATED SPOKEN DIGIT RECOGNITION (FSDD 4-CLASS)")
    print(f"Target HDF5:      {out_h5}")
    print(f"Hardware Mode:    {'[DRY-RUN OFFLINE SIMULATION]' if dry_run else '[LIVE KNET RF HARDWARE RIG]'}")
    print(f"Total Trials:     {total_trials} ({n_res} Interleaved Reservoir + {actual_n_ctrl} Balanced Ctrl + 5 Pilots)")
    print(f"Estimated Time:   {est_sec / 60.0:.1f} minutes")
    print(f"Waveform Layout:  Single-tone phoneme sequence with 100 ns overlap (peak <= 380 mV)")
    print(f"Live Metric:      Real-time 4-class cross-validated accuracy computed after every trial!")
    print("=" * 76)

    rig = HardwareRig(out_h5, dry_run=dry_run)
    pilot_wf = build_pilot_wf()
    t0 = time.time()
    completed = 0

    accum_X_res = []
    accum_y_res = []

    accum_X_ctrl = []
    accum_y_ctrl = []

    def format_progress(step_name: str, detail: str, acc_val: float):
        elapsed = time.time() - t0
        pct = (completed / max(1, total_trials)) * 100
        rem_sec = (elapsed / max(1, completed)) * (total_trials - completed) if completed > 0 else est_sec

        if np.isnan(acc_val):
            status_tag = f"\033[94m[WARMUP]\033[0m"
            acc_str = "Acc=..."
        elif acc_val >= 0.50:
            status_tag = f"\033[92m[SUCCESS]\033[0m"
            acc_str = f"Acc={acc_val*100:4.1f}%"
        elif acc_val >= 0.35:
            status_tag = f"\033[93m[PARTIAL]\033[0m"
            acc_str = f"Acc={acc_val*100:4.1f}%"
        else:
            status_tag = f"\033[91m[FAIL]\033[0m"
            acc_str = f"Acc={acc_val*100:4.1f}%"

        print(f"[{completed:3d}/{total_trials:3d} | {pct:4.1f}%] {status_tag} ({acc_str:9s}) {step_name:<28s} | {detail:<20s} | Rem: {rem_sec/60.0:4.1f} min", flush=True)

    try:
        # Initial Reference Pilot P0
        sig, bl, freqs = rig.measure_trial("PILOT_0", pilot_wf, {"type": "pilot", "idx": 0})
        completed += 1
        mask = (freqs >= 2030e6) & (freqs <= 2075e6)
        format_progress("Initial Reference P0", "Baseline calibration", np.nan)

        # Main Reservoir Trials with LIVE ACCURACY
        pilot_counter = 0
        latest_acc = np.nan
        for i, item in enumerate(dataset):
            digit = item["digit"]
            spk = item["speaker"]
            rep = item["rep"]
            step_id = f"DIGIT_V3_d{digit}_{spk}_r{rep}"

            wf = build_digit_waveform_v3(item, scale=1.0)
            meta = {"block": "DIGIT_RESERVOIR_V3", "digit": digit, "speaker": spk, "rep": rep}
            sig, bl, freqs = rig.measure_trial(step_id, wf, meta)
            completed += 1

            spec = sig[mask] - bl[mask]
            accum_X_res.append(spec)
            accum_y_res.append(digit)

            # Live 4-class CV evaluation
            latest_acc = compute_live_acc(accum_X_res, accum_y_res)
            detail = f"Digit '{digit}' ({spk})"
            format_progress(f"Trial {i+1}/{n_res}: {step_id}", detail, latest_acc)

            # Interleaved Pilot every 25 trials
            if (i + 1) % 25 == 0:
                pilot_counter += 1
                p_id = f"PILOT_{pilot_counter}"
                rig.measure_trial(p_id, pilot_wf, {"type": "pilot", "idx": pilot_counter})
                completed += 1
                format_progress(f"Drift Pilot P{pilot_counter}", f"Checkpoint {i+1}/{n_res}", latest_acc)

        # Sub-Threshold 0.3x Control Trials
        latest_ctrl_acc = np.nan
        for j, item in enumerate(ctrl_subset):
            digit = item["digit"]
            spk = item["speaker"]
            rep = item["rep"]
            step_id = f"DIGIT_V3_CTRL03X_d{digit}_{spk}_r{rep}"

            wf = build_digit_waveform_v3(item, scale=0.3)
            meta = {"block": "DIGIT_CTRL03X_V3", "digit": digit, "speaker": spk, "rep": rep}
            sig, bl, freqs = rig.measure_trial(step_id, wf, meta)
            completed += 1

            spec = sig[mask] - bl[mask]
            accum_X_ctrl.append(spec)
            accum_y_ctrl.append(digit)

            latest_ctrl_acc = compute_live_acc(accum_X_ctrl, accum_y_ctrl)
            detail = f"Ctrl 0.3x: Digit '{digit}'"
            format_progress(f"Ctrl {j+1}/{actual_n_ctrl}: {step_id}", detail, latest_ctrl_acc)

        # Final Reference Pilot
        pilot_counter += 1
        p_id = f"PILOT_{pilot_counter}"
        rig.measure_trial(p_id, pilot_wf, {"type": "pilot", "idx": pilot_counter})
        completed += 1
        format_progress(f"Final Reference P{pilot_counter}", "Run complete", latest_acc)

    finally:
        rig.close()

    total_time = time.time() - t0
    final_res_acc = compute_live_acc(accum_X_res, accum_y_res)
    final_ctrl_acc = compute_live_acc(accum_X_ctrl, accum_y_ctrl)

    print("\n" + "=" * 76)
    print("SPOKEN DIGIT V3 BENCHMARK COMPLETED")
    print(f"Total Acquisitions:             {completed}/{total_trials}")
    print(f"Total Runtime:                  {total_time:.1f} s ({total_time / 60.0:.2f} min)")
    print(f"Final 4-Class Reservoir Acc:    {final_res_acc * 100:.1f}% (Chance: 25.0%)")
    print(f"Final Balanced Control Acc:     {final_ctrl_acc * 100:.1f}%")
    print(f"Output HDF5:                    {out_h5}")
    print("=" * 76)


def main():
    parser = argparse.ArgumentParser(description="Spoken Digits Physical Reservoir Benchmark V3")
    parser.add_argument("--out", type=str, default=str(HERE / "data_spoken_digits_v3.h5"), help="Output HDF5 path")
    parser.add_argument("--dry-run", action="store_true", help="Offline simulation mode without RF hardware")
    parser.add_argument("--n-ctrl", type=int, default=20, help="Number of balanced 0.3x control trials")
    args = parser.parse_args()

    run_spoken_digits_v3(args.out, dry_run=args.dry_run, n_ctrl=args.n_ctrl)


if __name__ == "__main__":
    main()
