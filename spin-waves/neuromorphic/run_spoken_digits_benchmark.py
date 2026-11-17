"""
Hardware Execution Benchmark: Isolated Spoken Digit Recognition ("0", "1", "2", "3")
  - Uses 80 human recordings from Free Spoken Digit Dataset (FSDD)
  - 4 classes: Digits '0' (zero), '1' (one), '2' (two), '3' (three)
  - 4 distinct speakers (Jackson, Nicolas, Theo, Yweweler)
  - 5-step dual-formant RF pulse waveform:
      Tone A (4081 MHz) = low formant envelope (200-1200 Hz)
      Tone B (4125 MHz) = high formant envelope (1200-4000 Hz)
  - Drive: 325-365 mV supra-threshold memory regime
  - Controls: 20 sub-threshold 0.3x trials + 4 reference drift pilots
  - Real-time progress display with [SUCCESS] / [FAIL] status tags

Usage:
  python run_spoken_digits_benchmark.py --dry-run
  python run_spoken_digits_benchmark.py --out data_spoken_digits.h5
"""

import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np

from hw_core import HardwareRig, make_pulse_waveform

HERE = Path(__file__).resolve().parent
FEATURES_FILE = HERE / "spoken_digits_features.json"

F_A_MHZ = 4081.0  # Low formant carrier
F_B_MHZ = 4125.0  # High formant carrier
BASE_AMP_MV = 325.0
AMP_MOD_MV = 40.0   # Maps [0, 1] envelope to [325, 365] mV (calibrated 4-magnon memory window)
STEP_DUR_NS = 250.0  # 5 steps * 250 ns = 1.25 us << 20 us AWG limit
GAP_NS = 1500.0


def build_digit_waveform(trial: dict, scale: float = 1.0) -> np.ndarray:
    """
    Constructs a 5-step dual-tone RF pulse waveform from the dual-formant envelopes:
    In each step:
      Tone A (4081 MHz) amp = (BASE_AMP_MV + AMP_MOD_MV * low_formant[i]) * scale
      Tone B (4125 MHz) amp = (BASE_AMP_MV + AMP_MOD_MV * high_formant[i]) * scale
    """
    low_env = trial["low_formant"]
    high_env = trial["high_formant"]
    n_steps = len(low_env)

    syms = []
    t_cur = 0.0
    for i in range(n_steps):
        amp_a = (BASE_AMP_MV + AMP_MOD_MV * low_env[i]) * scale
        amp_b = (BASE_AMP_MV + AMP_MOD_MV * high_env[i]) * scale

        syms.append({
            "amp_mv": float(amp_a),
            "freq_mhz": float(F_A_MHZ),
            "start_ns": float(t_cur),
            "dur_ns": float(STEP_DUR_NS),
            "phase": 0.0,
        })
        syms.append({
            "amp_mv": float(amp_b),
            "freq_mhz": float(F_B_MHZ),
            "start_ns": float(t_cur),
            "dur_ns": float(STEP_DUR_NS),
            "phase": 0.0,
        })
        t_cur += STEP_DUR_NS

    return make_pulse_waveform(syms, gap_ns=GAP_NS)


def build_pilot_wf() -> np.ndarray:
    """Standard 2-tone pilot reference for drift tracking."""
    syms = [
        {"amp_mv": 350.0, "freq_mhz": F_A_MHZ, "start_ns": 0.0, "dur_ns": 400.0, "phase": 0.0},
        {"amp_mv": 350.0, "freq_mhz": F_B_MHZ, "start_ns": 200.0, "dur_ns": 400.0, "phase": 0.0},
    ]
    return make_pulse_waveform(syms, gap_ns=1500.0)


def run_spoken_digits_benchmark(out_h5: str, dry_run: bool = False, n_ctrl: int = 20):
    if not FEATURES_FILE.exists():
        print(f"[FAIL] Missing {FEATURES_FILE}. Run prepare_fsdd_audio.py first!")
        sys.exit(1)

    dataset = json.loads(FEATURES_FILE.read_text())
    n_res = len(dataset)
    total_trials = n_res + n_ctrl + 4  # 80 reservoir + 20 control + 4 pilots = 104 trials
    est_sec = total_trials * 4.8 if not dry_run else 0.5

    print("=" * 72)
    print("SOTA BENCHMARK: ISOLATED SPOKEN DIGIT RECOGNITION (FSDD 4-CLASS)")
    print(f"Target HDF5:      {out_h5}")
    print(f"Hardware Mode:    {'[DRY-RUN OFFLINE SIMULATION]' if dry_run else '[LIVE KNET RF HARDWARE RIG]'}")
    print(f"Total Trials:     {total_trials} ({n_res} Reservoir + {n_ctrl} Control + 4 Pilots)")
    print(f"Estimated Time:   {est_sec / 60.0:.1f} minutes")
    print("=" * 72)

    rig = HardwareRig(out_h5, dry_run=dry_run)
    pilot_wf = build_pilot_wf()
    t0 = time.time()
    completed = 0

    def print_progress(status: str, step_name: str, detail: str = ""):
        elapsed = time.time() - t0
        pct = (completed / max(1, total_trials)) * 100
        rem_sec = (elapsed / max(1, completed)) * (total_trials - completed) if completed > 0 else est_sec
        status_tag = f"\033[92m[SUCCESS]\033[0m" if status == "SUCCESS" else f"\033[91m[FAIL]\033[0m"
        print(f"[{completed:3d}/{total_trials:3d} | {pct:4.1f}%] {status_tag} {step_name:<28s} | {detail:<22s} | Rem: {rem_sec/60.0:4.1f} min", flush=True)

    try:
        # Initial Reference Pilot P0
        step_id = "PILOT_0"
        try:
            rig.measure_trial(step_id, pilot_wf, {"type": "pilot", "idx": 0})
            completed += 1
            print_progress("SUCCESS", "Initial Reference P0", "Baseline calibration")
        except Exception as e:
            print_progress("FAIL", "Initial Reference P0", str(e)[:20])
            raise

        # Main Reservoir Trials
        pilot_counter = 0
        for i, item in enumerate(dataset):
            digit = item["digit"]
            spk = item["speaker"]
            rep = item["rep"]
            step_id = f"DIGIT_d{digit}_{spk}_r{rep}"

            try:
                wf = build_digit_waveform(item, scale=1.0)
                meta = {"block": "DIGIT_RESERVOIR", "digit": digit, "speaker": spk, "rep": rep}
                rig.measure_trial(step_id, wf, meta)
                completed += 1
                detail = f"Digit '{digit}' ({spk})"
                print_progress("SUCCESS", f"Trial {i+1}/{n_res}: {step_id}", detail)
            except Exception as e:
                print_progress("FAIL", f"Trial {i+1}/{n_res}: {step_id}", str(e)[:20])
                raise

            # Interleaved Pilot every 25 trials
            if (i + 1) % 25 == 0:
                pilot_counter += 1
                p_id = f"PILOT_{pilot_counter}"
                try:
                    rig.measure_trial(p_id, pilot_wf, {"type": "pilot", "idx": pilot_counter})
                    completed += 1
                    print_progress("SUCCESS", f"Drift Pilot P{pilot_counter}", f"Checkpoint {i+1}/{n_res}")
                except Exception as e:
                    print_progress("FAIL", f"Drift Pilot P{pilot_counter}", str(e)[:20])
                    raise

        # Sub-Threshold 0.3x Control Trials
        ctrl_subset = dataset[:n_ctrl]
        for j, item in enumerate(ctrl_subset):
            digit = item["digit"]
            spk = item["speaker"]
            rep = item["rep"]
            step_id = f"DIGIT_CTRL03X_d{digit}_{spk}_r{rep}"

            try:
                wf = build_digit_waveform(item, scale=0.3)
                meta = {"block": "DIGIT_CTRL03X", "digit": digit, "speaker": spk, "rep": rep}
                rig.measure_trial(step_id, wf, meta)
                completed += 1
                detail = f"Control 0.3x: Digit '{digit}'"
                print_progress("SUCCESS", f"Ctrl {j+1}/{n_ctrl}: {step_id}", detail)
            except Exception as e:
                print_progress("FAIL", f"Ctrl {j+1}/{n_ctrl}: {step_id}", str(e)[:20])
                raise

        # Final Reference Pilot
        pilot_counter += 1
        p_id = f"PILOT_{pilot_counter}"
        try:
            rig.measure_trial(p_id, pilot_wf, {"type": "pilot", "idx": pilot_counter})
            completed += 1
            print_progress("SUCCESS", f"Final Reference P{pilot_counter}", "Run complete")
        except Exception as e:
            print_progress("FAIL", f"Final Reference P{pilot_counter}", str(e)[:20])
            raise

    finally:
        rig.close()

    total_time = time.time() - t0
    print("\n" + "=" * 72)
    print(f"\033[92m[SUCCESS]\033[0m SPOKEN DIGIT BENCHMARK COMPLETED SUCCESSFULLY!")
    print(f"Total Acquisitions: {completed}/{total_trials}")
    print(f"Total Runtime:      {total_time:.1f} s ({total_time / 60.0:.2f} min)")
    print(f"Output HDF5:        {out_h5}")
    print("=" * 72)


def main():
    parser = argparse.ArgumentParser(description="Spoken Digits Physical Reservoir Benchmark")
    parser.add_argument("--out", type=str, default=str(HERE / "data_spoken_digits.h5"), help="Output HDF5 path")
    parser.add_argument("--dry-run", action="store_true", help="Offline simulation mode without RF hardware")
    parser.add_argument("--n-ctrl", type=int, default=20, help="Number of 0.3x control trials")
    args = parser.parse_args()

    run_spoken_digits_benchmark(args.out, dry_run=args.dry_run, n_ctrl=args.n_ctrl)


if __name__ == "__main__":
    main()
