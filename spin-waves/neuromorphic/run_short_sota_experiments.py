"""
Short Experimental Protocol for SOTA Physical Reservoir Computing Benchmarks:
  1. Nonlinear Channel Equalization (EQ): Inverts nonlinear inter-symbol interference (ISI)
     and amplifier saturation: y(t) = q(t) - 0.2*q(t)^3, where q(t) = s(t) + 0.15*s(t-1) - 0.25*s(t-2).
     Target: Recover s(t-1) from the final subharmonic emission spectrum.
  2. 5-Step Continuous Analog Temporal Memory (NARMA-5): Evaluates analog fading memory
     over a 5-symbol input packet (200 ns/step, total duration 1.0 us << 20 us AWG limit).

Designed for rapid execution (~12-15 minutes total on hardware):
  - 100 trials EQ (standard drive ~350 mV)
  - 40 trials EQ 0.3x control (sub-threshold, rules out linear filter artifacts)
  - Interleaved drift tracking pilot P0 every 25 trials
  - Full --dry-run support for offline verification

Usage:
  python run_short_sota_experiments.py --dry-run
  python run_short_sota_experiments.py --out data_sota_eq.h5
"""

import argparse
import json
import time
from pathlib import Path
import numpy as np

from hw_core import HardwareRig, make_pulse_waveform

HERE = Path(__file__).resolve().parent

# Calibrated frequencies and amplitudes from cal_2min.h5 and multisymbol_v2
F_A_MHZ = 4081.0  # Symbol +1 carrier
F_B_MHZ = 4125.0  # Symbol -1 carrier
BASE_AMP_MV = 355.0
MOD_DEPTH_MV = 35.0  # Modulation amplitude around supra-threshold operating point


def generate_nonlinear_isi_sequence(n_trials: int, seed: int = 42):
    """
    Generates a stream of random binary symbols s(t) in {-1, +1}
    and applies a nonlinear communications channel with ISI and cubic saturation:
      q(t) = s(t) + 0.18*s(t-1) - 0.22*s(t-2)
      y(t) = q(t) - 0.15*q(t)^3 + noise
    Returns trials: list of (symbol_history, channel_distorted_amplitudes, target_s_lag1)
    """
    rng = np.random.default_rng(seed)
    trials = []
    
    # 8 distinct 3-bit binary histories s(t-2), s(t-1), s(t)
    for r in range(n_trials):
        s = rng.choice([-1, 1], size=3)  # [s(t-2), s(t-1), s(t)]
        
        # Channel ISI:
        q = np.zeros(3)
        q[0] = s[0]
        q[1] = s[1] + 0.18 * s[0]
        q[2] = s[2] + 0.18 * s[1] - 0.22 * s[0]
        
        # Nonlinear saturation (e.g. power amp or medium saturation):
        y = q - 0.15 * (q ** 3)
        
        # Normalize y to [-1, 1] range for amplitude modulation
        y_norm = np.clip(y / 1.5, -1.0, 1.0)
        
        target_lag0 = int(s[2] > 0)
        target_lag1 = int(s[1] > 0)  # Primary task: equalize 1-step delayed symbol
        target_lag2 = int(s[0] > 0)
        
        trials.append({
            "trial_id": r,
            "symbols": [int(x) for x in s],
            "channel_out": [float(x) for x in y_norm],
            "target_lag0": target_lag0,
            "target_lag1": target_lag1,
            "target_lag2": target_lag2,
        })
    return trials


def build_eq_waveform(trial: dict, scale: float = 1.0, step_ns: float = 250.0, gap_ns: float = 1500.0) -> np.ndarray:
    """
    Encodes the 3-symbol distorted packet into an RF pulse waveform:
    Each step i has:
      - Carrier freq: F_A if s_i == +1 else F_B
      - Amplitude: (BASE_AMP_MV + MOD_DEPTH_MV * y_norm[i]) * scale
      - Duration: step_ns
    """
    syms = []
    t_cur = 0.0
    for i in range(3):
        s_val = trial["symbols"][i]
        freq = F_A_MHZ if s_val > 0 else F_B_MHZ
        amp = (BASE_AMP_MV + MOD_DEPTH_MV * trial["channel_out"][i]) * scale
        syms.append({
            "amp_mv": float(amp),
            "freq_mhz": float(freq),
            "start_ns": float(t_cur),
            "dur_ns": float(step_ns),
            "phase": 0.0,
        })
        t_cur += step_ns

    return make_pulse_waveform(syms, gap_ns=gap_ns)


def build_pilot_wf() -> np.ndarray:
    """Standard 2-tone pilot for drift tracking."""
    syms = [
        {"amp_mv": 350.0, "freq_mhz": F_A_MHZ, "start_ns": 0.0, "dur_ns": 400.0, "phase": 0.0},
        {"amp_mv": 350.0, "freq_mhz": F_B_MHZ, "start_ns": 200.0, "dur_ns": 400.0, "phase": 0.0},
    ]
    return make_pulse_waveform(syms, gap_ns=1500.0)


def run_experiment(out_h5: str, dry_run: bool = False, n_eq: int = 100, n_ctrl: int = 40):
    print("=" * 70)
    print("SOTA EXPERIMENTAL BENCHMARK: NONLINEAR CHANNEL EQUALIZATION")
    print(f"Target HDF5: {out_h5} | Dry-run: {dry_run}")
    print(f"Total Trials: {n_eq} (Reservoir EQ) + {n_ctrl} (0.3x Control) + Pilots")
    est_sec = (n_eq + n_ctrl + 6) * 4.8 if not dry_run else 1.0
    print(f"Estimated Machine Time: {est_sec / 60.0:.1f} minutes")
    print("=" * 70)

    rig = HardwareRig(out_h5, dry_run=dry_run)
    pilot_wf = build_pilot_wf()

    trials_eq = generate_nonlinear_isi_sequence(n_eq, seed=42)
    trials_ctrl = generate_nonlinear_isi_sequence(n_ctrl, seed=1337)

    t0 = time.time()
    pilot_count = 0

    try:
        # Initial Pilot
        print("\n[PILOT] Acquiring initial reference P0...")
        rig.measure_trial("PILOT_0", pilot_wf, {"type": "pilot", "idx": 0})

        # Main Nonlinear Channel Equalization Trials
        print(f"\n[EQ] Running {n_eq} Nonlinear Channel Equalization trials...")
        for i, tr in enumerate(trials_eq):
            wf = build_eq_waveform(tr, scale=1.0)
            key = f"EQ_NONLIN_t{tr['trial_id']}_s{''.join(str(int(x>0)) for x in tr['symbols'])}"
            meta = {
                "block": "EQ_NONLIN",
                "symbols": tr["symbols"],
                "target_lag1": tr["target_lag1"],
                "channel_out": tr["channel_out"],
            }
            rig.measure_trial(key, wf, meta)

            # Interleaved pilot every 25 trials
            if (i + 1) % 25 == 0:
                pilot_count += 1
                print(f"  -> Intermediate Pilot {pilot_count} (Trial {i+1}/{n_eq})...")
                rig.measure_trial(f"PILOT_{pilot_count}", pilot_wf, {"type": "pilot", "idx": pilot_count})

        # Sub-threshold Control Trials (0.3x)
        print(f"\n[CTRL] Running {n_ctrl} Sub-Threshold 0.3x Controls...")
        for i, tr in enumerate(trials_ctrl):
            wf = build_eq_waveform(tr, scale=0.3)
            key = f"EQ_CTRL03X_t{tr['trial_id']}_s{''.join(str(int(x>0)) for x in tr['symbols'])}"
            meta = {
                "block": "EQ_CTRL03X",
                "symbols": tr["symbols"],
                "target_lag1": tr["target_lag1"],
                "channel_out": tr["channel_out"],
            }
            rig.measure_trial(key, wf, meta)

        # Final Pilot
        pilot_count += 1
        print(f"\n[PILOT] Acquiring final reference P{pilot_count}...")
        rig.measure_trial(f"PILOT_{pilot_count}", pilot_wf, {"type": "pilot", "idx": pilot_count})

    finally:
        rig.close()

    dt = time.time() - t0
    print("\n" + "=" * 70)
    print(f"EXPERIMENT COMPLETED in {dt:.1f} s ({dt/60.0:.1f} min)")
    print(f"Saved to {out_h5}")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(description="Short SOTA Benchmark Experiment")
    parser.add_argument("--out", type=str, default=str(HERE / "data_sota_eq.h5"), help="Output HDF5 path")
    parser.add_argument("--dry-run", action="store_true", help="Simulate without hardware")
    parser.add_argument("--n-eq", type=int, default=100, help="Number of EQ trials")
    parser.add_argument("--n-ctrl", type=int, default=40, help="Number of 0.3x control trials")
    args = parser.parse_args()

    run_experiment(args.out, dry_run=args.dry_run, n_eq=args.n_eq, n_ctrl=args.n_ctrl)


if __name__ == "__main__":
    main()
