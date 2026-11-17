"""
~2-Minute Calibration & Sanity Check (~45-50 acquisitions at ~2.4-4.8 s/trial).
Verifies:
  1. EX excitation monitor (4000-4200 MHz) for ZHL linearity & IMD3.
  2. Single-tone f/2 peaks & thresholds across [4035, 4071, 4095, 4125] MHz in 2000-2080 MHz window.
  3. Quick AB vs BA overlap/sequential contrast + 0.3x electronic control.

Usage:
  python spin-waves/neuromorphic/run_2min_calibration.py --out cal_2min.h5 [--dry-run]
"""

import argparse
import time
import numpy as np

from hw_core import HardwareRig, make_pulse_waveform, hampel_despike


def snr_summary(sig: np.ndarray, bl: np.ndarray, freqs_hz: np.ndarray):
    diff = hampel_despike(sig - bl)
    med = np.median(diff)
    sigma = 1.4826 * np.median(np.abs(diff - med)) + 1e-12
    pk_idx = int(np.argmax(diff))
    return float(freqs_hz[pk_idx] / 1e6), float((diff[pk_idx] - med) / sigma), float(np.sum(np.maximum(diff, 0.0)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=str, default="cal_2min.h5")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    t0 = time.time()
    rig = HardwareRig(args.out, dry_run=args.dry_run, script_filename=__file__)

    try:
        # 1. EX excitation monitor (2-tone)
        rig.configure_sa_excitation(center_hz=4100e6, span_hz=200e6)
        wf_ex = make_pulse_waveform([
            {"amp_mv": 260.0, "freq_mhz": 4071.0, "start_ns": 0.0, "dur_ns": 1000.0, "phase": 0.0},
            {"amp_mv": 280.0, "freq_mhz": 4125.0, "start_ns": 500.0, "dur_ns": 500.0, "phase": 0.0},
        ], gap_ns=2000.0)
        s_ex, b_ex, f_ex = rig.measure_trial("CAL_EX_2TONE", wf_ex, {"step": "ex"}, measure_baseline=False)
        pk_ex, _, p_ex = snr_summary(s_ex, b_ex, f_ex)
        print(f"[EX] peak={pk_ex:.2f} MHz | sum={p_ex:.3e}")

        # 2. Single-tone frequency & amplitude calibration (80 MHz span at f/2)
        rig.configure_sa_sample(center_hz=2040e6, span_hz=80e6, rbw_hz=500e3)
        print("\n[SINGLE-TONE CALIBRATION]")
        for f_mhz in [4035.0, 4071.0, 4095.0, 4125.0]:
            row = []
            for a_mv in [220.0, 260.0, 300.0, 340.0]:
                wf = make_pulse_waveform([{"amp_mv": a_mv, "freq_mhz": f_mhz, "start_ns": 0.0, "dur_ns": 1000.0, "phase": 0.0}], gap_ns=2000.0)
                s, b, fr = rig.measure_trial(f"CAL_1T_f{int(f_mhz)}_a{int(a_mv)}", wf, {"f": f_mhz, "a": a_mv})
                pk, snr, p_int = snr_summary(s, b, fr)
                row.append(f"{int(a_mv)}mV:{pk:.1f}MHz(SNR={snr:.1f})")
            print(f"  f={f_mhz:.0f} (f/2={f_mhz/2:.1f}): " + " | ".join(row))

        # 3. Quick 3-rep AB vs BA contrast + 0.3x control at calibrated supra-threshold (4095/4125 MHz, 350 mV)
        print("\n[QUICK ORDER CONTRAST (3 reps each)]")
        for label, scale in [("OVL_350mV", 1.0), ("CTRL_0.3x", 0.3)]:
            specs = {"AB": [], "BA": []}
            for r in range(3):
                for ord_name, (f_a, f_b) in [("AB", (4095.0, 4125.0)), ("BA", (4125.0, 4095.0))]:
                    wf = make_pulse_waveform([
                        {"amp_mv": 350.0 * scale, "freq_mhz": f_a, "start_ns": 0.0, "dur_ns": 1000.0, "phase": 0.0},
                        {"amp_mv": 350.0 * scale, "freq_mhz": f_b, "start_ns": 500.0, "dur_ns": 500.0, "phase": 0.0},
                    ], gap_ns=1500.0)
                    s, b, _ = rig.measure_trial(f"CAL_{label}_{ord_name}_r{r}", wf, {"label": label, "order": ord_name, "r": r})
                    specs[ord_name].append(hampel_despike(s - b))
            A, B = np.array(specs["AB"]), np.array(specs["BA"])
            j_val = float(np.sum((A.mean(0) - B.mean(0)) ** 2) / (np.sum(A.var(0) + B.var(0)) + 1e-12))
            print(f"  {label}: Fisher J(AB, BA) = {j_val:.3f}")

    finally:
        rig.close()
        print(f"\n[DONE] {rig.trial_count} trials in {time.time() - t0:.1f}s -> {args.out}")


if __name__ == "__main__":
    main()
