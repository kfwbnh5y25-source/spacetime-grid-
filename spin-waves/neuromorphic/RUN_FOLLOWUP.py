"""
One-Click PyCharm Follow-Up Runner (Step 1 + Step 2, ~60 min, ~750 trials).
Fixes E1 (integer-cycle G=0 control), E3 (true inter-symbol delay [A][d][B] at fixed T_L),
and runs E5 (Night-2 cross-session generalization of E2_ovl).

Output: spin-waves/neuromorphic/data_followup.h5 (+ .manifest.json)
"""

import argparse
import json
import time
from pathlib import Path
import numpy as np

from hw_core import HardwareRig, make_pulse_waveform, align_samples, FS_HZ

HERE = Path(__file__).resolve().parent
OUT_FILE = str(HERE / "data_followup.h5")


def build_pilot_wf() -> np.ndarray:
    syms = [
        {"amp_mv": 350.0, "freq_mhz": 4095.0, "start_ns": 0.0, "dur_ns": 500.0, "phase": 0.0},
        {"amp_mv": 350.0, "freq_mhz": 4125.0, "start_ns": 250.0, "dur_ns": 500.0, "phase": 0.0},
    ]
    return make_pulse_waveform(syms, gap_ns=1500.0)


def build_e1_fix(n_rep: int = 20) -> list[tuple[str, np.ndarray, dict]]:
    """Integer-cycle frequencies on 1.5 us (2 MHz grid -> 4096, 4110, 4126 MHz) so G=0 has zero wrap-around phase jump."""
    sA = {"freq_mhz": 4096.0, "amp_mv": 350.0, "dur_ns": 500.0}
    sB = {"freq_mhz": 4126.0, "amp_mv": 350.0, "dur_ns": 500.0}
    sC = {"freq_mhz": 4110.0, "amp_mv": 350.0, "dur_ns": 500.0}
    orders = [("ABC", [sA, sB, sC]), ("BCA", [sB, sC, sA]), ("CAB", [sC, sA, sB])]
    trials = []
    for g_ns in [0.0, 500.0, 2000.0]:
        for r in range(n_rep):
            for name, seq in orders:
                syms, t_cur = [], 0.0
                for s in seq:
                    syms.append({**s, "start_ns": t_cur, "phase": 0.0})
                    t_cur += s["dur_ns"]
                wf = make_pulse_waveform(syms, gap_ns=g_ns)
                trials.append((f"E1FIX_{name}_G{int(g_ns)}_r{r}", wf, {"block": "E1FIX", "order": name, "gap_ns": g_ns, "rep": r}))
    return trials


def build_e3_fix(n_rep: int = 15) -> list[tuple[str, np.ndarray, dict]]:
    """True inter-symbol delay [A][d][B] vs [B][d][A] at FIXED T_L = 3500 ns (constant duty cycle)."""
    delays_ns = [0.0, 10.0, 20.0, 40.0, 80.0, 150.0, 300.0]
    fixed_total_ns = 3500.0
    trials = []
    for amp_mv in [350.0, 380.0]:
        sA = {"freq_mhz": 4095.0, "amp_mv": amp_mv, "dur_ns": 500.0}
        sB = {"freq_mhz": 4125.0, "amp_mv": amp_mv, "dur_ns": 500.0}
        for d_ns in delays_ns:
            gap_ns = fixed_total_ns - (1000.0 + d_ns)
            for r in range(n_rep):
                for name, (s1, s2) in [("AB", (sA, sB)), ("BA", (sB, sA))]:
                    syms = [
                        {**s1, "start_ns": 0.0, "phase": 0.0},
                        {**s2, "start_ns": 500.0 + d_ns, "phase": 0.0},
                    ]
                    wf = make_pulse_waveform(syms, gap_ns=gap_ns)
                    trials.append((f"E3FIX_a{int(amp_mv)}_{name}_d{int(d_ns)}_r{r}", wf, {
                        "block": "E3FIX", "amp_mv": amp_mv, "order": name, "delay_ns": d_ns, "rep": r
                    }))
    return trials


def build_e5_session2(n_rep: int = 40, seed: int = 777) -> list[tuple[str, np.ndarray, dict]]:
    """E5: Cross-session generalization test of E2_ovl (+/-3% jitter, +/-3 MHz detuning)."""
    rng = np.random.default_rng(seed)
    trials = []
    for r in range(n_rep):
        amp_jit = 1.0 + rng.uniform(-0.03, 0.03)
        df_mhz = rng.uniform(-3.0, 3.0)
        for name, (f1, f2) in [("AB", (4095.0, 4125.0)), ("BA", (4125.0, 4095.0))]:
            syms = [
                {"amp_mv": 350.0 * amp_jit, "freq_mhz": f1 + df_mhz, "start_ns": 0.0, "dur_ns": 1000.0, "phase": 0.0},
                {"amp_mv": 350.0 * amp_jit, "freq_mhz": f2 + df_mhz, "start_ns": 500.0, "dur_ns": 500.0, "phase": 0.0},
            ]
            wf = make_pulse_waveform(syms, gap_ns=2000.0)
            trials.append((f"E5_ovl_{name}_r{r}", wf, {"block": "E5", "order": name, "rep": r}))
    return trials


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=str, default=OUT_FILE)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    rig = HardwareRig(args.out, dry_run=args.dry_run, script_filename=__file__)
    pilot_wf = build_pilot_wf()

    all_trials = build_e1_fix(20) + build_e3_fix(15) + build_e5_session2(40)
    np.random.default_rng(2027).shuffle(all_trials)

    manifest = []
    t0 = time.time()
    print(f"[FOLLOWUP START] {len(all_trials)} trials (~60 min) | dry_run={args.dry_run}")
    try:
        for idx, (key, wf, meta) in enumerate(all_trials):
            if idx % 10 == 0:
                rig.measure_trial(f"E0_PILOT_{idx // 10}", pilot_wf, {"block": "E0"})
            rig.measure_trial(key, wf, meta)
            manifest.append({"key": key, "t_s": round(time.time() - t0, 2), **meta})
            if (idx + 1) % 50 == 0:
                print(f"  Progress: {idx + 1}/{len(all_trials)} ({time.time() - t0:.1f}s)")
    except KeyboardInterrupt:
        print("[INTERRUPTED] Flushing data safely...")
    finally:
        rig.close()
        Path(args.out + ".manifest.json").write_text(json.dumps(manifest, indent=2))
        print(f"[DONE] Saved {len(manifest)} trials to {args.out}")


if __name__ == "__main__":
    main()
