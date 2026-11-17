"""
One-Click PyCharm Runner (v2) for Multi-Symbol Permutation (S_3), 3-Bit Words, and Phase-Gating.
Fixes from v1 audit:
  - Removes +/-2 MHz random detuning (which scrambled 4B relative phase by +/-2.4*pi).
  - Increases pulse duration (dur_ns=600, step_ns=450 -> 150 ns overlap) and equalizes
    symbol overdrive (4081 MHz at 375 mV, 4095 MHz at 355 mV, 4125 MHz at 350 mV; max peak <= 0.73 V).
  - Adds PHASE_GATE sweep (dphi in {0, pi/4, pi/2, 3pi/4, pi}) to directly measure 4B phase sensitivity.

Output: spin-waves/neuromorphic/data_multisymbol_v2.h5 (+ .manifest.json)
Duration: ~50 min (~610 trials + 61 pilots).
"""

import argparse
import itertools
import json
import time
from pathlib import Path
import numpy as np

from hw_core import HardwareRig, make_pulse_waveform

HERE = Path(__file__).resolve().parent
OUT_FILE = str(HERE / "data_multisymbol_v2.h5")

TONES = {
    "A": (4081.0, 365.0),
    "B": (4095.0, 355.0),
    "C": (4125.0, 350.0),
    "0": (4081.0, 365.0),
    "1": (4125.0, 350.0),
}


def build_pilot_wf() -> np.ndarray:
    syms = [
        {"amp_mv": 365.0, "freq_mhz": 4081.0, "start_ns": 0.0, "dur_ns": 500.0, "phase": 0.0},
        {"amp_mv": 350.0, "freq_mhz": 4125.0, "start_ns": 250.0, "dur_ns": 500.0, "phase": 0.0},
    ]
    return make_pulse_waveform(syms, gap_ns=1500.0)


def build_3sym_seq(tokens: tuple[str, ...], scale: float = 1.0, amp_jit: float = 1.0, dphi: float = 0.0, dur_ns: float = 600.0, step_ns: float = 450.0) -> np.ndarray:
    """3 symbols with 150 ns overlap between neighbors (max 2 tones overlap -> peak <= 745 mV) and controlled phase step dphi."""
    syms = [
        {
            "amp_mv": TONES[tok][1] * scale * amp_jit,
            "freq_mhz": TONES[tok][0],
            "start_ns": i * step_ns,
            "dur_ns": dur_ns,
            "phase": float(i * dphi),
        }
        for i, tok in enumerate(tokens)
    ]
    return make_pulse_waveform(syms, gap_ns=2000.0)


def build_trials(seed: int = 2029):
    rng = np.random.default_rng(seed)
    trials = []
    phases = [0.0, 0.25 * np.pi, 0.5 * np.pi, 0.75 * np.pi, np.pi]

    # Task A: 6-class S_3 permutations of (A, B, C) across varied relative phases
    perms = list(itertools.permutations(["A", "B", "C"]))
    for label, n_rep, scale in [("S3", 30, 1.0), ("S3_03x", 15, 0.3)]:
        for r in range(n_rep):
            amp_jit = 1.0 + rng.uniform(-0.02, 0.02)
            pi_idx = r % len(phases)
            dphi = phases[pi_idx]
            for p in perms:
                name = "".join(p)
                wf = build_3sym_seq(p, scale=scale, amp_jit=amp_jit, dphi=dphi)
                trials.append((f"{label}_{name}_p{pi_idx}_r{r}", wf, {"task": label, "word": name, "dphi": float(dphi), "rep": r}))

    # Task B: 8-class 3-bit binary words {0,1}^3 across varied relative phases
    words = list(itertools.product(["0", "1"], repeat=3))
    for label, n_rep, scale in [("BIT3", 25, 1.0), ("BIT3_03x", 10, 0.3)]:
        for r in range(n_rep):
            amp_jit = 1.0 + rng.uniform(-0.02, 0.02)
            pi_idx = r % len(phases)
            dphi = phases[pi_idx]
            for w in words:
                name = "".join(w)
                parity = sum(int(b) for b in w) % 2
                wf = build_3sym_seq(w, scale=scale, amp_jit=amp_jit, dphi=dphi)
                trials.append((f"{label}_{name}_p{pi_idx}_r{r}", wf, {"task": label, "word": name, "parity": parity, "dphi": float(dphi), "rep": r}))

    # Task C: Explicit 4B Phase-Gating Sweep (AB vs BA across 5 relative phases)
    for pi_idx, dphi in enumerate(phases):
        for r in range(6):
            for ord_name, (t1, t2) in [("AB", ("0", "1")), ("BA", ("1", "0"))]:
                syms = [
                    {"amp_mv": TONES[t1][1], "freq_mhz": TONES[t1][0], "start_ns": 0.0, "dur_ns": 800.0, "phase": 0.0},
                    {"amp_mv": TONES[t2][1], "freq_mhz": TONES[t2][0], "start_ns": 400.0, "dur_ns": 500.0, "phase": float(dphi)},
                ]
                wf = make_pulse_waveform(syms, gap_ns=2000.0)
                trials.append((f"PH_{ord_name}_p{pi_idx}_r{r}", wf, {"task": "PHASE_GATE", "order": ord_name, "dphi": float(dphi), "rep": r}))

    return trials


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=str, default=OUT_FILE)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    rig = HardwareRig(args.out, dry_run=args.dry_run, script_filename=__file__)
    pilot_wf = build_pilot_wf()
    all_trials = build_trials()
    np.random.default_rng(1001).shuffle(all_trials)

    manifest = []
    t0 = time.time()
    print(f"[MULTISYMBOL v2 START] {len(all_trials)} trials (~50 min) | dry_run={args.dry_run}")
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
