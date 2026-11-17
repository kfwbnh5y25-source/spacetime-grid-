"""
Experimental Verification Suite (Blocks E0 - E4) for Spin-Wave Neuromorphic Reservoir.
Executes overnight (~10,000 acquisitions, ~8 h) or in --dry-run mode for offline verification.

Usage:
  python run_e0_e4_verification.py --out data_e0_e4.h5 [--dry-run] [--blocks E0 E1 E2 E3 E4]
"""

import argparse
import json
import time
from pathlib import Path
import numpy as np

from hw_core import HardwareRig, make_pulse_waveform


# Calibrated symbols from cal_2min.h5 (4095 & 4125 MHz cross threshold cleanly at ~330-350 mV)
SYM_A = {"freq_mhz": 4095.0, "amp_mv": 350.0, "dur_ns": 500.0}
SYM_B = {"freq_mhz": 4125.0, "amp_mv": 350.0, "dur_ns": 500.0}
SYM_C = {"freq_mhz": 4110.0, "amp_mv": 350.0, "dur_ns": 500.0}
PROBE = {"freq_mhz": 4080.0, "amp_mv": 350.0, "dur_ns": 400.0}


def build_pilot_wf() -> np.ndarray:
    """Fixed 2-tone pilot P0 for minute/hour drift tracking."""
    syms = [
        {"amp_mv": 350.0, "freq_mhz": 4095.0, "start_ns": 0.0, "dur_ns": 500.0, "phase": 0.0},
        {"amp_mv": 350.0, "freq_mhz": 4125.0, "start_ns": 250.0, "dur_ns": 500.0, "phase": 0.0},
    ]
    return make_pulse_waveform(syms, gap_ns=1500.0)


def build_e1_trials(n_rep: int = 30) -> list[tuple[str, np.ndarray, dict]]:
    """E1: Cyclic shift & recovery gap control: ABC, BCA, CAB across G in {0, 250, 500, 1000, 2000, 4000, 8000} ns."""
    trials = []
    gaps_ns = [0.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0]
    orders = [("ABC", [SYM_A, SYM_B, SYM_C]), ("BCA", [SYM_B, SYM_C, SYM_A]), ("CAB", [SYM_C, SYM_A, SYM_B])]
    for g in gaps_ns:
        for r in range(n_rep):
            for name, seq in orders:
                syms = []
                t_cur = 0.0
                for s in seq:
                    syms.append({**s, "start_ns": t_cur, "phase": 0.0})
                    t_cur += s["dur_ns"]
                wf = make_pulse_waveform(syms, gap_ns=g)
                meta = {"block": "E1", "order": name, "gap_ns": g, "rep": r}
                trials.append((f"E1_{name}_G{int(g)}_r{r}", wf, meta))
    return trials


def build_e2_trials(n_rep: int = 80, seed: int = 42) -> list[tuple[str, np.ndarray, dict]]:
    """
    E2: Energy-matched binary order gate (AB vs BA) under +/-3% amp jitter and +/-3 MHz detuning.
    Tests both sequential (delta in {0, 50, 200} ns) and overlap (A -> A+B vs B -> A+B) layouts.
    """
    rng = np.random.default_rng(seed)
    trials = []
    for r in range(n_rep):
        amp_jit = 1.0 + rng.uniform(-0.03, 0.03)
        df_mhz = rng.uniform(-3.0, 3.0)

        # Sequential non-overlap
        for delta_ns in [0.0, 50.0, 200.0]:
            for name, (s1, s2) in [("AB", (SYM_A, SYM_B)), ("BA", (SYM_B, SYM_A))]:
                syms = [
                    {"amp_mv": s1["amp_mv"] * amp_jit, "freq_mhz": s1["freq_mhz"] + df_mhz, "start_ns": 0.0, "dur_ns": s1["dur_ns"], "phase": 0.0},
                    {"amp_mv": s2["amp_mv"] * amp_jit, "freq_mhz": s2["freq_mhz"] + df_mhz, "start_ns": s1["dur_ns"] + delta_ns, "dur_ns": s2["dur_ns"], "phase": 0.0},
                ]
                wf = make_pulse_waveform(syms, gap_ns=2000.0)
                meta = {"block": "E2", "layout": f"seq_d{int(delta_ns)}", "order": name, "amp_jit": float(amp_jit), "df_mhz": float(df_mhz), "rep": r}
                trials.append((f"E2_seq_d{int(delta_ns)}_{name}_r{r}", wf, meta))

        # Overlap layout: 500 ns lead + 500 ns simultaneous A+B
        for name, (s_lead, s_sec) in [("AB", (SYM_A, SYM_B)), ("BA", (SYM_B, SYM_A))]:
            syms = [
                {"amp_mv": s_lead["amp_mv"] * amp_jit, "freq_mhz": s_lead["freq_mhz"] + df_mhz, "start_ns": 0.0, "dur_ns": 1000.0, "phase": 0.0},
                {"amp_mv": s_sec["amp_mv"] * amp_jit, "freq_mhz": s_sec["freq_mhz"] + df_mhz, "start_ns": 500.0, "dur_ns": 500.0, "phase": 0.0},
            ]
            wf = make_pulse_waveform(syms, gap_ns=2000.0)
            meta = {"block": "E2", "layout": "overlap", "order": name, "amp_jit": float(amp_jit), "df_mhz": float(df_mhz), "rep": r}
            trials.append((f"E2_ovl_{name}_r{r}", wf, meta))
    return trials


def build_e3_trials(n_rep: int = 40) -> list[tuple[str, np.ndarray, dict]]:
    """E3: Pump-Probe Delayed-Order Memory Curve ([AB] - delay d - [Probe P])."""
    trials = []
    delays_ns = [0.0, 100.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0]
    for d_ns in delays_ns:
        for r in range(n_rep):
            for name, (s1, s2) in [("AB", (SYM_A, SYM_B)), ("BA", (SYM_B, SYM_A))]:
                t_probe = s1["dur_ns"] + s2["dur_ns"] + d_ns
                syms = [
                    {**s1, "start_ns": 0.0, "phase": 0.0},
                    {**s2, "start_ns": s1["dur_ns"], "phase": 0.0},
                    {**PROBE, "start_ns": t_probe, "phase": 0.0},
                ]
                wf = make_pulse_waveform(syms, gap_ns=2000.0)
                meta = {"block": "E3", "order": name, "delay_ns": d_ns, "rep": r}
                trials.append((f"E3_{name}_d{int(d_ns)}_r{r}", wf, meta))
    return trials


def build_e4_trials(n_rep: int = 40, seed: int = 99) -> list[tuple[str, np.ndarray, dict]]:
    """
    E4: Artifact Disproof Controls (all MUST yield chance <= 55-60%):
      - sub03x: 0.3x amplitude (well below parametric threshold -> tests RF/amplifier memory)
      - offres: >4.4 GHz off-resonance pump (same envelope, no spin-wave excitation)
      - randph: independent random carrier phases per trial
    """
    rng = np.random.default_rng(seed)
    trials = []
    for r in range(n_rep):
        for ctrl in ["sub03x", "offres", "randph"]:
            for name, (s1, s2) in [("AB", (SYM_A, SYM_B)), ("BA", (SYM_B, SYM_A))]:
                scale = 0.3 if ctrl == "sub03x" else 1.0
                f_shift = 400.0 if ctrl == "offres" else 0.0
                p1 = float(rng.uniform(0, 2 * np.pi)) if ctrl == "randph" else 0.0
                p2 = float(rng.uniform(0, 2 * np.pi)) if ctrl == "randph" else 0.0
                syms = [
                    {"amp_mv": s1["amp_mv"] * scale, "freq_mhz": s1["freq_mhz"] + f_shift, "start_ns": 0.0, "dur_ns": 1000.0, "phase": p1},
                    {"amp_mv": s2["amp_mv"] * scale, "freq_mhz": s2["freq_mhz"] + f_shift, "start_ns": 500.0, "dur_ns": 500.0, "phase": p2},
                ]
                wf = make_pulse_waveform(syms, gap_ns=2000.0)
                meta = {"block": "E4", "control": ctrl, "order": name, "rep": r}
                trials.append((f"E4_{ctrl}_{name}_r{r}", wf, meta))
    return trials


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=str, default="data_e0_e4.h5")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--blocks", nargs="+", default=["E1", "E2", "E3", "E4"])
    parser.add_argument("--pilot-every", type=int, default=10)
    parser.add_argument("--n-rep-scale", type=float, default=1.0, help="Scale repeats (e.g. 0.1 for quick smoke test)")
    args = parser.parse_args()

    rig = HardwareRig(args.out, dry_run=args.dry_run, script_filename=__file__)
    pilot_wf = build_pilot_wf()

    all_trials = []
    if "E1" in args.blocks:
        all_trials.extend(build_e1_trials(n_rep=max(1, int(round(30 * args.n_rep_scale)))))
    if "E2" in args.blocks:
        all_trials.extend(build_e2_trials(n_rep=max(1, int(round(80 * args.n_rep_scale)))))
    if "E3" in args.blocks:
        all_trials.extend(build_e3_trials(n_rep=max(1, int(round(40 * args.n_rep_scale)))))
    if "E4" in args.blocks:
        all_trials.extend(build_e4_trials(n_rep=max(1, int(round(40 * args.n_rep_scale)))))

    rng = np.random.default_rng(2026)
    rng.shuffle(all_trials)

    manifest = []
    t0 = time.time()
    print(f"[START] {len(all_trials)} randomized trials | pilot every {args.pilot_every} | dry_run={args.dry_run}")

    try:
        for idx, (key, wf, meta) in enumerate(all_trials):
            if idx % args.pilot_every == 0:
                p_idx = idx // args.pilot_every
                rig.measure_trial(f"E0_PILOT_{p_idx}", pilot_wf, {"block": "E0", "pilot_idx": p_idx})
                if p_idx % 5 == 0:
                    rig.configure_sa_excitation()
                    rig.measure_trial(f"E0_EX_{p_idx}", pilot_wf, {"block": "E0_EX", "pilot_idx": p_idx}, measure_baseline=False)
                    rig.configure_sa_sample()

            rig.measure_trial(key, wf, meta)
            manifest.append({"key": key, "t_s": round(time.time() - t0, 2), **meta})
            if (idx + 1) % 100 == 0:
                print(f"  Progress: {idx + 1}/{len(all_trials)} ({time.time() - t0:.1f}s)")
    except KeyboardInterrupt:
        print("[INTERRUPTED] Flushing data safely...")
    finally:
        rig.close()
        Path(args.out + ".manifest.json").write_text(json.dumps(manifest, indent=2))
        print(f"[DONE] Saved {len(manifest)} trials + manifest to {args.out}")


if __name__ == "__main__":
    main()
