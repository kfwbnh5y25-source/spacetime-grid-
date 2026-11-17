"""
meas3_depth_repeat.py - REPEAT-AVERAGED DEPTH AT THE MONEY OPERATING POINTS (Measurement #3).

GOAL: pin the phase-gating DEPTH (the >=50% money number) with error bars, and kill
the baseline-recalibration artifact. In the old 90-s-cached-baseline scripts, a mid-block
baseline step fakes large |c1| in ~9% of (delay, f1) phase blocks; here every trial gets
a FRESH baseline and repeats are interleaved.

Operating points (both CW-overlap, split-phase +-phi/2 convention):
    OP 0 (money): f1 = 4069, a1 = 324, f2 = 4111, a2 = 271 mV  (the +20.5 dB point)
    OP 1:         f1 = 4070, a1 = 327, f2 = 4135, a2 = 430 mV  (fig3/a_old strong zone)
60 phases x 8 trials per OP, fresh baseline EVERY trial. Estimated ~25-35 min total.

Usage: python meas3_depth_repeat.py [--dry-run] [--out PATH]
"""
import argparse
import os
import sys
import time

import numpy as np

# chassis is local to this folder: copy the whole measure/ directory to the setup PC
from hw_core_spectro import SAConfig, SpectroRig, phase_grid, two_pulse_overlap  # noqa: E402

OPS = [
    dict(name="MONEY", f1=4069.0, a1=324.0, f2=4111.0, a2=271.0),
    dict(name="STRONG", f1=4070.0, a1=327.0, f2=4135.0, a2=430.0),
]
N_PHI, TRIALS = 60, 8
RAMP_NS, OVERLAP_NS = 2000.0, 4000.0
F_MIN_MHZ, F_MAX_MHZ = 3500.0, 4500.0       # boss requirement: no out-of-band tones


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--out", default="meas3_depth_repeat.h5")
    ap.add_argument("--trials", type=int, default=TRIALS)
    args = ap.parse_args()

    for op in OPS:
        assert F_MIN_MHZ <= op["f1"] <= F_MAX_MHZ and F_MIN_MHZ <= op["f2"] <= F_MAX_MHZ, \
            f"tone frequency out of band: {op['f1']}/{op['f2']}"
        assert op["a1"] + op["a2"] <= 760.0, \
            f"OP {op['name']}: overlap peak {(op['a1'] + op['a2'])/1e3:.2f} V exceeds amplifier-compression ceiling"
    rig = SpectroRig(args.out, dry_run=args.dry_run,
                     sa_config=SAConfig(span_hz=5e6, rbw_hz=200e3, avg_count=4,
                                        sweep_time=0.25, ref_level_dbm=-38.8,
                                        scale_div=4.0))
    t0 = time.time()
    n_planned = len(OPS) * N_PHI * args.trials
    done = 0
    try:
        for op in OPS:
            rig.cfg.center_hz = op["f2"] / 2 * 1e6
            rig.configure_sa()
            phases = phase_grid(N_PHI, shuffle=True, seed=hash(op["name"]) % 2**31)
            order = list(range(N_PHI))
            import random
            random.seed(len(op["name"]))
            for t in range(args.trials):
                random.shuffle(order)  # interleave phases across trials
                for i in order:
                    ph = phases[i]
                    wf, spec = two_pulse_overlap(op["a1"], op["f1"], op["a2"], op["f2"], ph,
                                                 t0_ns=RAMP_NS, overlap_ns=OVERLAP_NS,
                                                 ramp_ns=RAMP_NS, split_phase=True)
                    key = f"DEP_{op['name']}_p{ph:.4f}_t{t}"
                    rig.register_waveform(key, wf, spec)
                    rig.measure(key, wf, meta={**{k: v for k, v in op.items() if k != "name"},
                                               "op": op["name"], "trial": t},
                                with_baseline=True)  # ALWAYS fresh baseline
                    done += 1
                    if done % 40 == 0:
                        el = time.time() - t0
                        print(f"[{done}/{n_planned}] {el/60:.1f} min, ETA "
                              f"{el/done*(n_planned-done)/60:.0f} min", flush=True)
    finally:
        rig.close()
        print(f"DONE {done} trials in {(time.time()-t0)/60:.1f} min -> {args.out}")


if __name__ == "__main__":
    main()
