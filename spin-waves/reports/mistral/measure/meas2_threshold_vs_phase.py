"""
meas2_threshold_vs_phase.py - THRESHOLD vs RELATIVE PHASE (Measurement #2).

GOAL: direct test of the two-tone threshold formula (GAMMA_EFF_MATH.md Eq. 7):
    h_th(Psi) = [ sqrt(Gamma^2 + delta^2 - W^2 sin^2(Psi)) - W cos(Psi) ] / V
by sweeping a2 across threshold at several programmed relative phases. Because W
scales with a1^2 (fixed here), the threshold shift should be linear in cos(Psi)
with slope 2W/V-ish; the fit gives BOTH W and V*Gamma product from spectroscopy.

Operating point (standard f1_4069_5-like): f1 = 4069 MHz, a1 = 324 mV fixed;
f2 = 4111 MHz; a2 swept 268 -> 312 mV (2 mV steps, 23 points) across the ~271-300 mV
threshold region. 16 relative phases on the uniform grid (split-half convention
+phi/2 on tone1, -phi/2 on tone2 -- f1_4069_5 convention).

Estimated runtime: ~25 min (23 a2 x 16 phases x 2 trials, sweep 0.25 s + baseline).

Usage: python meas2_threshold_vs_phase.py [--dry-run] [--out PATH]
"""
import argparse
import os
import sys
import time

import numpy as np

# chassis is local to this folder: copy the whole measure/ directory to the setup PC
from hw_core_spectro import (  # noqa: E402
    SAConfig, SpectroRig, phase_grid, two_pulse_overlap,
)

F1_MHZ, F2_MHZ = 4069.0, 4111.0
A1_MV = 324.0
F_MIN_MHZ, F_MAX_MHZ = 3500.0, 4500.0       # boss requirement: no out-of-band tones
A2_GRID = list(np.arange(268.0, 314.0, 2.0))   # 23 points straddling threshold
N_PHI = 16
TRIALS = 2
RAMP_NS, OVERLAP_NS = 2000.0, 4000.0            # long CW-overlap burst (fig3 style)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--out", default="meas2_threshold_vs_phase.h5")
    ap.add_argument("--trials", type=int, default=TRIALS)
    args = ap.parse_args()

    assert F_MIN_MHZ <= F1_MHZ <= F_MAX_MHZ and F_MIN_MHZ <= F2_MHZ <= F_MAX_MHZ, \
        f"tone frequency out of band: {F1_MHZ}/{F2_MHZ}"
    sa = SAConfig(center_hz=F2_MHZ / 2 * 1e6, span_hz=5e6, points=1001,
                  rbw_hz=200e3, avg_count=2, sweep_time=0.25, ref_level_dbm=-38.8,
                  scale_div=4.0)
    rig = SpectroRig(args.out, dry_run=args.dry_run, sa_config=sa)
    phases = phase_grid(N_PHI, shuffle=True, seed=42)
    a2s = list(A2_GRID)
    import random
    random.seed(11)
    random.shuffle(a2s)

    t0 = time.time()
    n_planned = len(a2s) * N_PHI * args.trials + 16
    done = 0
    try:
        # tone1-only reference sweep (threshold denominator / leakage check)
        for ph in phases[:16]:
            wf, spec = two_pulse_overlap(A1_MV, F1_MHZ, 0.0, F2_MHZ, ph,
                                         t0_ns=RAMP_NS, overlap_ns=OVERLAP_NS,
                                         ramp_ns=RAMP_NS, split_phase=True)
            key = f"REFT1_a2_0_p{ph:.4f}"
            rig.register_waveform(key, wf, spec)
            rig.measure(key, wf, meta={"a2_mv": 0.0, "ref": "tone1only"})
            done += 1
        for t in range(args.trials):
            for a2 in a2s:
                for ph in phases:
                    wf, spec = two_pulse_overlap(A1_MV, F1_MHZ, a2, F2_MHZ, ph,
                                                 t0_ns=RAMP_NS, overlap_ns=OVERLAP_NS,
                                                 ramp_ns=RAMP_NS, split_phase=True)
                    key = f"TVP_a2_{a2:.1f}_p{ph:.4f}_t{t}"
                    rig.register_waveform(key, wf, spec)
                    rig.measure(key, wf, meta={
                        "a2_mv": float(a2), "a1_mv": A1_MV, "Delta_mhz": F2_MHZ - F1_MHZ,
                        "psi_p_rel_t0": float(ph), "trial": t})
                    done += 1
                    if done % 50 == 0:
                        el = time.time() - t0
                        print(f"[{done}/{n_planned}] {el/60:.1f} min, ETA "
                              f"{el/done*(n_planned-done)/60:.0f} min", flush=True)
    finally:
        rig.close()
        print(f"DONE {done} trials in {(time.time()-t0)/60:.1f} min -> {args.out}")


if __name__ == "__main__":
    main()
