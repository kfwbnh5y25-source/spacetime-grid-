"""
meas1_longgap_memory.py - LONG-GAP DELAY SCAN AT LOW DETUNING (Measurement #1).

GOAL: settle the memory-time tension. Two-tone delay scans stopped at +15 ns; at
Delta = 50-58 MHz the 16-phase (delay_delta.h5) contrast was still alive at 7 ns,
and the small-signal theory (GAMMA_EFF_MATH.md Eq. 10) predicts a memory time
tau_mem = 1/(2*Gamma_k) = 51-61 ns. This scan extends t_gap to 200 ns.

Configs (a1 = 327 mV, a2 = 588 mV, pulse geometry copied from delay.h5:
tone1 tukey 300 samples = 4.6875 ns, tone2 tukey 400 samples = 6.25 ns):
    CFG 0: f1 = 4081, f2 = 4131  (Delta = 50 MHz, 60 phases)
    CFG 1: f1 = 4077, f2 = 4135  (Delta = 58 MHz, 60 phases)
    CFG 2: f1 = 4074, f2 = 4138  (Delta = 64 MHz, 16 phases)   <- NEW
    CFG 3: f1 = 4070, f2 = 4140  (Delta = 70 MHz, 16 phases)   <- NEW
Configs 2/3 (the 60-phase delay-recheck2 found the tau transition at
Delta = 62-70 MHz: no decay at Delta <= 61, tau ~ 1.5-4 ns at 62-70) extend
the long-gap scan across the transition - the t_gap > 15 ns behavior there is
UNKNOWN and is Fig 4's decay-with-delay axis. 16-phase grids keep c1 SNR
adequate (sigma_A1 ~ 2-3 uW vs plateau 50-110 uW) and bound the runtime.

GEOMETRY: timeline shifted so every start >= 0 while preserving the delay.h5
relative timing exactly (for t_gap >= 0 it is bit-identical to delay.h5; for
t_gap < 0 both pulses are shifted together). Phase convention: both tones
referenced to FIXED t0_ref = 0, phase phi on tone2 only -> the tone2 phase at its
own turn-on winds by 2*pi*Delta*t_gap across the scan (phi_eff IS in the data;
do not re-reference when analyzing |c1|).

DRIVER NOTES (verified against Keysight_AWG_M8195A.py):
- All waveforms share ONE fixed loop length (multiple of 1280 samples) so the
  segment is defined once, not re-downloaded per delay (send_waveform deletes +
  redefines the segment on length change).
- Peak amplitude: non-overlap = 588 mV; full-overlap worst case = |327|+|588| =
  915 mV (identical to the old delay.h5 data; loops stay in [-1, 1] as the
  driver requires. Use --peak-mv 760 to enforce the amplifier-compression
  ceiling for FRESH physics; 915 mV peaks are kept ONLY for continuity with
  the old scans).

Estimated runtime: ~3-4 h (full); `--trials 1` halves; `--config 0/1` one detuning.
Usage: python meas1_longgap_memory.py [--dry-run] [--out PATH] [--trials N] [--config 0|1]
"""
import argparse
import json
import os
import sys
import time

import numpy as np

# chassis is local to this folder: copy the whole measure/ directory to the setup PC
from hw_core_spectro import (  # noqa: E402
    AWG_BLOCK_SAMPLES, LoopSpec, Pulse, SAConfig, SpectroRig, build_loop,
    phase_grid,
)

F1_F2_NPHI = [(4081.0, 4131.0, 60), (4077.0, 4135.0, 60),
              (4074.0, 4138.0, 16), (4070.0, 4140.0, 16)]
A1_MV, A2_MV = 327.0, 588.0
D1_NS, D2_NS = 300 / 64.0, 400 / 64.0      # 4.6875 ns, 6.25 ns (delay.h5 geometry)
F_MIN_MHZ, F_MAX_MHZ = 3500.0, 4500.0       # boss requirement: no out-of-band tones

DELAYS_NS = ([-8.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.0]  # overlap refs
             + list(np.arange(0.0, 15.0, 1.0))                             # fine near onset
             + list(np.arange(15.0, 60.0, 2.5))                            # decay range
             + list(np.arange(60.0, 201.0, 10.0)))                         # long tail
N_PHI = 60
TRIALS = 1   # 60(16)-phase grids suffice: sigma_A1 ~ 1-3 uW vs plateau 50-110 uW
QUIET_TAIL_NS = 100.0
T_LOOP_NS = 2000.0  # June delay2.h5 loop - fixed for ALL delays (one segment).
                    # 2 us loop = 0.5 MHz pump comb spacing (comb-quantization
                    # safe; the 140 ns sanity loop's 7.14 MHz comb went flat).


T_START_NS = 1000.0  # June timeline: tone1 ends at t=1000 ns in the 2 us loop


def build_pulses(cfg, t_gap_ns, phi):
    f1, f2, _n = cfg
    # June delay2.h5 timeline: tone1 ends at T_START, tone2 starts T_START+t_gap
    s2_raw = T_START_NS + t_gap_ns          # can dip below tone1 start for t_gap<0
    shift = max(0.0, -(s2_raw - D1_NS))     # common shift keeps all starts >= 0
    s1 = T_START_NS - D1_NS + shift
    s2 = s2_raw + shift
    p1 = Pulse(A1_MV, f1, start_ns=s1, dur_ns=D1_NS, phase=0.0)
    p2 = Pulse(A2_MV, f2, start_ns=s2, dur_ns=D2_NS, phase=phi)
    end = max(s1 + D1_NS, s2 + D2_NS)
    return [p1, p2], end, f2


def waveforms_for(cfg, t_gap_ns, phi, peak_ceiling_mv):
    import hw_core_spectro as hc
    pulses, end, f2 = build_pulses(cfg, t_gap_ns, phi)
    spec = LoopSpec(pulses=pulses, t0_ref_ns=0.0, gap_ns=T_LOOP_NS - end)
    old, hc.MAX_PEAK_MV = hc.MAX_PEAK_MV, peak_ceiling_mv
    wf = build_loop(spec, edge_ns=1.0)     # boss-verified: June edges = 1 ns tukey
    hc.MAX_PEAK_MV = old
    # phi_eff at tone2 turn-on (the constant 2*pi*f2*T_START offset is common
    # to all delays -> only the slope 2*pi*Delta*t_gap matters for winding)
    phi2_at_on = 2 * np.pi * f2 * 1e6 * ((T_START_NS + max(t_gap_ns, 0.0)) * 1e-9) + phi
    return wf, spec, float(np.angle(np.exp(1j * phi2_at_on)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--out", default="meas1_longgap_memory.h5")
    ap.add_argument("--trials", type=int, default=TRIALS)
    ap.add_argument("--config", type=int, default=None, choices=(0, 1, 2, 3))
    ap.add_argument("--peak-mv", type=float, default=916.0,
                    help="peak-amplitude ceiling in mV. Default 916 mV = driver limit "
                         "(1000 mV) with margin, matching the old delay.h5 scans whose "
                         "full-overlap waveforms peak at 915 mV. Pass 760 to enforce the "
                         "amplifier-compression ceiling (this REJECTS the old-geometry "
                         "overlap waveforms).")
    args = ap.parse_args()

    for f1, f2, _ in F1_F2_NPHI:
        assert F_MIN_MHZ <= f1 <= F_MAX_MHZ and F_MIN_MHZ <= f2 <= F_MAX_MHZ, \
            f"tone frequency out of band: {f1}/{f2}"

    cfgs = F1_F2_NPHI if args.config is None else [F1_F2_NPHI[args.config]]
    sa = SAConfig(center_hz=cfgs[0][1] / 2 * 1e6, span_hz=30e6, points=1001,
                  rbw_hz=1e6, avg_count=2, sweep_time=0.5, ref_level_dbm=-53.8,
                  scale_div=0.5)
    rig = SpectroRig(args.out, dry_run=args.dry_run, sa_config=sa)
    t_start = time.time()
    n_planned = len(cfgs) * len(DELAYS_NS) * N_PHI * args.trials
    done = 0
    try:
        for ci, cfg in enumerate(cfgs):
            f1, f2, n_phi = cfg
            rig.cfg.center_hz = f2 / 2 * 1e6
            rig.configure_sa()
            phases = phase_grid(n_phi, shuffle=True, seed=100 + ci)
            delays = list(DELAYS_NS)
            import random
            random.seed(7 + ci)
            random.shuffle(delays)
            # tone2-only references at a few delays (idler baseline without seeding)
            for d in delays[::10]:
                wf, spec, phi2 = waveforms_for(cfg, d, 0.0, args.peak_mv)
                spec.pulses[0].amp_mv = 0.0
                wf = build_loop(spec, edge_ns=1.0)
                rig.register_waveform(f"REFTO2_c{ci}_d{d:.2f}", wf, spec)
                rig.measure(f"REFTO2_c{ci}_d{d:.2f}", wf,
                            meta={"delay_ns": d, "phi_eff_turn_on": phi2, "ref": "tone2only"})
                done += 1
            for t in range(args.trials):
                for d in delays:
                    for ph in phases:
                        wf, spec, phi2 = waveforms_for(cfg, d, ph, args.peak_mv)
                        key = f"LGM_c{ci}_d{d:+.3f}_p{ph:.4f}_t{t}"
                        rig.register_waveform(key, wf, spec)
                        rig.measure(key, wf, meta={
                            "t_gap_ns": float(d), "Delta_mhz": f2 - f1,
                            "phi_eff_turn_on": phi2,
                            "phi_programmed": float(ph), "trial": t})
                        done += 1
                        if done % 60 == 0:
                            el = time.time() - t_start
                            print(f"[{done}/{n_planned}] {el/60:.1f} min elapsed, "
                                  f"ETA {el/done*(n_planned-done)/60:.0f} min", flush=True)
    finally:
        rig.close()
        print(f"DONE {done} trials in {(time.time()-t_start)/60:.1f} min -> {args.out}")


if __name__ == "__main__":
    main()
