"""
meas4_duty.py - LOOP-PERIOD / DUTY CONTROL AT A FIXED OP (~15 min).

Purpose: discriminate between two candidate explanations for the flat burst
spectra in the first sanity run, both at the SAME operating point as today's
probe (f1 = 4075, f2 = 4149 MHz, a1 = 327 mV, a2 = 588 mV, SA center 2074.5):

  D0_june    : delay2.h5 geometry EXACTLY - single 4.7 + 6.25 ns burst pair
               in a 2000 ns loop (comb spacing 0.5 MHz, duty 0.55%),
               June tukey alpha = 64/300 (edge 0.5 ns).
  D1_burst140: today's failing geometry - same pulses, 140 ns loop
               (comb spacing 7.14 MHz, duty 8%).
  D3_CW      : f1_4069-style ramp 2 us + overlap 4 us + gap 2 us
               (duty 75%; already proven today: 322 uW line).

PREDICTION (comb-quantization mechanism):
  D0_june LINE PRESENT even at 0.55% duty   [else: pure-duty mechanism]
  D1_burst140 FLAT                          [today's observed result]
  D3_CW LINE PRESENT                        [today's observed result]

Decision rule:
  - D0 present, D1 flat -> loop period is the knob; overnight long-gap scan
    (meas1) keeps the 2000 ns June loop.
  - D0 flat, D3 present -> pure-duty story; meas1 uses a duty-carrying train.
16 phases x 2 trials per variant, interleaved, live verdict per trace.
"""
import argparse
import itertools
import random
import sys
import time

import numpy as np

sys.path.insert(0, "")
import hw_core_spectro as hc
from hw_core_spectro import (LoopSpec, Pulse, SAConfig, SpectroRig,
                             two_pulse_overlap, phase_grid)

F1, F2 = 4075.0, 4149.0
A1M, A2M = 327.0, 588.0
D1_JUNE_NS, D2_JUNE_NS = 300.0, 400.0    # boss-verified June burst lengths (ns)
D1_SHORT_NS, D2_SHORT_NS = 300 / 64.0, 400 / 64.0   # probe-era short bursts (control)
F_MIN, F_MAX = 3500.0, 4500.0
assert F_MIN <= F1 <= F_MAX and F_MIN <= F2 <= F_MAX
CENTER_HZ = F2 / 2 * 1e6
LINE_WINDOW_MHZ = 1.5
N_PHI = 16
TRIALS = 2
PHASES = phase_grid(N_PHI)
DUTY_STR = {"D0_june": "35%", "D1_burst140": "8%", "D3_CW": "75%"}


def wf_variant(variant, phi, peak_mv):
    """Returns (wf, spec, kind) or (None, None, reason)."""
    old, hc.MAX_PEAK_MV = hc.MAX_PEAK_MV, peak_mv
    try:
        if variant == "D0_june":
            spec = LoopSpec(pulses=[Pulse(A1M, F1, 1000.0 - D1_JUNE_NS, D1_JUNE_NS, 0.0),
                                    Pulse(A2M, F2, 1000.0, D2_JUNE_NS, phi)],
                            t0_ref_ns=0.0, gap_ns=2000.0 - 1000.0 - D2_JUNE_NS)
            wf = hc.build_loop(spec, edge_ns=1.0)
            return wf, spec, "ok"
        if variant == "D1_burst140":
            spec = LoopSpec(pulses=[Pulse(A1M, F1, 0.0, D1_SHORT_NS, 0.0),
                                    Pulse(A2M, F2, D1_SHORT_NS, D2_SHORT_NS, phi)],
                            t0_ref_ns=0.0, gap_ns=130.0)
            wf = hc.build_loop(spec, edge_ns=1.0)
            return wf, spec, "ok"
        if variant == "D3_CW":
            wf, spec = two_pulse_overlap(A1M, F1, A2M, F2, phi,
                                         t0_ns=1000.0, overlap_ns=4000.0,
                                         ramp_ns=2000.0, gap_ns=2000.0)
            wf = hc.build_loop(spec, edge_ns=1.0)
            return wf, spec, "ok"
        return None, None, f"unknown variant {variant}"
    except ValueError as e:
        return None, None, str(e)
    finally:
        hc.MAX_PEAK_MV = old


def line_power(f_ax, sig, bl):
    d = sig - bl
    win = np.abs(f_ax / 1e6 - CENTER_HZ / 1e6) <= LINE_WINDOW_MHZ
    return float(d[win].sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--out", default="meas4_duty.h5")
    ap.add_argument("--peak-mv", type=float, default=916.0,
                    help="916 keeps overlap peaks (worst 915 mV); 760 = strict "
                         "compression ceiling (D3 may skip)")
    ap.add_argument("--trials", type=int, default=TRIALS)
    args = ap.parse_args()

    rig = SpectroRig(args.out, dry_run=args.dry_run,
                     sa_config=SAConfig(center_hz=CENTER_HZ, span_hz=30e6,
                                        points=1001, rbw_hz=1e6, avg_count=2,
                                        sweep_time=0.5, ref_level_dbm=-53.8,
                                        scale_div=0.5))
    t_start = time.time()
    data = {}
    skipped = 0
    try:
        plan = [(v, ph, t) for v, ph, t in
                itertools.product(("D0_june", "D1_burst140", "D3_CW"),
                                  PHASES, range(args.trials))]
        random.seed(11)
        random.shuffle(plan)
        print(f"meas4: {len(plan)} planned measures, ETA ~= "
              f"{len(plan) * 1.5 / 60 * (0.5 if args.dry_run else 1):.0f} min")
        done = 0
        for v, ph, t in plan:
            wf, spec, note = wf_variant(v, ph, args.peak_mv)
            if wf is None:
                skipped += 1
                print(f"  [skip] {v} phi={ph:.3f}: {note}")
                continue
            key = f"DUTY_{v}_p{ph:.4f}_t{t}"
            rig.register_waveform(key, wf, spec)
            sig, bl, fr = rig.measure(key, wf, meta={
                "variant": v, "phi": float(ph), "trial": t})
            p = line_power(fr, sig, bl)
            data.setdefault(v, []).append((float(ph), p))
            done += 1
            print(f"  [{v} t{t}] phi={ph:5.2f}  line={p*1e6:9.1f} uW  "
                  f"({done}/{len(plan)-skipped}, {(time.time()-t_start)/60:.1f} min)",
                  flush=True)
    finally:
        rig.close()

    print(f"\nmeas4 summary ({(time.time()-t_start)/60:.1f} min)\n" + "-" * 66)
    for v in ("D0_june", "D1_burst140", "D3_CW"):
        rows = data.get(v, [])
        if len(rows) < 8:
            print(f"{v}: insufficient data ({len(rows)} points)")
            continue
        ph = np.array([r[0] for r in rows])
        pw = np.array([r[1] for r in rows])
        c1 = np.mean(pw * np.exp(-1j * ph))
        A1, A0 = float(2 * np.abs(c1)), float(np.mean(pw))
        m = A1 / max(A0, 1e-30)
        depth = 2 * A1 / (A0 + A1)
        verdict = "LINE PRESENT" if A0 > 30e-6 else ("weak" if A0 > 3e-6 else "FLAT")
        print(f"{v:12s} duty={DUTY_STR[v]:>5s}  A0={A0*1e6:8.1f} uW  "
              f"A1={A1*1e6:8.1f} uW  A1/A0={m:.2f}  mod-depth={depth*100:5.1f}%  "
              f"-> {verdict}")
    print("\npredictions (comb mechanism): D0 PRESENT, D1 FLAT, D3 PRESENT")


if __name__ == "__main__":
    main()
