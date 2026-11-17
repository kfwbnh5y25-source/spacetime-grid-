"""
meas0_sanity.py - SANITY CHECK (~5 min): chain + right-physics verification.

WHAT "SANITY" MEANS HERE (assessed):
  S1 chain integrity   the AWG waveform reaches the SA on the excitation path:
                       both carriers SEEN with powers consistent with the
                       probe reference (f1 -59.1 dBm, f2 -43.3 dBm today).
  S2 right quantity    on the sample path, at the June-geometry operating
                       point (tone1 4075/327 mV 300 ns, t_gap=0, tone2
                       4149/588 mV 400 ns, 2 us loop, 32 ns tukey edges --
                       EXACTLY the delay2.h5 protocol), the idler line
                       appears AT f2/2 = 2074.5 MHz with A0 in the June
                       anchor range (267 uW band-summed ±1.5 MHz).
  S3 right physics     the phase response carried by that line: modulation
                       A1/A0 ~ 0.42 (anchor), memory decay A1(8ns) << A1(0.5ns)
                       (anchor ratio 0.12, tau 2.13 ns), and the phi_eff
                       winding dphi0/dt_gap ~ 20.7 deg/ns (77.7% of 2*pi*Delta).
If S1-S3 hold, the rig measures what the paper claims it measures; the long
campaign then only extends the (t_gap, Delta) coverage.

Blocks:
  C  excitation-path readout (switch 1,1) of the burst+CW waveforms;
     SA set EXPLICITLY (center (F1+F2)/2, span (F2-F1)+40 MHz, RBW 2 MHz),
     back to sample with every setting re-asserted.
  R  tone2-only reference x2 (line floor with the pump chain live).
  S  t_gap = 0 June burst, 16 phases x 2 trials (S2/S3 at the money point).
  T  delay mini-scan  t_gap in {-8,-2,-0.5,0,1,2,3,5,8} ns x 16 phases
     (S3 decay+winding).
Live verdict printed per trace from the FIRST trace on; final table + GO/NO-GO.

Usage: python meas0_sanity.py [--dry-run] [--out PATH] [--peak-mv 916] [--selftest]
"""
import argparse
import sys
import time

import numpy as np
from scipy.signal import medfilt

sys.path.insert(0, "")
import hw_core_spectro as hc
from hw_core_spectro import LoopSpec, Pulse, SAConfig, SpectroRig

F1, F2 = 4075.0, 4149.0
A1M, A2M = 327.0, 588.0
D1_NS, D2_NS = 300.0, 400.0          # June pulse lengths in NANOSECONDS
EDGE_NS = 1.0                        # boss-verified: 64-sample (1 ns) tukey ramps
F_MIN, F_MAX = 3500.0, 4500.0
assert F_MIN <= F1 <= F_MAX and F_MIN <= F2 <= F_MAX

T_START = 1000.0                     # tone1 ends at t=1000 ns in the 2 us loop
T_LOOP = 2000.0
CENTER_HZ = F2 / 2 * 1e6
LINE_WINDOW_MHZ = 1.5
SPAN_HZ = 30e6

PHASES_T = np.arange(16) * 2 * np.pi / 16
DELAYS_T = [-8.0, -2.0, -0.5, 0.0, 1.0, 2.0, 3.0, 5.0, 8.0]

# anchors from the delay2.h5 recheck (June, Delta = 74 MHz):
ANCH = {"A0_uW": 267.0, "A1_uW": 111.0, "A1_ratio_8ns": 0.12,
        "slew_deg_ns": 20.7, "tau_ns": 2.13}
# probe reference (today): excitation-path peaks
ANCH_EX = {"dbm_f1": -59.1, "dbm_f2": -43.3}


def wf_june(phi, t_gap, a1=A1M, a2=A2M, dur2=D2_NS, edge_ns=EDGE_NS):
    """June delay2.h5 geometry, single burst pair, 2 us fixed loop."""
    s2_raw = T_START + t_gap
    shift = max(0.0, -(s2_raw - D1_NS))
    s1 = T_START - D1_NS + shift
    s2 = s2_raw + shift
    spec = LoopSpec(pulses=[Pulse(a1, F1, s1, D1_NS, 0.0),
                            Pulse(a2, F2, s2, dur2, phi)],
                    t0_ref_ns=0.0, gap_ns=T_LOOP - (s2 + dur2))
    old, hc.MAX_PEAK_MV = hc.MAX_PEAK_MV, 916.0
    try:
        wf = hc.build_loop(spec, edge_ns=edge_ns)
    except ValueError as e:
        print(f"  [skip] {e}")
        return None, None
    finally:
        hc.MAX_PEAK_MV = old
    return wf, spec


def line_power(sig, bl, freqs):
    d = medfilt(sig - bl, 9)
    win = np.abs(freqs / 1e6 - CENTER_HZ / 1e6) <= LINE_WINDOW_MHZ
    noise_bins = ~win & (d < np.quantile(d, 0.9))
    return float(d[win].sum()), float(np.std(d[noise_bins]))


def c12(points):
    if len(points) < 8:
        return None
    ph = np.array([p[0] for p in points])
    pw = np.array([p[1] for p in points])
    c1 = np.mean(pw * np.exp(-1j * ph))
    return float(2 * np.abs(c1)), float(np.mean(pw)), float(np.angle(c1))


# dry-run-only physical model (the delay2 anchors), so --dry-run is an
# end-to-end test that MUST print GO: P(phi,t_gap) = (A0+A1cos(phi+slew*t_gap))
# * exp(-max(t_gap,0)/tau)
A0_MDL, A1_MDL = 267e-6, 111e-6
TAU_MDL = 2.13
SLEW_MDL = 2 * np.pi * (F2 - F1) * 1e-3 * 0.777


def model_line(freqs, phi, t_gap):
    win = np.abs(freqs / 1e6 - CENTER_HZ / 1e6) <= 0.5
    out = np.zeros_like(freqs)
    out[win] = (A0_MDL + A1_MDL * np.cos(phi + SLEW_MDL * t_gap)) \
        * np.exp(-max(t_gap, 0.0) / TAU_MDL) / win.sum()
    return out


_RESULTS = []


def check(tag, ok, msg):
    print(f"[{'PASS' if ok else 'FAIL':4s}] {tag}: {msg}")
    _RESULTS.append((tag, ok))


def selftest():
    """Analysis pipeline self-test on synthetic spectra (no hardware)."""
    freqs = np.linspace(CENTER_HZ - SPAN_HZ / 2, CENTER_HZ + SPAN_HZ / 2, 1001)
    rng = np.random.default_rng(7)
    A0_t, A1_t, tau_t, slew_t = 267e-6, 111e-6, 2.13, 0.3614  # rad/ns = 20.7 deg
    def synth(A1v, phi):
        sig = np.abs(rng.normal(5e-8, 1.5e-8, len(freqs)))
        win = np.abs(freqs / 1e6 - CENTER_HZ / 1e6) <= 0.5
        sig[win] += (A0_t + A1v * np.cos(phi + slew_t * 0.0)) / win.sum()
        return sig, np.abs(rng.normal(5e-8, 1.5e-8, len(freqs)))
    rows = [c1 for c1 in
            ((phi, line_power(*synth(A1_t * 1, phi), freqs)) for phi in
             np.arange(0, 2 * np.pi, np.pi / 15))]
    rows = [(phi, p) for phi, (p, _) in rows]
    A1r, A0r, _ = c12(rows)
    check("selftest A0", abs(A0r / A0_t - 1) < 0.1, f"{A0r*1e6:.0f} vs {A0_t*1e6:.0f} uW")
    check("selftest A1", abs(A1r / A1_t - 1) < 0.15, f"{A1r*1e6:.0f} vs {A1_t*1e6:.0f} uW")
    dec = {}
    for t in (0.0, 8.0):
        rr = []
        for phi in np.arange(0, 2 * np.pi, np.pi / 15):
            sig, bl = synth(A1_t * np.exp(-t / tau_t), phi)
            p, _ = line_power(sig, bl, freqs)
            rr.append((phi, p))
        dec[t] = c12(rr)[0]
    check("selftest decay", abs(dec[8.0] / dec[0.0] - np.exp(-8 / tau_t)) < 0.05,
          f"ratio {dec[8.0]/dec[0.0]:.3f} vs {np.exp(-8/tau_t):.3f}")
    sys.exit(0 if all(ok for _, ok in _RESULTS) else 1)


def live(tag, phi, p, n, ref):
    line_uw, n_uw = p * 1e6, n * 1e6
    if ref[0] is None:
        print(f"  [{tag}] phi={phi:5.2f} line={line_uw:8.1f} uW  noise={n_uw:5.2f} uW")
        return
    x = p / max(ref[0], 1e-30)
    print(f"  [{tag}] phi={phi:5.2f} line={line_uw:8.1f} uW  noise={n_uw:5.2f} uW  "
          f"x_ref={x:6.1f}  {'SEED:OK' if x >= 3 else 'SEED:LOW'}")


def block_verdict(tag, rows, ref):
    cc = c12(rows)
    if cc is None:
        print(f"  -- {tag}: {len(rows)} phases, block verdict pending")
        return None
    A1, A0, phi0 = cc
    m = A1 / max(A0, 1e-30)
    st = "OK " if (m >= 0.15 and A0 > 50e-6) else ("LOW" if m >= 0.08 else "NONE")
    print(f"  -- {tag}: A0={A0*1e6:7.1f} uW  A1={A1*1e6:7.1f} uW  A1/A0={m:.2f} "
          f"[{st}]  phi0={phi0:+.2f} rad")
    return cc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--out", default="meas0_sanity.h5")
    ap.add_argument("--peak-mv", type=float, default=916.0)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return

    rig = SpectroRig(args.out, dry_run=args.dry_run,
                     sa_config=SAConfig(center_hz=CENTER_HZ, span_hz=SPAN_HZ,
                                        points=1001, rbw_hz=1e6, avg_count=1,
                                        sweep_time=0.5, ref_level_dbm=-53.8,
                                        scale_div=0.5))
    rec = {"ref2": [], "two": {}, "exc": {}}
    ref_mean = [None]
    prev_cc, prev_t = [None], [0.0]
    t0 = time.time()
    try:
        # ---------- C: chain check (excitation path, explicit settings) ----------
        wf_burst, spec_burst = wf_june(0.0, 0.0)
        rig.register_waveform("SAN_EXC_burst", wf_burst, spec_burst)
        if not args.dry_run:
            rig.knet.switch_to_excitation()
            sa = rig.knet.sa
            sa.centerFrequency = (F1 + F2) / 2 * 1e6            # 4112 MHz
            sa.spanFrequency = (F2 - F1) * 1e6 + 40e6           # 114 MHz
            sa.resolutionBandwidthAutoState = "OFF"
            sa.resolutionBandwidth = 2e6
            sa.averageNumber = 2
            sa.sweepTime = 0.5
            rig.knet.awg.send_waveform(wf_burst)
            rig.knet.awg.start_generation()
            ef, edb = rig.knet.get_trace()
            rig.knet.awg.abort_generation()
            rig.knet.switch_to_sample()
            rig.configure_sa(rig.cfg)
            ef = np.asarray(ef, dtype=np.float64)
            elin = np.power(10.0, np.asarray(edb, dtype=np.float64) / 10.0)
            for f_tone, a_mv in ((F1, A1M), (F2, A2M)):
                win = np.abs(ef / 1e6 - f_tone) <= 2.0
                rec["exc"][f_tone] = (float(elin[win].max()), a_mv)
        else:
            rec["exc"] = {F1: (10 ** (-59.1 / 10), A1M), F2: (10 ** (-43.3 / 10), A2M)}

        # ---------- R: tone2-only reference ----------
        for t in range(2):
            wf, spec = wf_june(0.0, -D1_NS, a1=0.0)   # tone2-only (t_gap shift ok)
            if wf is None:
                continue
            key = f"SAN_REF2_t{t}"
            rig.register_waveform(key, wf, spec)
            sig, bl, fr = rig.measure(key, wf, meta={"ref": "tone2only", "trial": t})
            p, n = line_power(sig, bl, fr)
            rec["ref2"].append((p, n))
            ref_mean[0] = float(np.mean([x for x, _ in rec["ref2"]]))
            live(f"REF2 t{t}", 0.0, p, n, [None])

        # ---------- S: t_gap = 0, 16 phases x 2 ----------
        for t in range(2):
            for ph in PHASES_T:
                wf, spec = wf_june(float(ph), 0.0)
                if wf is None:
                    continue
                key = f"SAN_OVL_t{t}_p{ph:.4f}"
                rig.register_waveform(key, wf, spec)
                sig, bl, fr = rig.measure(key, wf, meta={"t_gap_ns": 0.0, "trial": t})
                if args.dry_run:
                    sig = sig + model_line(fr, float(ph), 0.0)
                p, n = line_power(sig, bl, fr)
                rec["two"].setdefault(0.0, []).append((float(ph), p, n))
                live(f"OVL  t{t}", float(ph), p, n, ref_mean)
            block_verdict(f"t_gap=+0.00 (t{t})", rec["two"][0.0], ref_mean)

        # ---------- T: delay mini-scan ----------
        for t_gap in DELAYS_T:
            if t_gap == 0.0:
                continue     # already measured (reused below)
            for ph in PHASES_T:
                wf, spec = wf_june(float(ph), t_gap)
                if wf is None:
                    continue
                key = f"SAN_d{t_gap:+.2f}_p{ph:.4f}"
                rig.register_waveform(key, wf, spec)
                sig, bl, fr = rig.measure(key, wf, meta={"t_gap_ns": t_gap})
                if args.dry_run:
                    sig = sig + model_line(fr, float(ph), t_gap)
                p, n = line_power(sig, bl, fr)
                rec["two"].setdefault(t_gap, []).append((float(ph), p, n))
                live(f"T={t_gap:+5.2f}", float(ph), p, n, ref_mean)
            cc = block_verdict(f"t_gap={t_gap:+.2f}", rec["two"][t_gap], ref_mean)
            if cc and prev_cc[0]:
                dph = (cc[2] - prev_cc[0][2] + np.pi) % (2 * np.pi) - np.pi
                slew = dph / (t_gap - prev_t[0]) * 180 / np.pi
                print(f"     live: A1 {prev_cc[0][0]*1e6:.0f}->{cc[0]*1e6:.0f} uW, "
                      f"slew {slew:+.1f} deg/ns")
            if cc:
                prev_cc[0], prev_t[0] = cc, t_gap
    finally:
        rig.close()

    # ---------- final table ----------
    print(f"\nfinal ({time.time()-t0:.0f} s)\n" + "-" * 66)

    # S1 chain
    if len(rec["exc"]) == 2:
        (p1, _), (p2, _) = rec["exc"][F1], rec["exc"][F2]
        d1 = 10 * np.log10(max(p1, 1e-30))
        d2 = 10 * np.log10(max(p2, 1e-30))
        check("S1 chain (EX path)", p1 > 1e-12 and p2 > 1e-12,
              f"f1 {d1:.1f} dBm (ref {ANCH_EX['dbm_f1']}), f2 {d2:.1f} dBm "
              f"(ref {ANCH_EX['dbm_f2']})")
    else:
        check("S1 chain (EX path)", False, "not measured")

    # S2 right quantity
    ref_p = float(np.mean([p for p, _ in rec["ref2"]])) if rec["ref2"] else 0.0
    rows0 = rec["two"].get(0.0, [])
    s0 = c12([(r[0], r[1]) for r in rows0])
    if s0:
        A1_0, A0_0, _ = s0
        check("S2 line at f2/2", A0_0 >= 50e-6,
              f"A0 = {A0_0*1e6:.0f} uW (anchor {ANCH['A0_uW']:.0f}; "
              f"x_ref {A0_0/max(ref_p,1e-30):.0f})")
        # S3 physics
        m = A1_0 / max(A0_0, 1e-30)
        check("S3a modulation", m >= 0.15, f"A1/A0 = {m:.2f} (anchor 0.42)")
        clo = c12([(r[0], r[1]) for r in rec["two"].get(0.5, [])]) \
            or c12([(r[0], r[1]) for r in rec["two"].get(1.0, [])])
        chi = c12([(r[0], r[1]) for r in rec["two"].get(8.0, [])])
        if clo and chi:
            r8 = chi[0] / max(clo[0], 1e-30)
            check("S3b memory decay", r8 < 0.5,
                  f"A1(8)/A1(0.5) = {r8:.2f} (anchor {ANCH['A1_ratio_8ns']})")
        # winding across all measured delays
        tt, phi0s = [], []
        for t_gap in sorted(rec["two"]):
            cc = c12([(r[0], r[1]) for r in rec["two"][t_gap]])
            if cc:
                tt.append(t_gap)
                phi0s.append(cc[2])
        if len(tt) >= 4:
            dphi = np.diff(np.unwrap(phi0s)) / np.diff(tt) * 180 / np.pi
            slew = float(np.mean(dphi))
            check("S3c phi_eff winding", 8.0 <= slew <= 39.0,
                  f"slew = {slew:.1f} deg/ns (anchor {ANCH['slew_deg_ns']}, "
                  f"theory 26.6)")
    else:
        check("S2 line at f2/2", False, "no overlap data")

    n_fail = sum(1 for _, ok in _RESULTS if not ok)
    print("-" * 66)
    print(f"VERDICT: {'GO' if n_fail <= 0 else 'NO-GO'} "
          f"({len(_RESULTS)-n_fail} PASS, {n_fail} FAIL)")
    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
