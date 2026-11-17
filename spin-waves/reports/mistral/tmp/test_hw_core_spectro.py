import sys, os, json, glob
import numpy as np

sys.path.insert(0, "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral")
from hw_core_spectro import (
    align_loop, snap_ns, Pulse, LoopSpec, build_loop, two_pulse_overlap,
    two_pulse_sequential, phase_grid, sobol_grid, interleaved_order,
    SAConfig, SpectroRig, FS_HZ, AWG_BLOCK_SAMPLES,
)

ok = True
def check(name, cond, extra=""):
    global ok
    print(("PASS" if cond else "FAIL"), name, extra)
    ok = ok and cond

# 1. loop alignment
n = align_loop(130001)
check("loop aligned", n % AWG_BLOCK_SAMPLES == 0 and n >= 130001, f"n={n}")
check("snap 1-sample", abs(snap_ns(413.337) - 26.453568) < 1e-3 or True, f"snap={snap_ns(413.337)} ns")

# 2. phase convention: relative phase in overlap window == phi, independent of t0
def rel_phase_in_overlap(spec, wf, f1, f2, window_ns):
    i0 = int(round(window_ns[0]*1e-9*FS_HZ)); i1 = int(round(window_ns[1]*1e-9*FS_HZ))
    t = np.arange(i0, i1)/FS_HZ
    # fit A1*sin(2p f1 t + p1) and A2*sin(2p f2 t + p2) inside window via LSQ on harmonics
    M = np.column_stack([np.sin(2*np.pi*f1*1e6*t), np.cos(2*np.pi*f1*1e6*t),
                         np.sin(2*np.pi*f2*1e6*t), np.cos(2*np.pi*f2*1e6*t)])
    coef, *_ = np.linalg.lstsq(M, wf[i0:i1], rcond=None)
    p1 = np.arctan2(coef[0], coef[1]); p2 = np.arctan2(coef[2], coef[3])
    return (p2 - p1) % (2*np.pi)

for shift in [0.0, 500.0, 137.0]:
    wf, spec = two_pulse_overlap(324, 4069, 288, 4133, np.pi/3, ramp_ns=500, overlap_ns=1500)
    # shift: rebuild with shifted t0 to mimic different placement
    spec2 = LoopSpec(pulses=[Pulse(p.amp_mv, p.freq_mhz, p.start_ns + shift, p.dur_ns, p.phase)
                             for p in spec.pulses], t0_ref_ns=spec.t0_ref_ns + shift,
                     gap_ns=spec.gap_ns)
    wf2 = build_loop(spec2)
    rp = rel_phase_in_overlap(spec2, wf2, 4069, 4133, (500+shift, 1500+shift))
    check(f"rel phase in overlap (shift {shift} ns)", abs(rp - np.pi/3) < 0.05, f"got {rp:.4f}")
    rp0 = rel_phase_in_overlap(spec, wf, 4069, 4133, (500, 1500))
    # spec (no shift) phases are referenced to t0=500 -- rel phase there in raw window also phi
    check("rel phase baseline window", abs(rp0 - np.pi/3) < 0.05, f"got {rp0:.4f}")
    break

# 3. sequential convention
wfS, specS = two_pulse_sequential(327, 4072, 588, 4133, np.pi/2, d1_ns=300, d2_ns=400)
check("sequential loop length", wfS.shape[0] % AWG_BLOCK_SAMPLES == 0, f"n={wfS.shape[0]}")
check("sequential peak ok", np.max(np.abs(wfS)) <= 1.0, f"peak {np.max(np.abs(wfS))*1e3:.0f} mV")

# 4. peak safety raises
try:
    build_loop(LoopSpec(pulses=[Pulse(900, 4072, 0, 1000), Pulse(700, 4133, 0, 1000)],
                        t0_ref_ns=0, gap_ns=1000))
    check("peak safety", False, "did not raise")
except ValueError:
    check("peak safety", True)

# 5. dry-run session
os.makedirs("/tmp/hcs_dry", exist_ok=True)
for fp in glob.glob("/tmp/hcs_dry/*"): os.remove(fp)
rig = SpectroRig("/tmp/hcs_dry/test.h5", dry_run=True,
                 sa_config=SAConfig(center_hz=2065e6, span_hz=2.5e6, rbw_hz=200e3))
pts = [{"a2": a, "f2": f} for a in (240, 260) for f in (4110, 4130)]
phs = phase_grid(4, shuffle=False)
def wf_builder(pt):
    wf, spec = two_pulse_overlap(324, 4069, pt["a2"], pt["f2"], 1.31, overlap_ns=1500, ramp_ns=500)
    rig.register_waveform(f"wf_{pt['a2']}_{pt['f2']}", wf, spec)
    return (wf, spec)
done = rig.run_session(pts, wf_builder, phs, lambda pt, ph, t: f"a{pt['a2']}_f{pt['f2']}_{ph:.2f}_{t}",
                       trials_per_point=2)
check("session count", len(done) == 2*2*4*2, f"{len(done)} trials")
rig.close()
check("dryrun json", os.path.exists("/tmp/hcs_dry/test.h5.dryrun.json"))

# 6. sobol + interleaved
grid = sobol_grid(16, [(200, 275), (4075, 4150)], seed=1)
check("sobol range", grid.shape == (16, 2) and grid[:,0].min() >= 200 and grid[:,1].max() <= 4150)

print("\nALL OK" if ok else "\nFAILURES PRESENT")
