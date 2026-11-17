import sys
import numpy as np
sys.path.insert(0, "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral")
from hw_core_spectro import Pulse, LoopSpec, build_loop, two_pulse_overlap, FS_HZ

def rel_phase_in_overlap(spec, wf, f1, f2, window_ns):
    i0 = int(round(window_ns[0]*1e-9*FS_HZ)); i1 = int(round(window_ns[1]*1e-9*FS_HZ))
    t = np.arange(i0, i1)/FS_HZ
    t = t - t[0]  # origin at overlap window start: physical relative phase
    M = np.column_stack([np.sin(2*np.pi*f1*1e6*t), np.cos(2*np.pi*f1*1e6*t),
                         np.sin(2*np.pi*f2*1e6*t), np.cos(2*np.pi*f2*1e6*t)])
    coef, *_ = np.linalg.lstsq(M, wf[i0:i1], rcond=None)
    return (np.arctan2(coef[2], coef[3]) - np.arctan2(coef[0], coef[1])) % (2*np.pi)

phi = 1.31
for shift in [0.0, 100.0, 500.0, 137.3]:
    wf, spec = two_pulse_overlap(324, 4069, 288, 4133, phi, ramp_ns=500, overlap_ns=1500)
    spec2 = LoopSpec(pulses=[Pulse(p.amp_mv, p.freq_mhz, p.start_ns+shift, p.dur_ns, p.phase)
                             for p in spec.pulses], t0_ref_ns=spec.t0_ref_ns+shift, gap_ns=spec.gap_ns)
    wf2 = build_loop(spec2)
    rp = rel_phase_in_overlap(spec2, wf2, 4069, 4133, (500+shift, 1500+shift))
    print(f"shift {shift:7.1f} ns -> rel phase {rp:.4f} rad (want {phi:.2f}) error {abs(rp-phi):.5f}")
