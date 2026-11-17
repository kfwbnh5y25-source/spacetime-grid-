import sys
import numpy as np

sys.path.insert(0, "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/measure")
import hw_core_spectro as hc
from hw_core_spectro import LoopSpec, Pulse, build_loop

hc.MAX_PEAK_MV = 1e6  # scan only
F1, F2 = 4075.0, 4149.0
A1M, A2M = 327.0, 588.0
D1_NS, D2_NS = 300 / 64.0, 400 / 64.0


def wf_two_tone(phi, t_gap):
    s2 = D1_NS + t_gap
    sh = max(0.0, -s2)
    spec = LoopSpec(pulses=[Pulse(A1M, F1, sh, D1_NS, 0.0),
                            Pulse(A2M, F2, s2 + sh, D2_NS, phi)],
                    t0_ref_ns=0.0, gap_ns=20.0)
    return build_loop(spec, edge_ns=1.0)


phases = np.arange(0, 2 * np.pi, np.pi / 60)
peaks = []
for t_gap in (-2.0, -0.5, 0.0, 1.0, 2.0, 5.0, 8.0):
    for phi in phases:
        peaks.append((wf_two_tone(phi, t_gap).max() * 1e3, t_gap, phi))
peaks.sort(reverse=True)
print("worst 8 (peak mV, t_gap, phi):")
for p in peaks[:8]:
    print(f"  {p[0]:7.1f}  {p[1]:+.2f} ns  phi={p[2]:.3f}")
print("phi=0 peaks:", [round(wf_two_tone(0.0, t).max() * 1e3, 1) for t in (0.0, 1.0)])
