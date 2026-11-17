import sys, numpy as np
sys.path.insert(0, "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/neuromorphic")
from hw_core import align_samples, make_pulse_waveform, hampel_despike, HardwareRig

# 1. dry-run full cycle
rig = HardwareRig("/tmp/dryrun_test.h5", dry_run=True)
rig.configure_sa_sample()
sym = [
    {"amp_mv": 327, "freq_mhz": 4072, "start_ns": 0, "dur_ns": 500, "phase": 0.0},
    {"amp_mv": 288, "freq_mhz": 4133, "start_ns": 1000, "dur_ns": 500, "phase": 1.31},
]
wf = make_pulse_waveform(sym, gap_ns=2000)
print("wf shape", wf.shape, "=> loop period us", wf.shape[0]/64e9*1e6)
print("align_samples(413ns) =", align_samples(413), "samples (should be 20ns multiples)")
print("align_samples(20ns) =", align_samples(20), "samples")
sig, bl, fr = rig.measure_trial("test__1", wf, {"trial":1}, measure_baseline=True)
print("dry trial ok; sig shape", sig.shape, "delta-sum %.4g" % float(np.sum(sig-bl)))
rig.close()
print("dry run closed")

# 2. energy check: overlap nonlinearity sanity - can two tukey pulses overlap?
over = [
    {"amp_mv": 327, "freq_mhz": 4072, "start_ns": 0, "dur_ns": 1000, "phase": 0.0},
    {"amp_mv": 288, "freq_mhz": 4133, "start_ns": 500, "dur_ns": 1000, "phase": -1.31/2},
]
wf2 = make_pulse_waveform(over, gap_ns=1000)
print("overlap wf peak mV %.1f" % (np.max(np.abs(wf2))*1e3))

# 3. hampel
sp = np.zeros(1001); sp[500] = 100.0; sp += np.abs(np.random.default_rng(0).normal(0,1e-3,1001))
d = hampel_despike(sp)
print("hampel spike removed:", abs(d[500]) < 1.0)
