import argparse
import time

import numpy as np

from kNET import kNETInterface

np.set_printoptions(suppress=True)
from kNET import HDF5Buffer
from scipy.signal.windows import tukey

# 16-phase CW overlap sweep: phase-response curve for Fig 3 (overlapping pulses).
# Same 3-segment AWG sequence as a2_266_5_f2_a1_f1_sweep.py:
#   seg A (1 us): f1 always on (ramped), f2 turns on in the last 1 ns
#   seg B (15 us): both tones on continuously  (overlap region)
#   seg C (1 us): both tones removed (1 ns ramp down)
# Relative carrier phase = phase (applied symmetrically as +phase/2, -phase/2).

knet = kNETInterface()
knet.switch_to_sample()

awg = knet.awg
awg.output_voltage = 1.0
awg.clear_data()

sa = knet.sa
sa.write(":FORM REAL,32")
awg.sampling_freq = 64e9
awg.abort_generation()
sa.points = 1001
sa.averageNumber = 2
sa.sweepTime = 1 / 5
sa.spanFrequency = 250e3 * 10
sa.resolutionBandwidthAutoState = False
sa.resolutionBandwidth = 200e3
sa.referenceLevel = -38.8
sa.scaleDiv = 4

_cached_bl, _cached_bl_time = None, None
_BASELINE_MAX_AGE = 90

ramp_up = tukey(128, 1)[:64]
ramp_down = np.flip(ramp_up)

parser = argparse.ArgumentParser(description="16-phase CW overlap phase-response sweep")
parser.add_argument("--out", default="phase_sweep_overlap.h5")
parser.add_argument("--a1-list", default="240,280,320,340",
                    help="comma-separated a1 values (mV); pick below/at/above the two-tone knee")
parser.add_argument("--a2", type=float, default=266.0)
parser.add_argument("--f1-min", type=float, default=4050.0)
parser.add_argument("--f1-max", type=float, default=4095.0)
parser.add_argument("--f1-step", type=float, default=5.0)
parser.add_argument("--f2-list", default="4102,4104,4117,4124,4136,4144")
parser.add_argument("--n-phase", type=int, default=16)
parser.add_argument("--reps", type=int, default=2)
parser.add_argument("--force", action="store_true")
args = parser.parse_args()


def get_baseline():
    global _cached_bl, _cached_bl_time
    now = time.time()
    if _cached_bl is None or (now - _cached_bl_time) > _BASELINE_MAX_AGE:
        bl = np.power(10.0, sa.trigAndWait() / 10)
        _cached_bl = bl
        _cached_bl_time = now
    return _cached_bl


hdf_buf = HDF5Buffer(
    args.out,
    script_filename=__file__,
    settings=sa.get_settings(),
    append=True,
    buffer_size=20,
)
keys = hdf_buf.get_keys()
print(len(keys))


def waveforms_long(a1, a2, f1, f2, phase):
    contact_point = 2.999e-6
    t = np.arange(3 * 64000) / 64e9

    wfs = [np.empty(64, dtype=np.float32), np.empty(64, dtype=np.float32)]
    amps = [a1, a2]
    freqs = [f1, f2]
    phases = [phase / 2, -phase / 2]

    for i in range(2):
        wfs[i] = (
            amps[i]
            / 1000
            * np.sin(2 * np.pi * freqs[i] * 1e6 * (t - contact_point) + phases[i])
        )

    envA_f1 = np.ones_like(t)
    envA_f1[:64] = ramp_up

    envA_f2 = np.zeros_like(t)
    envA_f2[-64:] = ramp_up

    envC = np.zeros_like(t)
    envC[:64] = ramp_down

    wfA = envA_f1 * wfs[0] + envA_f2 * wfs[1]
    wfB = wfs[0] + wfs[1]
    wfC = envC * wfs[0] + envC * wfs[1]

    if np.abs(wfA).max() > 1:
        raise AssertionError
    if np.abs(wfB).max() > 1:
        raise AssertionError
    if np.abs(wfC).max() > 1:
        raise AssertionError

    return wfA, wfB, wfC


def measure_long(a1, a2, f1, f2, phase, force=False, label=1):
    a1 = np.round(a1)
    a2 = np.round(a2)
    f1 = np.round(f1)
    f2 = np.round(f2)

    KKEY = f"{a1:.0f}_{a2:.0f}_{f1:.0f}_{f2:.0f}_{np.round(phase, 3):.3f}_{label}"
    if KKEY in keys and not force:
        return

    wfA, wfB, wfC = waveforms_long(a1, a2, f1, f2, phase)

    bl = get_baseline()
    awg.send_waveform(wfA, segment_id=1)
    awg.send_waveform(wfB, segment_id=2)
    awg.send_waveform(wfC, segment_id=3)

    awg.start_sequence(
        (1, 2, 3),
        (1, 15, 1),
    )
    awg.start_generation()
    spec, sa_freq = np.power(10.0, sa.trigAndWait() / 10) - bl, sa.get_frequencies()

    awg.abort_generation()

    hdf_buf.append_data(KKEY, sa_freq, spec)
    print(KKEY, spec.mean() * 1e10)
    keys.append(KKEY)
    return


def main():
    global _cached_bl, _cached_bl_time

    _cached_bl = None
    _cached_bl_time = None

    a1_vals = np.array([float(v) for v in args.a1_list.split(",")])
    f2_vals = np.array([float(v) for v in args.f2_list.split(",")])
    f1_vals = np.arange(args.f1_min, args.f1_max + 1, args.f1_step)
    phases = np.arange(0, 2 * np.pi, 2 * np.pi / args.n_phase)

    n_pts = len(a1_vals) * len(f1_vals) * len(f2_vals) * args.n_phase * args.reps
    est_min = n_pts * (sa.sweepTime * sa.averageNumber + 2.0) / 60.0
    print(f"Planned points: {n_pts} (~{est_min:.0f} min)")

    # Interleave: for each (f2, a1, f1) point sweep all phases in random order
    # so slow drift does not alias onto a fixed phase.
    rng = np.random.default_rng(2026)
    jobs = []
    for f2 in f2_vals:
        for a1 in a1_vals:
            for f1 in f1_vals:
                for ph in phases:
                    for rep in range(1, args.reps + 1):
                        jobs.append((a1, f2, f1, ph, rep))
    rng.shuffle(jobs)

    from tqdm import tqdm

    for a1, f2, f1, ph, rep in tqdm(jobs, total=len(jobs)):
        cfl = sa.centerFrequency
        want = f2 / 2 * 1e6
        if abs(cfl - want) > 1.0:
            sa.centerFrequency = want
            _cached_bl = None
            _cached_bl_time = None
        measure_long(a1, args.a2, f1, f2, ph, force=args.force, label=rep)

    hdf_buf.flush_to_hdf5()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt as e:
        print(e)
        hdf_buf.flush_to_hdf5()
        exit()
