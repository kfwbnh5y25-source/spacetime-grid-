import time

import numpy as np

from kNET import kNETInterface

np.set_printoptions(suppress=True)
from kNET import HDF5Buffer
from scipy.signal.windows import tukey

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


def get_baseline():
    global _cached_bl, _cached_bl_time

    now = time.time()

    # If baseline older than 1 hour â†’ remeasure
    if _cached_bl is None or (now - _cached_bl_time) > _BASELINE_MAX_AGE:
        bl = np.power(10.0, sa.trigAndWait() / 10)

        _cached_bl = bl
        _cached_bl_time = now

    return _cached_bl


hdf_buf = HDF5Buffer(
        "f1_4069_5_a1_a2_f2_sweep.h5",
        script_filename=__file__,
        settings=sa.get_settings(),
        append=True,
        buffer_size=20,
)
print(8)
keys = hdf_buf.get_keys()
print(9)

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

    freqs0 = [4069, 4127]
    amps0 = [324, 266]
    phases = [4.55, 1.31]

    f1 = 4069
    a1 = [200, 250, 268, 279, 286, 298, 309, 324]
    a1=[324]
    a2s = []
    f2s = []

    a2_min, a2_max = 200, 275
    a2_min, a2_max = 210, 250

    f2_min, f2_max = 4075, 4150

    a2 = np.arange(a2_min, a2_max + 1)
    f2 = np.arange(f2_min, f2_max + 1)
    np.random.shuffle(a2)
    np.random.shuffle(f2)

    from tqdm import trange

    for i in trange(len(f2)):
        _cached_bl, _cached_bl_time = None, None
        sa.centerFrequency = f2[i] * 1e6 / 2

        for k in trange(len(a2)):
            for j in trange(len(a1)):

                for l in range(2):
                    measure_long(a1[j], a2[k], f1, f2[i], phases[l], force=True)
            measure_long(0, a2[k], 0, f2[i], 0, force=True)

    hdf_buf.flush_to_hdf5()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt as e:
        print(e)
        hdf_buf.flush_to_hdf5()
        exit()
