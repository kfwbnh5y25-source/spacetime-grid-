import time

import numpy as np

from kNET import kNETInterface

np.set_printoptions(suppress=True)
from kNET import HDF5Buffer

sa_mask = (2.110e9, 2.170e9)

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
sa.averageNumber = 5
sa.sweepTime = 1 / 5
sa.spanFrequency = 250e3 * 10
sa.resolutionBandwidthAutoState = False
sa.resolutionBandwidth = 200e3
sa.referenceLevel = -38.8
sa.scaleDiv = 4

_cached_bl, _cached_bl_time = None, None
_BASELINE_MAX_AGE = 90
from scipy.signal.windows import tukey

ramp_up = tukey(128, 1)[:64]
ramp_down = np.flip(ramp_up)


# 4108 350
def get_baseline():
    global _cached_bl, _cached_bl_time

    now = time.time()

    # If baseline older than 1 hour → remeasure
    if _cached_bl is None or (now - _cached_bl_time) > _BASELINE_MAX_AGE:
        bl = np.power(10.0, sa.trigAndWait() / 10)
        sa_freq = sa.get_frequencies()
        # bl[(sa_freq > sa_mask[0]) & (sa_freq < sa_mask[1])] = 0
        _cached_bl = bl
        _cached_bl_time = now

    return _cached_bl


hdf_buf = HDF5Buffer(
    "data_broad.h5",
    script_filename=__file__,
    settings=sa.get_settings(),
    append=True,
    buffer_size=50,
)
keys = hdf_buf.get_keys()


def waveforms_long(a2, f2, phase, f1_on, f2_on):
    contact_point = 2.999e-6
    t = np.arange(3 * 64000) / 64e9

    wfs = [np.empty(64, dtype=np.float32), np.empty(64, dtype=np.float32)]
    amps = [a2, 266]
    freqs = [f2, 4130]
    phases = [phase / 2, -phase / 2]
    on = [f1_on, f2_on]

    for i in range(2):
        wfs[i] = (
            on[i]
            * amps[i]
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


def measure_long(a2, f2, phase, onlyf2=False, force=False):

    label = 0
    f1_on = True
    f2_on = True

    if onlyf2:
        f1_on = False
        label = 2

    KKEY = f"{np.round(a2,1):.1f}_{np.round(f2,1):.1f}_{np.round(phase,3):.3f}_{label}_4130"
    if KKEY in keys and not force:
        return

    wfA, wfB, wfC = waveforms_long(a2, f2, phase, f1_on, f2_on)

    bl = get_baseline()
    awg.send_waveform(wfA, segment_id=1)
    awg.send_waveform(wfB, segment_id=2)
    awg.send_waveform(wfC, segment_id=3)

    awg.start_sequence(
        (1, 2, 3),
        (1, 1, 1),
    )
    awg.start_generation()
    spec, sa_freq = np.power(10.0, sa.trigAndWait() / 10) - bl, sa.get_frequencies()

    awg.abort_generation()

    hdf_buf.append_data(KKEY, sa_freq, spec)
    print(KKEY, spec.mean() * 1e10)
    keys.append(KKEY)
    return spec.mean() * 1e10


def main():
    global _cached_bl, _cached_bl_time
    from scipy.stats import qmc

    _cached_bl = None
    _cached_bl_time = None

    a2_min1, a2_max1 = 280, 350
    f2_min1, f2_max1 = 4050, 4090

    a2_min, a2_max = 220, 350
    f2_min, f2_max =3900,4140

    m = 14
    n_points = 2**m

    sampler = qmc.Sobol(d=2, scramble=True, seed=42)
    samples = sampler.random_base2(m=m)

    a2 = a2_min + samples[:, 0] * (a2_max - a2_min)
    f2 = f2_min + samples[:, 1] * (f2_max - f2_min)

    f2 = np.round(f2)

    a2 = np.round(a2)

    a22 = a2_min1 + samples[:, 0] * (a2_max1 - a2_min1)
    f22 = f2_min1 + samples[:, 1] * (f2_max1 - f2_min1)

    a22 = np.round(a22)
    f22 = np.round(f22 )

    a2 = np.concatenate((a2[:2000], a22))
    f2 = np.concatenate((f2[:2000], f22))
    a2[1000:2000]-=50
    n_points = len(f2)
    phases = np.arange(0, 2 * np.pi, 2 * np.pi / 4)



    from tqdm import trange

    sa.centerFrequency = 4130 * 1e6 / 2

    for j in range(10):
        for jjj in range(len(phases)):
            measure_long(324, 4069, phases[jjj], onlyf2=False, force=True)
        measure_long(324, 4069, 0, onlyf2=True, force=True)

    _cached_bl, _cached_bl_time = None, None
    sa.centerFrequency =4130 * 1e6 / 2
    for i in trange(4981, n_points):


        np.random.shuffle(phases)

        measure_long(0, f2[i], 0)



        measure_long(a2[i], f2[i], 0, onlyf2=True)


        for j in range(len(phases)):
            measure_long(a2[i], f2[i], phases[j])



if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt as e:
        print(e)
        hdf_buf.flush_to_hdf5()
        exit()
