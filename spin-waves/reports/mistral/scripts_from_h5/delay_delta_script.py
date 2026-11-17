import time

import numpy as np
from tqdm import tqdm

from kNET import kNETInterface

np.set_printoptions(suppress=True)
from scipy.signal.windows import tukey
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
awg.clear_data()
sa.points = 1001
sa.averageNumber = 1
sa.sweepTime = 1 / 2
sa.spanFrequency = 30e6
sa.resolutionBandwidthAutoState = True
sa.resolutionBandwidth = 1e6
sa.referenceLevel = -53.8
sa.scaleDiv = 0.5
print(sa.resolutionBandwidth)

_cached_bl, _cached_bl_time = None, None
_BASELINE_MAX_AGE = 90


def get_baseline():
    global _cached_bl, _cached_bl_time

    now = time.time()

    # If baseline older than 1 hour â†’ remeasure
    if _cached_bl is None or (now - _cached_bl_time) > _BASELINE_MAX_AGE:
        bl = np.power(10.0, sa.trigAndWait() / 10)
        sa_freq = sa.get_frequencies()
        bl[(sa_freq > sa_mask[0]) & (sa_freq < sa_mask[1])] = 0
        _cached_bl = bl
        _cached_bl_time = now
        print(f"[Baseline] Recalculated at {time.ctime(now)}")

    return _cached_bl


hdf_buf = HDF5Buffer(
    r"D:\Python_GIT\RF_probestation_1\experiments\BiYIG_3mum_F7_NEW\40G\experiments_dedalo_02062025\2_adapted\delay_delta.h5",
    script_filename=__file__,
    settings=sa.get_settings(),
    append=True,
)
keys = hdf_buf.get_keys()
print(len(keys))

alpha = 0.01

amps = np.array([327, 588])


def waveform(delay, f1, f2, phase):

    t_full = np.arange(128000) / 64e9
    lf1_steps = round(300 * 64)
    lf2_steps = round(400 * 64)
    delay_steps = round(delay * 64)

    wf = np.zeros_like(t_full)
    window1, window2 = np.zeros_like(t_full), np.zeros_like(t_full)
    if lf1_steps != 0:
        window1[64000 - lf1_steps : 64000] = tukey(lf1_steps, 64 / lf1_steps)
    if lf2_steps != 0:
        window2[64000 + delay_steps : 64000 + lf2_steps + delay_steps] = tukey(
            lf2_steps, 64 / lf2_steps
        )
    windows = [window1, window2]
    phases = [0, phase]
    freqs_ = [f1, f2]

    for i in range(len(amps)):
        wf += (amps[i] / 1e3 * windows[i]) * np.sin(
            2 * np.pi * freqs_[i] * 1e6 * t_full + phases[i]
        )
    if np.abs(wf).max() > 1:
        raise AssertionError

    return wf


def measure(delay, f1, f2, phase, readoff, force_eval=False):

    KKEY = f"{delay:.4f}_{round(f2-f1)}_{phase:.3f}_{readoff}"

    if KKEY in keys and not force_eval:
        return

    wf = waveform(delay, f1, f2, phase)
    bl = get_baseline()

    awg.send_waveform(wf)

    awg.start_generation()
    spec, sa_freq = np.power(10.0, sa.trigAndWait() / 10), sa.get_frequencies()
    spec[(sa_freq > sa_mask[0]) & (sa_freq < sa_mask[1])] = 0

    spec -= bl
    print(delay, f1, f2, phase, round(spec.mean() * 1e10))
    awg.abort_generation()
    awg.clear_data()

    hdf_buf.append_data(KKEY, sa_freq, spec)
    keys.append(KKEY)
    return spec.mean() * 1e9


def main():
    global _cached_bl, _cached_bl_time
    _cached_bl = None
    _cached_bl_time = None
    sigma = 1 / 2 * (4070 + 4142)
    deltas = np.arange(50, 90, 2)
    delays = np.arange(-750 / 64, 7, 5 / 64)[::-1]
    delays = delays[delays > -5]

    phases = np.arange(0, 2 * np.pi, 2 * np.pi / 16)
    np.random.shuffle(phases)
    np.random.shuffle(delays)
    np.random.shuffle(deltas)
    f1s = sigma - deltas / 2
    f2s = sigma + deltas / 2

    for i, f2 in tqdm(enumerate(f2s), total=len(f2s)):
        sa.centerFrequency = f2 / 2 * 1e6
        assert np.isclose(1 / 2 * (f1s[i] + f2s[i]), sigma, rtol=1e-5)
        for delay in delays:
            for phase in phases:
                measure(delay, f1s[i], f2s[i], phase, 1)

    _cached_bl = None
    hdf_buf.flush_to_hdf5()
    exit()


if __name__ == "__main__":
    main()
