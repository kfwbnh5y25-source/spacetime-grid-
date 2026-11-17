import time

import numpy as np
from tqdm import tqdm

from kNET import kNETInterface

np.set_printoptions(suppress=True)
from scipy.signal.windows import tukey
from kNET import HDF5Buffer

sa_mask = (2.110e9, 2.170e9)

knet = kNETInterface()
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

_cached_bl, _cached_bl_time = None, None
_BASELINE_MAX_AGE = 90


def get_baseline():
    global _cached_bl, _cached_bl_time

    now = time.time()

    # If baseline older than 1 hour â†’ remeasure
    if _cached_bl is None or (now - _cached_bl_time) > _BASELINE_MAX_AGE:
        bl = np.power(10.0, sa.trigAndWait() / 10)
        sa_freq = sa.get_frequencies()
        # bl[(sa_freq > sa_mask[0]) & (sa_freq < sa_mask[1])] = 0
        _cached_bl = bl
        _cached_bl_time = now
        print(f"[Baseline] Recalculated at {time.ctime(now)}")

    return _cached_bl


hdf_buf = HDF5Buffer(
    r"D:\Python_GIT\RF_probestation_1\experiments\BiYIG_3mum_F7_NEW\40G\experiments_dedalo_02062025\2_adapted\detune_fine031126.h5",
    script_filename=__file__,
    settings=sa.get_settings(),
    append=False,
)
keys = hdf_buf.get_keys()
print(len(keys))

alpha = 0.01

amps = np.array([390, 588])
freqs = np.array([4072, 4133])


def waveform(f1, f2, phase):

    t_full = np.arange(128000) / 64e9
    lf1_steps = round(300 * 64)
    lf2_steps = round(400 * 64)

    wf = np.zeros_like(t_full)
    window1, window2 = np.zeros_like(t_full), np.zeros_like(t_full)
    if lf1_steps != 0:
        window1[64000 - lf1_steps : 64000] = tukey(lf1_steps, 64 / lf1_steps)
    if lf2_steps != 0:
        window2[64000 : 64000 + lf2_steps] = tukey(lf2_steps, 64 / lf2_steps)
    windows = [window1, window2]
    phases = [0, phase]
    freqs = [f1, f2]

    for i in range(len(amps)):
        wf += (amps[i] / 1e3 * windows[i]) * np.sin(
            2 * np.pi * freqs[i] * 1e6 * t_full + phases[i]
        )
    if np.abs(wf).max() > 1:
        raise AssertionError

    return wf


def measure(f1, f2, phase, readoff):

    KKEY = f"{f1:.2f}_{f2:.2f}_{phase:.3f}_{readoff}"

    if KKEY in keys:
        pass

    wf = waveform(f1, f2, phase)
    bl = get_baseline()

    awg.send_waveform(wf)
    awg.start_generation()
    spec, sa_freq = np.power(10.0, sa.trigAndWait() / 10), sa.get_frequencies()
    # spec[(sa_freq > sa_mask[0]) & (sa_freq < sa_mask[1])] = 0

    spec -= bl
    print(f1, f2, phase, round(spec.mean() * 1e10))
    awg.abort_generation()
    awg.clear_data()

    hdf_buf.append_data(KKEY, sa_freq, spec)
    keys.append(KKEY)
    return spec.mean() * 1e9


def main():
    global _cached_bl, _cached_bl_time
    _cached_bl = None
    _cached_bl_time = None

    freqs = np.array([4072, 4133])

    f1s = np.arange(4060, 4085)
    f2s = np.arange(4120, 4155)



    phases = np.array(np.arange(0, 2 * np.pi, 2 * np.pi / 16).tolist())
    # np.random.shuffle(phases)
    np.random.shuffle(f1s)
    np.random.shuffle(f2s)

    knet.switch_to_sample()
    # readoff =0

    for ll in range(10):
        for f2 in tqdm(f2s, total=len(f2s)):

            sa.centerFrequency = f2 / 2 * 1e6
            _cached_bl = None

            for f1 in tqdm(f1s, total=len(f1s)):
                for phase in phases:
                    measure(f1, f2, phase, 1)

    exit()

    _cached_bl = None
    hdf_buf.flush_to_hdf5()
    min_freq = f1s.min()
    max_delta = (f2s[None] - f1s[:, None]).max()
    max_freq = f2s.max() + max_delta
    sa.spanFrequency = (max_freq - min_freq) / 2 * 1e6 + 40e6
    sa.centerFrequency = (max_freq + min_freq) / 4 * 1e6

    print(
        min_freq,
        sa.centerFrequency / 1e9,
        sa.spanFrequency / 1e6,
        (f2s.max() - f1s.min()) / 2,
    )
    sa.resolutionBandwidth = 10e6
    _cached_bl = None

    for f2 in tqdm(f2s, total=len(f2s)):

        for f1 in tqdm(f1s, total=len(f1s)):
            for phase in phases:
                measure(f1, f2, phase, "TOT")
    exit()
    sa.centerFrequency = freqs[0] / 2 * 1e6

    for f1 in tqdm(lf1s, total=len(lf1s)):
        for f2 in lf2s:
            measure(f1, f2, 0)
    _cached_bl = None

    _cached_bl, _cached_bl_time = None, None
    awg.abort_generation()
    knet.switch_to_excitation()
    sa.averageNumber = 10
    sa.sweepTime = 10e-3
    sa.spanFrequency = 60e6
    sa.resolutionBandwidth = 1e6
    sa.centerFrequency = freqs[1] * 1e6
    sa.referenceLevel = -20
    sa.scaleDiv = 5

    for f1 in tqdm(lf1s, total=len(lf1s)):
        for f2 in lf2s:
            measure(f1, f2, "1EX")
    _cached_bl = None
    sa.centerFrequency = freqs[0] * 1e6
    for f1 in tqdm(lf1s, total=len(lf1s)):
        for f2 in lf2s:
            measure(f1, f2, "0EX")

    hdf_buf.flush_to_hdf5()


if __name__ == "__main__":
    main()
