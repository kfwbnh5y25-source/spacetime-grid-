"""
Shared Hardware & Waveform Core for Spin-Wave Neuromorphic Reservoir Experiments.
Enforces all hardware & measurement-physics safeguards from NEUROMORPHIC_REPORT.md:
  - Forced SA RMS detector (DET:RAV), RMS averaging, manual RBW, 80 MHz span (2000-2080 MHz)
  - Never caches baseline >1 trial; stores RAW signal and RAW baseline separately in HDF5
  - Enforces 64 GSa/s, 1280-sample (20 ns) alignment, ZHL peak-crest protection (<= -5 dBm input)
  - Supports --dry-run offline simulation/verification mode without kNET hardware
"""

import time
import json
import numpy as np
from scipy.signal.windows import tukey

FS_HZ = 64e9
BLOCK_SAMPLES = 1280  # 20 ns AWG sequencer / memory quantum
MAX_PEAK_MV = 760.0   # Calibrated against threshold_reduction/run.py (<= 810 mV, <0.1 dB ZHL compression)


def align_samples(duration_ns: float) -> int:
    """Round duration in ns to nearest 1280-sample (20 ns) multiple."""
    n = int(round((duration_ns * 1e-9) * FS_HZ))
    return max(BLOCK_SAMPLES, int(round(n / BLOCK_SAMPLES)) * BLOCK_SAMPLES)


def make_pulse_waveform(
    symbols: list[dict],
    gap_ns: float = 2000.0,
    edge_ns: float = 1.0,
) -> np.ndarray:
    """
    Construct a periodic looped waveform from a list of symbol dicts:
      {'amp_mv': float, 'freq_mhz': float, 'start_ns': float, 'dur_ns': float, 'phase': float}
    Total loop length = max(end_ns) + gap_ns, aligned to 1280-sample (20 ns) blocks.
    Returns float64 array in Volts ([-1.0, 1.0]) normalized for awg.output_voltage = 1.0 V.
    """
    end_ns = max((s["start_ns"] + s["dur_ns"]) for s in symbols) if symbols else 200.0
    total_samples = align_samples(end_ns + gap_ns)
    t_s = np.arange(total_samples, dtype=np.float64) / FS_HZ
    wf = np.zeros(total_samples, dtype=np.float64)

    for s in symbols:
        amp_v = float(s["amp_mv"]) / 1e3
        if amp_v == 0.0:
            continue
        freq_hz = float(s["freq_mhz"]) * 1e6
        phase = float(s.get("phase", 0.0))
        i0 = int(round(s["start_ns"] * 1e-9 * FS_HZ))
        n_dur = int(round(s["dur_ns"] * 1e-9 * FS_HZ))
        i1 = min(total_samples, i0 + n_dur)
        if i1 <= i0:
            continue

        n_actual = i1 - i0
        alpha = min(1.0, (2.0 * edge_ns * 1e-9 * FS_HZ) / max(1, n_actual))
        env = tukey(n_actual, alpha=alpha)
        wf[i0:i1] += amp_v * env * np.sin(2.0 * np.pi * freq_hz * t_s[i0:i1] + phase)

    peak_mv = np.max(np.abs(wf)) * 1e3
    if peak_mv > MAX_PEAK_MV or np.max(np.abs(wf)) > 1.0:
        raise ValueError(f"Waveform peak {peak_mv:.1f} mV exceeds ZHL/DAC safety ceiling ({MAX_PEAK_MV} mV)")
    return wf


def hampel_despike(spec: np.ndarray, half_win: int = 7, n_sig: float = 5.0) -> np.ndarray:
    """15-bin Hampel median filter to remove narrowband 5-15 kHz RF spurs."""
    out = spec.copy()
    n = len(spec)
    for i in range(n):
        lo = max(0, i - half_win)
        hi = min(n, i + half_win + 1)
        win = spec[lo:hi]
        med = np.median(win)
        mad = np.median(np.abs(win - med))
        sigma = 1.4826 * mad
        if sigma > 0 and abs(spec[i] - med) > n_sig * sigma:
            out[i] = med
    return out


class HardwareRig:
    def __init__(self, h5_path: str, dry_run: bool = False, script_filename: str = __file__):
        self.dry_run = dry_run
        self.h5_path = h5_path
        self.trial_count = 0

        if not dry_run:
            from kNET import kNETInterface, HDF5Buffer
            self.knet = kNETInterface()
            self.knet.switch_to_sample()
            self.awg = self.knet.awg
            self.sa = self.knet.sa

            self.awg.abort_generation()
            self.awg.clear_data()
            self.awg.sampling_freq = FS_HZ
            self.awg.output_voltage = 1.0

            self.configure_sa_sample()
            self.hdf_buf = HDF5Buffer(
                h5_path,
                script_filename=script_filename,
                settings=self.sa.get_settings(),
                append=True,
                buffer_size=25,
            )
        else:
            self._dry_records = []
            self._sa_freq = np.linspace(2000e6, 2080e6, 1001)

    def configure_sa_sample(self, center_hz: float = 2040e6, span_hz: float = 80e6, rbw_hz: float = 500e3):
        if self.dry_run:
            self._sa_freq = np.linspace(center_hz - span_hz / 2, center_hz + span_hz / 2, 1001)
            return
        self.knet.switch_to_sample()
        self.sa.write(":FORM REAL,32")
        self.sa.centerFrequency = center_hz
        self.sa.spanFrequency = span_hz
        self.sa.points = 1001
        self.sa.resolutionBandwidthAutoState = "OFF"
        self.sa.resolutionBandwidth = rbw_hz
        self.sa.detectorTypeAutoState = "OFF"
        self.sa.detectorType = "RAV"
        self.sa.averageTypeAutoState = "OFF"
        self.sa.averageType = "RMS"
        self.sa.averageNumber = 3
        self.sa.sweepTime = 0.25
        self.sa.referenceLevel = -50.0
        self.sa.scaleDiv = 5.0

    def configure_sa_excitation(self, center_hz: float = 4100e6, span_hz: float = 200e6):
        if self.dry_run:
            self._sa_freq = np.linspace(center_hz - span_hz / 2, center_hz + span_hz / 2, 1001)
            return
        self.knet.switch_to_excitation()
        self.sa.centerFrequency = center_hz
        self.sa.spanFrequency = span_hz
        self.sa.points = 1001
        self.sa.resolutionBandwidthAutoState = "OFF"
        self.sa.resolutionBandwidth = 500e3
        self.sa.detectorTypeAutoState = "OFF"
        self.sa.detectorType = "RAV"
        self.sa.averageTypeAutoState = "OFF"
        self.sa.averageType = "RMS"
        self.sa.averageNumber = 1
        self.sa.sweepTime = 0.25

    def acquire_raw_trace(self) -> tuple[np.ndarray, np.ndarray]:
        """Returns (linear_power_mW, freq_Hz)."""
        if self.dry_run:
            rng = np.random.default_rng(self.trial_count + 1337)
            noise = np.abs(rng.normal(5e-8, 1.5e-8, size=1001))
            return noise, self._sa_freq.copy()
        spec_dbm = self.sa.trigAndWait()
        freqs = self.sa.get_frequencies()
        return np.power(10.0, spec_dbm / 10.0), freqs

    def measure_trial(self, key: str, wf: np.ndarray, meta: dict, measure_baseline: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Executes a single trial:
          1. Acquires fresh baseline (AWG output off) — never cached across trials.
          2. Uploads & starts waveform on M8195A.
          3. Acquires RMS spectrum.
          4. Aborts & clears AWG.
          5. Saves RAW signal and RAW baseline separately.
        """
        self.trial_count += 1
        if measure_baseline:
            if not self.dry_run:
                self.awg.abort_generation()
            bl_lin, freqs = self.acquire_raw_trace()
        else:
            bl_lin = np.zeros(1001, dtype=np.float64)
            freqs = self._sa_freq.copy() if self.dry_run else self.sa.get_frequencies()

        if not self.dry_run:
            self.awg.send_waveform(wf)
            self.awg.start_generation()
            sig_lin, freqs = self.acquire_raw_trace()
            self.awg.abort_generation()
            self.awg.clear_data()
        else:
            # Synthetic signal proportional to waveform energy for dry-run verification
            rng = np.random.default_rng(self.trial_count + 42)
            p_wf = float(np.mean(wf**2))
            sig_lin = bl_lin + 1e-4 * p_wf * np.exp(-((freqs - 2035e6) / 4e6) ** 2) * (1.0 + rng.normal(0, 0.03, 1001))

        if not self.dry_run:
            self.hdf_buf.append_data(f"{key}__SIG", freqs, sig_lin)
            if measure_baseline:
                self.hdf_buf.append_data(f"{key}__BL", freqs, bl_lin)
        else:
            self._dry_records.append({"key": key, "meta": meta, "sig_sum": float(np.sum(sig_lin - bl_lin))})

        return sig_lin, bl_lin, freqs

    def close(self):
        if not self.dry_run:
            self.awg.abort_generation()
            self.awg.clear_data()
            self.hdf_buf.flush_to_hdf5()
        else:
            with open(self.h5_path + ".dryrun.json", "w") as f:
                json.dump({"trials": len(self._dry_records), "sample": self._dry_records[:10]}, f, indent=2)
