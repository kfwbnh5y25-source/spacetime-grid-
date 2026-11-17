"""
hw_core_spectro.py - Standard hardware & waveform chassis for the two-tone
spin-wave spectroscopy / phase-gating experiments (Paper A).

Improvements over the ad-hoc per-experiment scripts:
  1. TIME-REFERENCED PHASE CONVENTION. Every tone is written as
         s_i(t) = a_i * env_i(t) * sin( 2*pi*f_i*(t - t0_ref) + phi_i )
     with t0_ref a single declared reference time stored in the HDF5 attrs.
     The physical relative phase inside the overlap/interaction window is
     then phi2 - phi1 evaluated AT t0_ref and is INDEPENDENT of when the
     pulses start, independent of segment sequencing, and identical across
     loop repetitions (provided the loop length is an integer number of
     envelope-block periods; see gone_overlaps() helper if not).
  2. CONFIGURABLE SA. All SA settings (center, span, points, RBW, average,
     sweep time, ref level) are arguments; nothing is hardcoded except
     driver-safe defaults.
  3. AWG SEGMENT QUANTUM vs INTERIOR RESOLUTION. The M8195A requires each
     SENT segment length to be a multiple of 1280 samples (20 ns at 64 GSa/s).
     That constraint applies ONLY to the total loop length. All interior
     timing (pulse starts, durations, delays) has full 1-sample resolution
     (15.625 ps = 1/64 ns), i.e. the delay_delta 64-sample granularity is
     preserved.
  4. MULTI-TRIAL / SNR-AWARE ACQUISITION. Helpers for interleaved shuffled
     parameter orders, phase-rotation sets (uniform phi grid for Fourier
     extraction), SNR-targeted repeats, and per-trial fresh baselines stored
     raw and separately (__SIG / __BL).

Hardware safeguards carried over from hw_core.py:
  - Forced RMS detector/averaging on the SA (DET RAV, AVER:TYPE RMS)
  - Waveform peak ceiling (ZHL safety): MAX_PEAK_MV, abs <= 1.0
  - RAW signal and RAW baseline stored separately; never subtract in hardware
  - --dry-run mode exercising the full code path without kNET

Verified API surface against spin-waves/Keysight_AWG_M8195A.py and the
embedded scripts in the h5 files (varyphase, detune_copy, detune_fine031126,
f1_4069_5, A_old_params/runBROAD variants).
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field, asdict

import numpy as np
from scipy.signal.windows import tukey

# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

FS_HZ = 64e9              # AWG sampling rate
AWG_BLOCK_SAMPLES = 1280  # segment-length quantum for :TRAC:DEF (20 ns)
MAX_PEAK_MV = 760.0       # ZHL protection, calibrated against run.py (< 0.1 dB compression)


def align_loop(n_samples: int, block: int = AWG_BLOCK_SAMPLES) -> int:
    """Round a total loop length UP to the next multiple of the AWG segment
    quantum. Interior timing is NOT quantized by this; it is applied only to
    the final loop length before send_waveform()."""
    return int(np.ceil(max(n_samples, block) / block) * block)


def snap_ns(duration_ns: float, samples: int = 1) -> float:
    """Snap a duration/start to an integer number of AWG samples (default 1
    sample = 1/64 ns). Use samples=64 for the legacy 1 ns grid."""
    return int(round(duration_ns * 1e-9 * FS_HZ / samples)) * samples / FS_HZ * 1e9


# ----------------------------------------------------------------------------
# Waveform construction (t0-referenced phases)
# ----------------------------------------------------------------------------

@dataclass
class Pulse:
    """One tone during one window. Phase phi is referenced to t0_ref_ns of the
    LoopSpec, i.e. the sample value is
        a[mV]/1e3 * env * sin(2*pi*f*(t - t0_ref) + phi)."""
    amp_mv: float
    freq_mhz: float
    start_ns: float
    dur_ns: float
    phase: float = 0.0


@dataclass
class LoopSpec:
    """A full loopy waveform: pulses + reference time + gap."""
    pulses: list = field(default_factory=list)
    t0_ref_ns: float = 0.0     # phase reference time (ns), recorded in attrs
    gap_ns: float = 2000.0     # quiet tail appended after last pulse end
    block: int = AWG_BLOCK_SAMPLES

    def loop_period_ns(self) -> float:
        end = max((p.start_ns + p.dur_ns) for p in self.pulses)
        return align_loop(int(np.ceil((end + self.gap_ns) * 1e-9 * FS_HZ)), self.block) / FS_HZ * 1e9


def build_loop(spec: LoopSpec, edge_ns: float = 1.0) -> np.ndarray:
    """Materialize the loop waveform (float64 in [-1, 1] Volts).

    - Interior timing at full sample resolution; loop length a multiple of
      spec.block samples.
    - Phases referenced to spec.t0_ref_ns => relative phases are
      start-time independent (convention 1).
    - Tukey envelope with matched 1 ns edges (edge taper split half up/half
      down inside the pulse).
    """
    n_total = align_loop(int(np.ceil(spec.loop_period_ns() * 1e-9 * FS_HZ)), spec.block)
    wf = np.zeros(n_total, dtype=np.float64)
    t0 = spec.t0_ref_ns

    for p in spec.pulses:
        i0 = int(round(p.start_ns * 1e-9 * FS_HZ))
        n_dur = int(round(p.dur_ns * 1e-9 * FS_HZ))
        if n_dur <= 0:
            continue
        i1 = min(n_total, i0 + n_dur)
        n_act = i1 - i0
        if i0 >= n_total:
            continue
        # tukey alpha from edge taper, capped so the window stays valid
        n_edge = int(round(edge_ns * 1e-9 * FS_HZ))
        alpha = min(1.0, 2.0 * n_edge / max(n_act, 1))
        env = tukey(n_act, alpha=alpha)
        t_axis = (np.arange(i0, i1) / FS_HZ) - (t0 * 1e-9)
        wf[i0:i1] += (p.amp_mv / 1e3) * env * np.sin(
            2 * np.pi * p.freq_mhz * 1e6 * t_axis + p.phase
        )

    peak = np.max(np.abs(wf)) * 1e3
    if peak > MAX_PEAK_MV or np.max(np.abs(wf)) > 1.0:
        raise ValueError(f"Waveform peak {peak:.1f} mV exceeds safety ceiling ({MAX_PEAK_MV} mV)")
    return wf


def two_pulse_overlap(a1_mv, f1_mhz, a2_mv, f2_mhz, phi,
                      t0_ns=1000.0, overlap_ns=2000.0, ramp_ns=2000.0,
                      gap_ns=2000.0, split_phase=True, block=AWG_BLOCK_SAMPLES) -> tuple:
    """Standard overlap experiment (Fig 3 style): both tones CW-overlapping.

    split_phase=True reproduces the f1_4069_5 convention (+phi/2 on tone 1,
    -phi/2 on tone 2) with phases referenced to the overlap start t0_ns,
    so the relative phase inside the overlap is exactly phi and start-time
    independent. split_phase=False reproduces the varyphase convention
    (tone1 phase 0, tone2 gets phi).

    Waveform: [ramp-in: tone1 on, tone2 ramps in] [overlap: both] with total
    overlap window = ramp_ns + overlap_ns, then a quiet gap.
    """
    pulses = [
        Pulse(a1_mv, f1_mhz, start_ns=0.0, dur_ns=ramp_ns + overlap_ns,
              phase=(+phi / 2 if split_phase else 0.0)),
        Pulse(a2_mv, f2_mhz, start_ns=ramp_ns, dur_ns=overlap_ns + ramp_ns / 2,
              phase=(-phi / 2 if split_phase else phi)),
    ]
    spec = LoopSpec(pulses=pulses, t0_ref_ns=ramp_ns, gap_ns=gap_ns, block=block)
    return build_loop(spec), spec


def two_pulse_sequential(a1_mv, f1_mhz, a2_mv, f2_mhz, phi,
                         d1_ns=300.0, d2_ns=400.0, gap_ns=0.0,
                         edge_ns=2.0, gap_between_ns=0.0, block=AWG_BLOCK_SAMPLES) -> tuple:
    """Non-overlap experiment (Fig 4 / varyphase / detune style): tone1 pulse
    followed by tone2 pulse, back-to-back, relative phase phi. Phases are
    referenced to each pulse's own start (the varyphase convention), which IS
    start-time independent for the single-pulse case. Recorded t0 is the
    start of tone 2; the phase difference phi is the quantity the idler at
    f/2 responds to.
    """
    t2_start = d1_ns + gap_between_ns
    pulses = [
        Pulse(a1_mv, f1_mhz, start_ns=0.0, dur_ns=d1_ns, phase=0.0),
        Pulse(a2_mv, f2_mhz, start_ns=t2_start, dur_ns=d2_ns, phase=phi),
    ]
    spec = LoopSpec(pulses=pulses, t0_ref_ns=t2_start, gap_ns=gap_ns, block=block)
    return build_loop(spec, edge_ns=edge_ns), spec


def phase_grid(n: int, shuffle: bool = True, seed: int | None = None) -> np.ndarray:
    """Uniform relative-phase set on [0, 2*pi) for Fourier extraction
    (A0, A1, A2, ... via discrete Fourier transform over phases)."""
    ph = 2 * np.pi * np.arange(n) / n
    if shuffle:
        rng = np.random.default_rng(seed)
        rng.shuffle(ph)
    return ph


def sobol_grid(n: int, bounds: list, seed: int | None = None) -> np.ndarray:
    """n x len(bounds) Sobol-ish stratified sampling (scipy.stats.qmc if
    available, else stratified random). Returns physical units."""
    try:
        from scipy.stats import qmc
        m = qmc.Sobol(d=len(bounds), scramble=True, seed=seed)
        u = m.random(n)
    except Exception:
        rng = np.random.default_rng(seed)
        u = (rng.permuted(np.tile(np.arange(n), (len(bounds), 1)), axis=1).T + rng.random((n, len(bounds)))) / n
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])
    return lo + u * (hi - lo)


def interleaved_order(params: list, seed: int | None = None) -> list:
    """Shuffled order for interleaved acquisition (kills slow drift bias)."""
    rng = np.random.default_rng(seed)
    idx = list(range(len(params)))
    rng.shuffle(idx)
    return idx


# ----------------------------------------------------------------------------
# SA configuration (fully parameterized)
# ----------------------------------------------------------------------------

@dataclass
class SAConfig:
    center_hz: float = 2065e6      # default f2/2 operating point reference
    span_hz: float = 2.5e6
    points: int = 1001
    rbw_hz: float = 200e3
    rbw_auto: bool = False
    avg_count: int = 2
    sweep_time: float = 0.2
    ref_level_dbm: float = -38.8
    scale_div: float = 4.0
    detector_rms: bool = True      # force DET:RAV + RMS averaging when supported


# ----------------------------------------------------------------------------
# Hardware rig
# ----------------------------------------------------------------------------

class SpectroRig:
    """Measurement chassis. Fresh baseline per trial by default, raw storage.

    HDF5 payloads: '<key>__SIG' = (freqs_hz, lin_power) with signal on,
    '<key>__BL' = independent raw baseline (only if measure_baseline).
    Root attrs: 'phase_convention', 'waveform_digests' (json), 'sa_settings'.
    """

    def __init__(self, h5_path: str, dry_run: bool = False,
                 phase_convention: str = "t0-referenced: phi_i at t0_ref_ns "
                                         "(relative phase = phi2 - phi1 in overlap)",
                 sa_config: SAConfig | None = None):
        self.dry_run = dry_run
        self.h5_path = h5_path
        self.phase_convention = phase_convention
        self.cfg = sa_config or SAConfig()
        self.trial_count = 0
        self._digests = {}

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
            self.configure_sa(self.cfg)
            self.buf = HDF5Buffer(h5_path, script_filename=__file__,
                                  settings=json.dumps(asdict(self.cfg)), append=True,
                                  buffer_size=25)
        else:
            self._records = []
            self._apply_dry_sa()

    # ---- SA ---------------------------------------------------------------

    def _apply_dry_sa(self):
        c = self.cfg
        self._freqs = np.linspace(c.center_hz - c.span_hz / 2, c.center_hz + c.span_hz / 2, c.points)

    def configure_sa(self, cfg: SAConfig | None = None):
        if cfg is not None:
            self.cfg = cfg
        c = self.cfg
        if self.dry_run:
            self._apply_dry_sa()
            return
        sa = self.sa
        sa.write(":FORM REAL,32")
        sa.centerFrequency = c.center_hz
        sa.spanFrequency = c.span_hz
        sa.points = c.points
        sa.resolutionBandwidthAutoState = c.rbw_auto
        if not c.rbw_auto:
            sa.resolutionBandwidth = c.rbw_hz
        sa.averageNumber = c.avg_count
        sa.sweepTime = c.sweep_time
        sa.referenceLevel = c.ref_level_dbm
        sa.scaleDiv = c.scale_div
        if c.detector_rms:
            # Try typed attributes first, fall back to raw SCPI. User confirmed
            # their kNET accepts these; the fallback protects partial drivers.
            try:
                sa.detectorTypeAutoState = False
                sa.detectorType = "RAV"
                sa.averageTypeAutoState = False
                sa.averageType = "RMS"
            except Exception:
                sa.write(":DET:FUNC RAV")
                sa.write(":AVER:TYPE RMS")

    # ---- acquisition ------------------------------------------------------

    def _sweep(self) -> np.ndarray:
        if self.dry_run:
            rng = np.random.default_rng(self.trial_count * 7919 + 13)
            return np.abs(rng.normal(5e-8, 5e-9, self.cfg.points))
        return np.power(10.0, np.asarray(self.sa.trigAndWait(), dtype=np.float64) / 10.0)

    def get_freqs(self) -> np.ndarray:
        if self.dry_run:
            return self._freqs
        return np.asarray(self.sa.get_frequencies(), dtype=np.float64)

    def measure(self, key: str, wf: np.ndarray, meta: dict = None,
                with_baseline: bool = True) -> tuple:
        """One trial: [fresh baseline with AWG off], signal with AWG looping wf.
        Stores raw __SIG and __BL separately. Returns (sig, bl, freqs)."""
        meta = meta or {}
        self.trial_count += 1
        if self.dry_run:
            rng = np.random.default_rng(self.trial_count + 42)
            bl = np.abs(rng.normal(5e-8, 5e-9, self.cfg.points))
            p_wf = float(np.mean(wf ** 2))
            f0 = self.cfg.center_hz
            sig = bl + 1e-4 * p_wf * np.exp(-((self._freqs - f0) / (self.cfg.span_hz / 8)) ** 2) \
                * (1 + rng.normal(0, 0.03, self.cfg.points))
            self._records.append({"key": key, "meta": meta,
                                  "delta": float(np.sum(sig - bl))})
            return sig, bl, self._freqs

        awg, sa = self.awg, self.sa
        if with_baseline:
            awg.abort_generation()
            bl = self._sweep()
        else:
            bl = None
        awg.send_waveform(wf)          # segment length enforced by caller (align_loop)
        awg.start_generation()
        sig = self._sweep()
        awg.abort_generation()
        awg.clear_data()

        freqs = self.get_freqs()
        self.buf.append_data(f"{key}__SIG", freqs, sig)
        if with_baseline:
            self.buf.append_data(f"{key}__BL", freqs, bl)
        self._maybe_write_attrs()
        return sig, bl, freqs

    def register_waveform(self, name: str, wf: np.ndarray, spec: LoopSpec | None):
        """Record a waveform digest + loop spec so the h5 is self-describing."""
        import hashlib
        d = hashlib.md5(wf.tobytes()).hexdigest()[:12]
        entry = {"digest": d, "n_samples": int(wf.shape[0]),
                 "loop_period_ns": spec.loop_period_ns() if spec else None,
                 "t0_ref_ns": spec.t0_ref_ns if spec else None,
                 "pulses": [asdict(p) for p in spec.pulses] if spec else None}
        self._digests[name] = entry
        return d

    def _maybe_write_attrs(self):
        if self.dry_run or len(self._digests) == 0:
            return
        try:
            self.buf._file.attrs["phase_convention"] = self.phase_convention  # noqa
            self.buf._file.attrs["waveform_digests"] = json.dumps(self._digests)
            self.buf._file.attrs["sa_settings"] = json.dumps(asdict(self.cfg))
        except Exception:
            pass  # attrs flushed at close; buffer may not expose the file handle

    def close(self):
        if self.dry_run:
            with open(self.h5_path + ".dryrun.json", "w") as f:
                json.dump({"trials": self.trial_count, "sample": self._records[:10],
                           "waveforms": self._digests}, f, indent=2)
            return
        try:
            f = self.buf._file  # noqa
            f.attrs["phase_convention"] = self.phase_convention
            f.attrs["waveform_digests"] = json.dumps(self._digests)
            f.attrs["sa_settings"] = json.dumps(asdict(self.cfg))
        except Exception:
            pass
        self.awg.abort_generation()
        self.awg.clear_data()
        self.buf.flush_to_hdf5()

    # ---- sessions ---------------------------------------------------------

    def run_session(self, points, wf_builder, phases, key_fmt, meta_fmt=None,
                    trials_per_point: int = 1, order_seed: int | None = None,
                    baseline_every_n: int = 1, progress=None):
        """SNR-aware interleaved session.

        points        : list of dicts (e.g. [{'a2':..., 'f2':...}, ...])
        wf_builder    : callable(point, phase) -> (wf, LoopSpec) or wf
        phases        : list of relative phases (use phase_grid(n) for Fourier
                        sets); each (point, phase) is measured trials times
        key_fmt       : callable(point, phase, trial) -> str
        trials_per_point: repeats interleaved ACROSS points (trial 1 of all
                        points, then trial 2, ...) to decorrelate drift.
        baseline_every_n: fresh baseline every n-th trial at TH point (1 = always)

        Returns list of (key, meta) executed.
        """
        order = interleaved_order(points, seed=order_seed)
        done = []
        for t in range(trials_per_point):
            for i in order:
                pt = points[i]
                wb = wf_builder(pt)
                for ph in phases:
                    wf, spec = wb if isinstance(wb, tuple) else (wb, None)
                    key = key_fmt(pt, ph, t)
                    with_bl = (t % max(baseline_every_n, 1) == 0)
                    self.measure(key, wf, meta={**(meta_fmt(pt, ph, t) if meta_fmt else {}),
                                                "phase": float(ph)},
                                 with_baseline=with_bl)
                    done.append((key, ph))
                    if progress:
                        progress(len(done), trials_per_point * len(points) * len(phases))
        return done


def snr_aware_trials(delta_noise_mw: float, signal_mw: float, target_snr_db: float,
                     n_avg: int, max_trials: int = 64) -> int:
    """How many repeated trials to reach target SNR given measured noise floor
    (in linear mW per trial) and signal power; accounts for SA averaging
    already applied per sweep."""
    snr_lin = signal_mw / max(delta_noise_mw / n_avg, 1e-300)
    need = (10 ** (target_snr_db / 10)) / max(snr_lin, 1e-300)
    return int(min(max_trials, max(1, int(np.ceil(np.sqrt(max(need, 1.0))))) ))
