"""
V2: Extracts Log-Decibel Dual-Formant Envelopes from Free Spoken Digit Dataset:
  - 4 classes: Digits '0', '1', '2', '3' ("zero", "one", "two", "three")
  - 4 speakers: 'jackson', 'nicolas', 'theo', 'yweweler' (80 audio files)
  - Uses independent log-mel dB normalization (30 dB dynamic range floor)
    so high-frequency sibilants/consonants are not drowned out by low vowel fundamentals!
  - Saves to: spin-waves/neuromorphic/spoken_digits_features_v2.json
"""

import glob
import json
import os
import wave
from pathlib import Path
import numpy as np
from scipy.signal import spectrogram

HERE = Path(__file__).resolve().parent
AUDIO_DIR = HERE / "audio_cache"
OUT_FEATURES = HERE / "spoken_digits_features_v2.json"

DIGITS = [0, 1, 2, 3]
SPEAKERS = ["jackson", "nicolas", "theo", "yweweler"]


def extract_log_features(wav_path: str, n_steps: int = 5):
    """
    Computes log-magnitude spectrogram and independent dB-normalized formant envelopes:
      - Low Formant (200-1200 Hz): vowel/vocal fundamental energy
      - High Formant (1200-4000 Hz): consonant/sibilant/plosive friction energy
    """
    with wave.open(wav_path, "rb") as wf:
        n_channels = wf.getnchannels()
        sampwidth = wf.getsampwidth()
        fs = wf.getframerate()
        n_frames = wf.getnframes()
        raw = wf.readframes(n_frames)

    dtype = np.int16 if sampwidth == 2 else np.uint8
    audio = np.frombuffer(raw, dtype=dtype).astype(float)
    if n_channels > 1:
        audio = audio.reshape(-1, n_channels).mean(axis=1)

    audio = audio / (np.max(np.abs(audio)) + 1e-12)

    # Trim leading/trailing silence (2% energy threshold)
    win = int(0.01 * fs)
    if win > 0 and len(audio) > win:
        rolling = np.convolve(audio**2, np.ones(win)/win, mode="same")
        thresh = 0.02 * np.max(rolling)
        active = np.where(rolling > thresh)[0]
        if len(active) > 200:
            audio = audio[active[0]:active[-1]]

    # 25 ms windowed spectrogram
    nperseg = min(len(audio), int(0.025 * fs))
    if nperseg < 16:
        nperseg = 16
    f_ax, t_ax, Sxx = spectrogram(audio, fs=fs, nperseg=nperseg, noverlap=nperseg // 2)

    mask_lo = (f_ax >= 200) & (f_ax <= 1200)
    mask_hi = (f_ax > 1200) & (f_ax <= 4000)

    e_lo = np.sum(Sxx[mask_lo, :], axis=0) if np.any(mask_lo) else np.ones(Sxx.shape[1])
    e_hi = np.sum(Sxx[mask_hi, :], axis=0) if np.any(mask_hi) else np.ones(Sxx.shape[1])

    # Log-energy in decibels
    log_lo = 10.0 * np.log10(e_lo + 1e-10)
    log_hi = 10.0 * np.log10(e_hi + 1e-10)

    # Independent normalization with 30 dB floor
    def norm_db(v):
        floor = np.max(v) - 30.0
        v_clip = np.clip(v, floor, np.max(v))
        return (v_clip - np.min(v_clip)) / (np.max(v_clip) - np.min(v_clip) + 1e-6)

    n_lo = norm_db(log_lo)
    n_hi = norm_db(log_hi)

    # Segment into n_steps equal time frames
    splits = np.array_split(np.arange(len(t_ax)), n_steps)
    step_lo = [float(np.mean(n_lo[idx])) for idx in splits]
    step_hi = [float(np.mean(n_hi[idx])) for idx in splits]

    return step_lo, step_hi


def main():
    wav_files = sorted(glob.glob(str(AUDIO_DIR / "*.wav")))
    if not wav_files:
        print(f"[FAIL] No audio files found in {AUDIO_DIR}")
        return

    dataset = []
    for f in wav_files:
        base = os.path.basename(f)
        parts = base.replace(".wav", "").split("_")
        if len(parts) < 3:
            continue
        digit, spk, rep = int(parts[0]), parts[1], int(parts[2])
        if digit not in DIGITS or spk not in SPEAKERS:
            continue

        step_lo, step_hi = extract_log_features(f, n_steps=5)
        dataset.append({
            "digit": digit,
            "speaker": spk,
            "rep": rep,
            "low_formant": step_lo,
            "high_formant": step_hi,
        })

    OUT_FEATURES.write_text(json.dumps(dataset, indent=2))
    print("=" * 68)
    print(f"[SUCCESS] Extracted V2 log-decibel features for {len(dataset)} audio files.")
    print(f"Saved to: {OUT_FEATURES}")
    digits = [d["digit"] for d in dataset]
    print("Class counts:", {d: digits.count(d) for d in DIGITS})
    print("=" * 68)


if __name__ == "__main__":
    main()
