"""
Fetches 80 authentic spoken digit audio recordings from the Free Spoken Digit Dataset (FSDD):
  - 4 classes: Digits '0', '1', '2', '3' ("zero", "one", "two", "three")
  - 4 speakers: 'jackson', 'nicolas', 'theo', 'yweweler'
  - 5 utterances per speaker per digit = 20 trials per class = 80 total trials
Extracts 5-step temporal dual-formant energy envelopes:
  - Low Formant (200-1200 Hz): mapped to Tone A (4081 MHz)
  - High Formant (1200-4000 Hz): mapped to Tone B (4125 MHz)
Saves:
  - spin-waves/neuromorphic/audio_cache/ (cached WAVs)
  - spin-waves/neuromorphic/spoken_digits_features.json
"""

import json
import urllib.request
import wave
from pathlib import Path
import numpy as np
from scipy.signal import spectrogram

HERE = Path(__file__).resolve().parent
AUDIO_DIR = HERE / "audio_cache"
AUDIO_DIR.mkdir(parents=True, exist_ok=True)

DIGITS = [0, 1, 2, 3]
SPEAKERS = ["jackson", "nicolas", "theo", "yweweler"]
REPS_PER_SPEAKER = 5  # 5 reps * 4 speakers = 20 per digit = 80 total trials

GITHUB_BASE = "https://raw.githubusercontent.com/Jakobovski/free-spoken-digit-dataset/master/recordings"


def download_audio_files():
    manifest = []
    print("Fetching authentic human audio recordings from FSDD repository...")
    for digit in DIGITS:
        for spk in SPEAKERS:
            for rep in range(REPS_PER_SPEAKER):
                fname = f"{digit}_{spk}_{rep}.wav"
                fpath = AUDIO_DIR / fname
                if not fpath.exists():
                    url = f"{GITHUB_BASE}/{fname}"
                    try:
                        with urllib.request.urlopen(url, timeout=10) as resp:
                            fpath.write_bytes(resp.read())
                    except Exception as e:
                        print(f"Error fetching {url}: {e}")
                        continue
                manifest.append({"digit": digit, "speaker": spk, "rep": rep, "path": str(fpath)})
    print(f"Total audio files ready: {len(manifest)}")
    return manifest


def extract_dual_formant_envelope(wav_path: str, n_steps: int = 5):
    """
    Reads a WAV file and computes the normalized 5-step dual-formant acoustic energy:
      - Low Formant (200-1200 Hz): vowel/vocal cord fundamental energy
      - High Formant (1200-4000 Hz): consonant/sibilant/fricative energy
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

    # Normalize audio
    audio = audio / (np.max(np.abs(audio)) + 1e-12)

    # Trim leading/trailing silence (energy threshold 5%)
    win = int(0.01 * fs)
    if win > 0 and len(audio) > win:
        rolling_sq = np.convolve(audio**2, np.ones(win)/win, mode="same")
        thresh = 0.03 * np.max(rolling_sq)
        active = np.where(rolling_sq > thresh)[0]
        if len(active) > 200:
            audio = audio[active[0]:active[-1]]

    # Compute Spectrogram
    nperseg = min(len(audio), int(0.025 * fs))  # 25 ms window
    if nperseg < 16:
        nperseg = 16
    f_axis, t_axis, Sxx = spectrogram(audio, fs=fs, nperseg=nperseg, noverlap=nperseg // 2)

    # Define formant masks
    mask_low = (f_axis >= 200) & (f_axis <= 1200)
    mask_high = (f_axis > 1200) & (f_axis <= 4000)

    e_low = np.sum(Sxx[mask_low, :], axis=0) if np.any(mask_low) else np.ones(Sxx.shape[1])
    e_high = np.sum(Sxx[mask_high, :], axis=0) if np.any(mask_high) else np.ones(Sxx.shape[1])

    # Downsample / average into n_steps equal temporal segments
    idx_splits = np.array_split(np.arange(len(t_axis)), n_steps)
    low_envelope = np.array([float(np.mean(e_low[idx])) if len(idx) > 0 else 0.0 for idx in idx_splits])
    high_envelope = np.array([float(np.mean(e_high[idx])) if len(idx) > 0 else 0.0 for idx in idx_splits])

    # Normalize each to [0, 1] relative to peak
    max_e = max(np.max(low_envelope), np.max(high_envelope), 1e-12)
    low_norm = (low_envelope / max_e).tolist()
    high_norm = (high_envelope / max_e).tolist()

    return low_norm, high_norm


def main():
    manifest = download_audio_files()
    dataset = []
    for item in manifest:
        try:
            low_env, high_env = extract_dual_formant_envelope(item["path"], n_steps=5)
            dataset.append({
                "digit": item["digit"],
                "speaker": item["speaker"],
                "rep": item["rep"],
                "low_formant": low_env,
                "high_formant": high_env,
            })
        except Exception as e:
            print(f"Error processing {item['path']}: {e}")

    out_file = HERE / "spoken_digits_features.json"
    out_file.write_text(json.dumps(dataset, indent=2))
    print(f"\nExtracted dual-formant envelopes for {len(dataset)} trials.")
    print(f"Saved to {out_file}")

    # Inspect class balance
    digits = [d["digit"] for d in dataset]
    print("Class distribution:", {d: digits.count(d) for d in DIGITS})


if __name__ == "__main__":
    main()
