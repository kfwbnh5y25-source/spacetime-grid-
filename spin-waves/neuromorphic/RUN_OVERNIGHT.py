"""
One-Click PyCharm Overnight Runner for Spin-Wave Neuromorphic Verification (E0 - E4).
Simply right-click -> Run in PyCharm.

Output: spin-waves/neuromorphic/data_e0_e4.h5 (+ .manifest.json)
Duration: ~6.2 hours at N_REP_SCALE = 2.0 (~4,640 trials).
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import run_e0_e4_verification as exp

# --- Overnight Configuration ---
OUT_FILE = str(HERE / "data_e0_e4.h5")
BLOCKS = ["E1", "E2", "E3", "E4"]
N_REP_SCALE = 2.0      # 2.0 -> ~6.2 h (~4,640 trials); 2.5 -> ~7.8 h (~5,800 trials)
PILOT_EVERY = 10
DRY_RUN = False

# Optional: path to another Python script to run automatically when finished (or None)
POST_RUN_SCRIPT = None  # e.g. r"C:\path\to\your_other_script.py"
POST_RUN_ARGS = []      # optional CLI args for the post-run script

if __name__ == "__main__":
    sys.argv = [
        __file__,
        "--out", OUT_FILE,
        "--blocks", *BLOCKS,
        "--n-rep-scale", str(N_REP_SCALE),
        "--pilot-every", str(PILOT_EVERY),
    ] + (["--dry-run"] if DRY_RUN else [])
    exp.main()

    if POST_RUN_SCRIPT:
        print(f"[POST-RUN] Launching {POST_RUN_SCRIPT}...")
        subprocess.run([sys.executable, str(POST_RUN_SCRIPT), *POST_RUN_ARGS], check=False)
