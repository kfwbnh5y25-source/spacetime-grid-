"""Full script_content + fast line-position scan across delay2.h5 (June drift)."""
import h5py
import numpy as np

BASE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/fig1_extract/"

with h5py.File(BASE + "delay2.h5", "r") as f:
    root = f.attrs["script_content"]
    sc = root.decode() if isinstance(root, bytes) else str(root)
    print("===== delay2.h5 script_content (full) =====")
    print(sc)
