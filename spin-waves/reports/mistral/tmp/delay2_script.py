"""Extract script_content + waveform/loop structure from delay.h5 / delay2.h5."""
import h5py
import numpy as np

BASE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/fig1_extract/"

for fn in ("delay2.h5", "delay.h5"):
    print("=" * 70)
    print(fn)
    with h5py.File(BASE + fn, "r") as f:
        root = f.attrs.get("script_content") or f.attrs.get("script")
        if root is None:
            # search groups for attrs
            for k in list(f.keys())[:3]:
                g = f[k]
                for ak in g.attrs:
                    print(" group", k, "attr:", ak)
            root = None
        if root:
            sc = root.decode() if isinstance(root, bytes) else str(root)
            print(f"--- script_content: {len(sc)} chars; AWG/waveform lines:")
            for line in sc.splitlines():
                s = line.strip()
                if any(w in s for w in ("send_waveform", "np.tile", "resize",
                                        "for ", "arange", "zeros", "segment",
                                        "1280", "sin(", "where", "while", "sleep",
                                        "baseline", "avg", "AVER", "SWEEP")):
                    print("   ", s[:120])
        # grab one dataset: waveform info embedded?
        key0 = list(f.keys())[0]
        g = f[key0]
        print("first key:", key0, "members:", list(g.keys()))
        for ak, av in g.attrs.items():
            s = str(av)
            print("   attr", ak, "=", s[:150])
