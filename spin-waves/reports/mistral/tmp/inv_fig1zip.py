import h5py, numpy as np
BASE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/fig1_extract"
for name in ("delay.h5", "delay2.h5"):
    with h5py.File(f"{BASE}/{name}", "r") as f:
        keys = sorted(f.keys())
        sc = f.attrs.get("script_content", "")
        print(f"===== {name}: {len(keys)} keys =====")
        if sc:
            # print only the sweep-loop part of the script
            i = sc.find("def main")
            j = sc.find("sa.centerFrequency", i)
            print(sc[i:i+1400])
        pr = [k.split("_") for k in keys[:2000]]
        lens = set(len(p) for p in pr)
        print("field counts:", lens)
        if 4 in lens:
            fields = {}
            for idx in range(4):
                vals = sorted(set(p[idx] for p in pr if len(p) == 4))
                fields[idx] = (vals[:6], vals[-3:], len(vals))
            for idx, (lo, hi, n) in fields.items():
                print(f"  field {idx}: n={n} first={lo} last={hi}")
        if keys:
            g = f[keys[0]]
            mem = list(g.keys())
            print("  group members:", mem)
            if "trace" in g:
                print("  trace shape:", g["trace"].shape)
