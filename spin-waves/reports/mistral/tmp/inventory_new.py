"""Inventory of the two new files added today (2026-11-17):
phase gating/delay_f1_JULY_varyf1a1.h5  (428 MB)
phase gating/fig1.zip                   (876 MB)
Prints key structures, parameter grids, and (for the h5) a quick idler-line
harmonic sanity value for a few keys.
"""
import io
import re
import json
import zipfile
import numpy as np
import h5py

BASE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating"
OUT = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp"

print("========== delay_f1_JULY_varyf1a1.h5 ==========")
with h5py.File(f"{BASE}/delay_f1_JULY_varyf1a1.h5", "r") as f:
    attrs = dict(f.attrs)
    sc = attrs.get("script_content", "")
    print("root attrs:", [k for k in attrs])
    if sc:
        print("---- embedded script (first 2600 chars) ----")
        print(sc[:2600])
    keys = sorted(f.keys())
    print(f"\n{n_keys} keys" if (n_keys := len(keys)) else "")
    print("first keys:", keys[:5])
    print("last keys:", keys[-5:])
    # parse
    rec = {}
    ok = 0
    for k in keys:
        pr = k.split("_")
        if len(pr) != 4:
            continue
        rec.setdefault((pr[0], pr[1]), set()).add(pr[2])
        ok += 1
    delays = sorted(set(float(a) for a, _ in rec))
    f1s = sorted(set(float(b) for _, b in rec))
    nph = sorted(set(len(v) for v in rec.values()))
    print(f"parsed {ok}/{len(keys)}; delay level count={len(delays)}, f1 level count={len(f1s)}, phases-per-combo={nph[:3]}")
    if delays:
        print(f"delay range: {min(delays)} .. {max(delays)}")
        print("delays:", np.round(delays, 3).tolist())
    if f1s:
        print(f"f1 range: {min(f1s)} .. {max(f1s)}")
    if keys:
        g = f[keys[0]]
        print("example group members:", list(g.keys()))
        if "trace" in g:
            print("trace shape/dtype:", g["trace"].shape, g["trace"].dtype)

print("\n========== fig1.zip ==========")
with zipfile.ZipFile(f"{BASE}/fig1.zip") as z:
    names = z.namelist()
    print(f"{len(names)} members")
    for nm in names[:40]:
        info = z.getinfo(nm)
        print(f"  {nm}  {info.file_size/1e6:.1f} MB" if info.file_size > 1e5 else f"  {nm}  {info.file_size} B")
