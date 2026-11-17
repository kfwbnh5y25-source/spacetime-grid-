import h5py, numpy as np, re
from collections import Counter

p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/detune_fine.h5"
with h5py.File(p, 'r') as f:
    keys = list(f.keys())
    print("n_keys", len(keys))
    print("samples", keys[:3], keys[-2:])
    pat = re.compile(r'(?P<f1>[-\d.]+)_(?P<f2>[-\d.]+)_(?P<phase>[-\d.]+)_(?P<label>\w+)')
    ok, bad = [], []
    for k in keys:
        m = pat.match(k)
        (ok.append(m.groupdict()) if m else bad.append(k))
    print("bad", len(bad), bad[:3])
    f1s = sorted(set(float(d['f1']) for d in ok))
    f2s = sorted(set(float(d['f2']) for d in ok))
    phs = sorted(set(float(d['phase']) for d in ok))
    labs = Counter(d['label'] for d in ok)
    print("f1", f1s[0], "-", f1s[-1], len(f1s))
    print("f2", f2s[0], "-", f2s[-1], len(f2s))
    print("phases", len(phs), phs[:6], "...", phs[-3:])
    print("labels", labs)
    import json
    sc = f.attrs.get('script_content')
    if sc:
        sc = sc.decode() if isinstance(sc, bytes) else str(sc)
        open("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/detune_fine_script.py","w").write(sc)
        print("script saved", len(sc))
    st = f.attrs.get('settings')
    print("settings:", (st[:300] if isinstance(st, (str, bytes)) else st))
    k0 = ok and keys[0]
    tr = f[k0]['trace']
    print("trace shape", tr.shape)
