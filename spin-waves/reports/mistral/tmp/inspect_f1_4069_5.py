import h5py, numpy as np, re
from collections import Counter

p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
f = h5py.File(p, 'r')
keys = list(f.keys())
print("n_keys", len(keys))
print("sample", keys[:3], keys[-2:])

pat = re.compile(r'(?P<a1>[\d.]+)_(?P<a2>[\d.]+)_(?P<f1>[\d.]+)_(?P<f2>[\d.]+)_(?P<phase>[-\d.]+)_(?P<label>\w+)')
parsed = []
bad = []
for k in keys:
    m = pat.match(k)
    if m:
        parsed.append(m.groupdict())
    else:
        bad.append(k)
print("bad", len(bad), bad[:5])

def col(name):
    return np.array([float(d[name]) for d in parsed])

a1s = sorted(set(col('a1').tolist()))
a2s = sorted(set(col('a2').tolist()))
f1s = sorted(set(col('f1').tolist()))
f2s = sorted(set(col('f2').tolist()))
phs = sorted(set(col('phase').tolist()))
labs = Counter(d['label'] for d in parsed)
print("a1", a1s)
print("a2", a2s)
print("f1", f1s)
print("f2 n=%d range %.1f..%.1f step~%.2f" % (len(f2s), f2s[0], f2s[-1], (f2s[-1]-f2s[0])/max(len(f2s)-1,1)))
print("phases", phs)
print("labels", labs)

k0 = keys[0]
ds = f[k0]
print("entry type:", type(ds))
if hasattr(ds, 'keys'):
    for kk in list(ds.keys())[:6]:
        sub = ds[kk]
        print(" ", kk, getattr(sub, 'shape', None), getattr(sub, 'dtype', None))
    if 'script_content' in ds.attrs or 'script' in dict(ds.attrs):
        pass
print("attrs of entry:", {k: (v[:80] if isinstance(v, bytes) else str(v)[:80]) for k, v in ds.attrs.items()} if ds.attrs else {})
print("root attrs keys:", list(f.attrs.keys()))
f.close()
