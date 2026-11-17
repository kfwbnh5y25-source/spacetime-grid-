import h5py, re, collections, json
path = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
with h5py.File(path,'r') as f:
    keys = sorted(f.keys())
print("n keys:", len(keys))
for k in keys[:5]: print(k)
print("...")
for k in keys[-5:]: print(k)
pat = re.compile(r'^(\d+)_(\d+)_(\d+)_(\d+)_([\d.]+)_(.+)$')
groups = collections.Counter(); phases = collections.Counter(); f1s = collections.Counter(); a2s=collections.Counter(); f2s=collections.Counter(); labels=collections.Counter()
bad=[]
for k in keys:
    m = pat.match(k)
    if not m: bad.append(k); continue
    a1,a2,f1,f2,ph,lab = m.groups()
    groups[int(a1)] += 1; phases[float(ph)] += 1; f1s[int(f1)] += 1
    a2s[int(a2)] += 1; f2s[int(f2)] += 1; labels[lab] += 1
print("bad:", bad[:10], len(bad))
print("a1 groups:", sorted(groups.items()))
print("phases:", dict(phases))
print("f1s:", dict(sorted(f1s.items())[:40]))
print("a2s:", sorted(a2s.items()))
print("f2s:", sorted(f2s.items()))
print("labels:", dict(labels))
with h5py.File(path,'r') as f:
    k0 = '0_0_0_0' if '0_0_0_0' in keys else keys[0]
    for k in keys[:3]:
        obj = f[k]
        print(k, type(obj))
        try: print(" attrs:", dict(obj.attrs))
        except Exception as e: print(" attrs err", e)
        try: print(" keys:", list(obj.keys()), "trace shape:", obj['trace'].shape, obj['trace'].dtype)
        except Exception as e: print(" trace err:", e)
