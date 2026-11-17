import h5py, numpy as np, collections, re

p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/detune_fine031126.h5"
pairs_by_phase = collections.defaultdict(set)
phases_by_pair = collections.defaultdict(set)
with h5py.File(p, "r") as f:
    for k in f.keys():
        m = re.match(r"([\d.]+)_([\d.]+)_([\d.]+)_(\S+)", k)
        f1, f2, ph = float(m.group(1)), float(m.group(2)), float(m.group(3))
        pairs_by_phase[ph].add((f1, f2))
        phases_by_pair[(f1, f2)].add(ph)
print("total pairs:", len(phases_by_pair))
print("pairs per phase:", {ph: len(s) for ph, s in sorted(pairs_by_phase.items())})
nph = np.array([len(v) for v in phases_by_pair.values()])
print("phases-per-pair histogram:", collections.Counter(nph))
print("f1 range", sorted(set(f1 for (f1, f2) in phases_by_pair))[:5], "...",
      sorted(set(f1 for (f1, f2) in phases_by_pair))[-3:])
print("f2 range", sorted(set(f2 for (f1, f2) in phases_by_pair))[:5], "...",
      sorted(set(f2 for (f1, f2) in phases_by_pair))[-3:])
