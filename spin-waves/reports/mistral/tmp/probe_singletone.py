import h5py, re, collections
for p in ["/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/Figure!/Figure1_arnold_plot/data.h5",
          "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/Figure!/arnold_tongue_zoom/data2.h5"]:
    with h5py.File(p, "r") as f:
        keys = list(f.keys())
    KEY_RE = re.compile(r"^(\d+)_(\d+)_(half|full|halfaverage)_(\d{8})_(\d{6})_(\d+)_(\d+)$")
    ch = collections.Counter(KEY_RE.match(k).group(3) for k in keys if KEY_RE.match(k))
    print(p.split("/")[-2], "keys:", len(keys), "matched:", dict(ch),
          "unmatched:", len(keys) - sum(ch.values()))
