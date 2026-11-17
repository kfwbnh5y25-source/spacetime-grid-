# Salvage truncated delay_fine.h5 via superblock EOF patch (tmp only).
import struct, os
import h5py

fn = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/salvage_delay_fine.h5"
size = os.path.getsize(fn)
with open(fn, "r+b") as fh:
    fh.seek(0); sig = fh.read(8)
    print("sig ok:", sig == b"\x89HDF\r\n\x1a\n")
    fh.seek(8); v = fh.read(8)
    print("superblock v:", v[0], "offset/length sizes:", v[5], v[6])
    fh.seek(40); stored_eof = struct.unpack("<Q", fh.read(8))[0]
    print("stored eof:", stored_eof, "-> patching to", size)
    fh.seek(40); fh.write(struct.pack("<Q", size))

f = h5py.File(fn, "r")
keys = sorted(f.keys())
print("n_keys:", len(keys))
print("first 10:", keys[:10])
print("last 10:", keys[-10:])
from collections import defaultdict
nf = defaultdict(int)
for k in keys:
    nf[len(k.split("_"))] += 1
print("field counts:", dict(nf))
print("attrs:", list(f.attrs.keys()))
tr = f[keys[10]]["trace"]
print("trace shape:", tr.shape, tr.dtype)
print("ch0:", tr[0, 0, 0], tr[0, 0, 500], tr[0, 0, -1])
sc = f.attrs.get("script_content")
if sc is not None:
    sc = sc.decode() if isinstance(sc, bytes) else str(sc)
    open("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/delay_fine_script.py", "w").write(sc)
    print("script saved, len", len(sc))
