import h5py, os, re, glob

base = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating"
out_dir = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/scripts_from_h5"
os.makedirs(out_dir, exist_ok=True)

files = [
    os.path.join(base, "varyphase.h5"),
    os.path.join(base, "detune copy.h5"),
    os.path.join(base, "detune_fine031126.h5"),
]

def sanitize(name):
    name = re.sub(r'[^A-Za-z0-9._-]+', '_', name)
    return name[:80]

def extract(path, tag):
    try:
        with h5py.File(path, 'r') as f:
            sc = f.attrs.get('script_content')
            if sc is None:
                print(tag, "NO script_content")
                return
            if isinstance(sc, bytes):
                sc = sc.decode(errors='replace')
            fn = os.path.join(out_dir, sanitize(tag + ".py"))
            with open(fn, "w") as out:
                out.write(sc)
            print("wrote", fn, len(sc))
    except Exception as e:
        print(tag, "ERR", e)

for p in files:
    extract(p, os.path.basename(p))

# delay_delta.zip contents: extract h5 members to tmp first is big; instead try direct reading of zipped h5 via h5py is not possible
# unzip the two scripts region: use zipfile to list, then extract only needed
import zipfile
zp = os.path.join(base, "delay_delta.zip")
if os.path.exists(zp):
    with zipfile.ZipFile(zp) as z:
        print("zip members:", z.namelist())

print("--- done loose files ---")
