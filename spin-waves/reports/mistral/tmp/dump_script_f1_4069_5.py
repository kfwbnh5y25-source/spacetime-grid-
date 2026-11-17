import h5py

p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
f = h5py.File(p, 'r')
sc = f.attrs['script_content']
if isinstance(sc, bytes):
    sc = sc.decode()
with open("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/f1_4069_5_script.py", "w") as out:
    out.write(sc)
print("script len", len(sc))
print("settings type", type(f.attrs['settings']))
st = f.attrs['settings']
if isinstance(st, bytes):
    st = st.decode()
print(str(st)[:2000])
f.close()
