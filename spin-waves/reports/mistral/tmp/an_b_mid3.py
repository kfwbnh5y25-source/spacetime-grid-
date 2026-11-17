import h5py
with h5py.File('phase gating/delay_f1_JULY_varyf1a1.h5','r') as f:
    sc = f.attrs['script_content']
    print(sc[3560:4470])
