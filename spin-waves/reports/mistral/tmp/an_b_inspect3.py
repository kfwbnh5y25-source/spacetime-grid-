import h5py, numpy as np
fn = 'phase gating/delay_f1_JULY_varyf1a1.h5'
with h5py.File(fn,'r') as f:
    sc = f.attrs['script_content']
    print(sc[1900:2950])
