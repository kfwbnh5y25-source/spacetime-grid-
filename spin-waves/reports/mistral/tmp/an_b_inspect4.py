import h5py, numpy as np
fn = 'phase gating/delay_f1_JULY_varyf1a1.h5'
with h5py.File(fn,'r') as f:
    sc = f.attrs['script_content']
    print(sc[1450:1700])
    print('=== settings ===')
    print(f.attrs['settings'])
    g = f['256.0000_4020_0.000_1']
    tr = g['trace'][0]
    print('fields[0][:6]:', tr[0][:6], '...[6]:', tr[0][-3:])
    print('fields[1][:6]:', tr[1][:6])
    ts = g['timestamp'][0,0]
    print('ts', ts)
    ks = list(f.keys())
    print('ts last', f[ks[-1]]['timestamp'][0,0])
