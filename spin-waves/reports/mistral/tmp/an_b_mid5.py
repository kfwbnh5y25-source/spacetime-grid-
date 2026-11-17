import h5py, re
with h5py.File('phase gating/delay_f1_JULY_varyf1a1.h5','r') as f:
    sc = f.attrs['script_content']
    for m in re.finditer(r'.*(a1s|f1s|phases)\s*=.*', sc):
        print(m.group(0).strip())
    i = sc.find('def main')
    print(sc[i:i+900])
