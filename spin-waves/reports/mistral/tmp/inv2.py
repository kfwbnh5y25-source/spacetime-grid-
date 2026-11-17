import h5py
import re

for name in ("delay.h5", "delay2.h5"):
    with h5py.File(f"fig1_extract/{name}", "r") as f:
        sc = f.attrs.get("script_content", "")
        print(f"===== {name} =====")
        for m in re.finditer(r"^(amps|freqs|sigma)\s*=.*$", sc, re.M):
            print("  ", m.group(0))
        i = sc.find("def waveform")
        print(sc[i:i + 560].replace("\r\n", "\n"))
        j = sc.find("def measure")
        print(sc[j:j + 380].replace("\r\n", "\n"))
