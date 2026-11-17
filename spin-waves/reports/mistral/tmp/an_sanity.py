"""Analyze the real meas0_sanity.h5 against the delay2.h5 anchors.

Storage format = hw_core_spectro HDF5Buffer: '<key>__SIG' and '<key>__BL'
datasets with (freqs, lin_power). We re-derive the sanity verdicts from raw
spectra so a discrepancy with the live verdicts becomes debuggable.
"""
import numpy as np
from scipy.signal import medfilt

path = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/measure/meas0_sanity.h5"
F2 = 4149.0
CENTER = F2 / 2   # SA center on the sample path, MHz

import h5py

f = h5py.File(path, "r")
print("attrs:", dict(f.attrs))
keys = list(f.keys())
print(f"{len(keys)} datasets")
sigs = sorted(k[:-5] for k in keys if k.endswith("__SIG"))
print(sigs[:6], "...", sigs[-3:])


def line_power(sig, bl, freqs):
    d = sig - bl
    d = medfilt(d, 9)
    win = np.abs(freqs / 1e6 - CENTER) <= 1.5
    noise_bins = ~win & (d < np.quantile(d, 0.9))
    return float(d[win].sum()), float(np.std(d[noise_bins]))


def parse(key):
    # SAN_d{+2.00}_p{1.2345} / SAN_OVL_t1_p0.0000 / SAN_REF2_t0
    d = np.nan
    ph = np.nan
    t = "t0"
    for tok in key.split("_"):
        if tok.startswith("d"):
            try:
                d = float(tok[1:])
            except ValueError:
                pass
        elif tok.startswith("p"):
            try:
                ph = float(tok[1:])
            except ValueError:
                pass
        elif tok.startswith("t") and any(c.isdigit() for c in tok):
            try:
                tt = float(tok[1:])
                if tok[0] == "t" and tt in (0.0, 1.0):
                    t = tok
                    if "OVL" in key:
                        d = 0.0
            except ValueError:
                pass
    return d, ph


# group by block
blocks = {}
refs = {}
for k in sigs:
    fr = np.array(f[k + "__SIG"]["trace"][:, 0])
    sig = np.array(f[k + "__SIG"]["trace"][:, 1])
    if k + "__BL" in f:
        bl = np.array(f[k + "__BL"]["trace"][:, 1])
    else:
        bl = np.zeros_like(sig)
    p, n = line_power(sig, bl, fr)
    if "REF2" in k:
        refs[k] = (p, n)
        continue
    d, ph = parse(k)
    blocks.setdefault(d, []).append((ph, p, n, k))

ref_p = np.mean([p for p, _ in refs.values()]) if refs else np.nan
print(f"\nREF2: mean line = {ref_p*1e6:.2f} uW, noise = "
      f"{np.mean([n for _, n in refs.values()])*1e6:.2f} uW"
      f"  ({len(refs)} refs)")
print(f"{'t_gap':>8s} {'n':>3s} {'A0_uW':>8s} {'A1_uW':>8s} {'A1/A0':>6s} "
      f"{'phi0':>7s} {'max/ref':>8s}")


def c12(rows):
    ph = np.array([r[0] for r in rows])
    pw = np.array([r[1] for r in rows])
    if len(rows) < 8 or np.ptp(ph) < 1:
        return None
    c1 = np.mean(pw * np.exp(-1j * ph))
    return float(2 * np.abs(c1)), float(np.mean(pw)), float(np.angle(c1))


for d in sorted(blocks):
    rows = blocks[d]
    cc = c12(rows)
    xref = max(r[1] for r in rows) / ref_p if refs else np.nan
    if cc:
        A1, A0, phi0 = cc
        print(f"{d:+8.2f} {len(rows):3d} {A0*1e6:8.1f} {A1*1e6:8.1f} {A1/A0:6.2f} "
              f"{phi0:+7.2f} {xref:8.1f}")
    else:
        print(f"{d:+8.2f} {len(rows):3d}  (fewer than 8 phases or degenerate grid)")

# spectrum shape sanity: where is the line actually?
k0 = [k for k in sigs if "OVL" in k][0]
sig = np.array(f[k0 + "__SIG"]["trace"][:, 1])
bl = np.array(f[k0 + "__BL"]["trace"][:, 1]) if k0 + "__BL" in f else np.zeros_like(sig)
fr = np.array(f[k0 + "__SIG"]["trace"][:, 0])
d = medfilt(sig - bl, 9)
top = np.argsort(d)[-5:]
print(f"\nexample '{k0}': center {fr[500]/1e6:.2f} MHz, span {np.ptp(fr)/1e6:.1f} MHz")
print("top-5 bins (f_MHz, d_mW):", [(round(fr[i]/1e6, 3), round(d[i]*1e3, 4)) for i in top])
f.close()
