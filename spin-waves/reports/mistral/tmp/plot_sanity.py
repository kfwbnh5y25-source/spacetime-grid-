"""Dump + plot sanity-check traces: raw SIG, raw BL, difference spectra."""
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

path = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/measure/meas0_sanity.h5"
OUT = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/"

f = h5py.File(path, "r")
keys = sorted(k[:-5] for k in f.keys() if k.endswith("__SIG"))


def read(key):
    tr = np.array(f[key + "__SIG"]["trace"])          # (rep, 2, 1001)
    freqs = np.mean(tr[:, 0, :], axis=0)
    sig = np.mean(tr[:, 1, :], axis=0).astype(float)
    if key + "__BL" in f:
        trb = np.array(f[key + "__BL"]["trace"])
        bl = np.mean(trb[:, 1, :], axis=0).astype(float)
    else:
        bl = np.zeros_like(sig)
    return freqs, sig, bl


sel = {"ref": [k for k in keys if "REF2" in k][:2],
       "ovl": [k for k in keys if "OVL" in k],
       "delays": [k for k in keys if k.startswith("SAN_d")]}

delays = sorted({float(k.split("_p")[0].split("SAN_d")[1]) for k in sel["delays"]})
print("delays present:", delays)

# ---- fig 1: raw reference + overlap spectra --------------------------------
fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
for i, k in enumerate(sel["ref"]):
    fr, s, b = read(k)
    ax[0].plot(fr / 1e6, s * 1e6, lw=0.7, label=f"{k} SIG")
    ax[0].plot(fr / 1e6, b * 1e6, lw=0.7, ls="--", label=f"{k} BL")
    ax[0].set_title("tone2-only reference (line window should be flat)")
fr, s, b = read(sel["ovl"][0])
ax[1].plot(fr / 1e6, s * 1e6, lw=0.7, label=f"{sel['ovl'][0]} SIG")
ax[1].plot(fr / 1e6, b * 1e6, lw=0.7, ls="--", label="BL")
ax[1].set_title("overlap t_gap=0, phi=0  (expect line at 2074.5 MHz)")
for a in ax:
    a.set_ylabel("power (uW)")
    a.legend(fontsize=7)
ax[1].set_xlabel("frequency (MHz)")
fig.tight_layout()
fig.savefig(OUT + "sanity_trace_ref_ovl.png", dpi=120)
print("saved sanity_trace_ref_ovl.png")

# ---- fig 2: difference spectra vs delay -----------------------------------
fig, ax = plt.subplots(1, 1, figsize=(11, 5))
cmap = plt.cm.viridis(np.linspace(0, 1, max(len(delays), 2)))
for c, d in zip(cmap, delays):
    ks = [k for k in sel["delays"] if f"SAN_d{d:+.2f}_" in k]
    if not ks:
        continue
    fr, s, b = read(ks[0])
    diff = s - b
    ax.plot(fr / 1e6, diff * 1e6, lw=0.7, color=c,
            label=f"t_gap={d:+.1f} ns" if len(delays) <= 12 else None)
ax.axvspan(2073, 2076, alpha=0.15, color="r", label="line window")
ax.set_xlabel("frequency (MHz)")
ax.set_ylabel("SIG - BL (uW)")
title = "difference spectra, first phase per delay"
ax.set_title(title)
if len(delays) <= 12:
    ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(OUT + "sanity_diff_vs_delay.png", dpi=120)
print("saved sanity_diff_vs_delay.png")

# ---- fig 3: line power and c1 vs delay ------------------------------------
res = {}
for d in delays:
    rows_p, rows_ph = [], []
    for k in [k for k in keys if f"SAN_d{d:+.2f}_" in k]:
        fr, s, b = read(k)
        diff = s - b
        ph = float(k.split("_p")[1])
        win = np.abs(fr / 1e6 - 2074.5) <= 1.5
        rows_ph.append(ph)
        rows_p.append(diff[win].sum())
    if len(rows_p) >= 8:
        ph = np.array(rows_ph)
        pw = np.array(rows_p)
        c1 = np.mean(pw * np.exp(-1j * ph))
        res[d] = (float(2 * np.abs(c1)), float(np.mean(pw)), float(np.angle(c1)))

fig, (a1x, a2x) = plt.subplots(1, 2, figsize=(11, 4))
if res:
    ds = sorted(res)
    A1 = [res[d][0] * 1e6 for d in ds]
    A0 = [res[d][1] * 1e6 for d in ds]
    PH = [np.degrees(res[d][2]) for d in ds]
    a1x.plot(ds, A1, "o-", label="A1 = |c1|*2")
    a1x.plot(ds, A0, "s-", label="A0 = mean")
    a1x.set_xlabel("t_gap (ns)")
    a1x.set_ylabel("power (uW)")
    a1x.legend()
    a1x.set_title("phase-modulation amplitude vs delay")
    a2x.plot(ds, PH, "o-")
    a2x.set_ylabel("phi0 (deg)")
    a2x.set_xlabel("t_gap (ns)")
    a2x.set_title("phi0 unwound (expect ~2*pi*Delta*t_gap)")
    print("\ndelay sweep:")
    for d in ds:
        print(f"  t_gap={d:+7.2f}  A0={res[d][1]*1e6:8.1f} uW  "
              f"A1={res[d][0]*1e6:8.1f} uW  A1/A0={res[d][0]/max(res[d][1],1e-30):.2f}  "
              f"phi0={res[d][2]:+.2f} rad")
else:
    a1x.text(0.5, 0.5, "no block with >=8 phases", ha="center")
fig.tight_layout()
fig.savefig(OUT + "sanity_A1_phi0_vs_delay.png", dpi=120)
print("saved sanity_A1_phi0_vs_delay.png")
f.close()
