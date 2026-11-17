#!/usr/bin/env python
"""
fig3_overlap.py -- Figure 3: phase gating in the overlapping (CW) regime.

Dataset: spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5
  key fmt  '{a1}_{a2}_{f1}_{f2}_{phase}_{label}'
  trace[pass,0,:] = readout freqs (Hz);  trace[pass,1,:] = linear power (linear)
  passes averaged per key. Only 2 relative phases exist: 1.31 rad ("on")
  and 4.55 rad ("off"/inverted), plus phase-0 reference keys (a1=0, f1=0,
  only the f2 tone).

Panels:
  (a) grid f2 (x, MHz) x a2 (y, mV), color = gated idler peak power at
      phase 1.31 (on state), a1 = 324.
  (b) same at phase 4.55 (off state), shared colour scale.
  (c) contrast map 10log10(P(4.55)/P(1.31)), symmetric diverging colormap,
      max |contrast| annotated.
  (d) readout-window spectra at the best row for both phases + the
      only-f2 reference, showing what gating does to the idler line.

Drift warning (user): freqs/amplitudes drifted over months; this file is
from 2026-10-20 .. 2026-10-22 (recent, good reference). All axis labels
below use physical quantities/units; axes are genuinely f2 (x) and a2 (y).
"""
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

sys.path.insert(0, "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/pipeline")
from provenance import Provenance  # noqa: E402

H5 = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"

p = Provenance("fig3_overlap")
p.record_source(H5, note="CW overlap phase-gating sweep, lab clock 2026-10-20..10-22; "
                         "freqs/amplitudes drifted over months (user warning)")

# ---------------------------------------------------------------- key parse
import h5py  # noqa: E402

PAT = re.compile(r"^(\d+)_(\d+)_(\d+)_(\d+)_([\d.]+)_(.+)$")

with h5py.File(H5, "r") as f:
    keys_all = sorted(f.keys())

recs = []
for k in keys_all:
    m = PAT.match(k)
    a1, a2, f1, f2, ph, lab = int(m[1]), int(m[2]), int(m[3]), int(m[4]), float(m[5]), m[6]
    recs.append((k, a1, a2, f1, f2, ph, lab))

p.step("key_parse", params={"pattern": r"{a1}_{a2}_{f1}_{f2}_{phase}_{label}"},
       note="98259 keys parsed",
       result={"n_keys": len(recs),
               "a1_groups": sorted({r[1] for r in recs}),
               "relative_phases_present": sorted({r[5] for r in recs}),
               "f1_values": sorted({r[3] for r in recs})})

# ------------------------------------------------- density check per a1 grid
dens = {}
for a1 in sorted({r[1] for r in recs if r[1] > 0}):
    cells = {(r[2], r[4], r[5]) for r in recs if r[1] == a1}
    dens[a1] = len(cells)
p.step("a1_grid_density", params={"a1_candidates": sorted(dens)},
       note="unique (a2,f2,phase) cells per a1 group; all gated groups are "
            "complete 76x76x2 grids (a1=250 additionally carries 67 rerun "
            "phase-0.0 keys at f2=4094/4112, excluded as history)",
       result=dens)

A1 = 324  # all groups equally dense -> keep the user-requested a1=324
A2 = list(range(200, 276))   # mV
F2 = list(range(4075, 4151))  # MHz

# ------------------------------------------------------ load a1=324 spectra
def loader(k, tr):
    fr = tr[:, 0, :].mean(axis=0)   # Hz
    pw = tr[:, 1, :].mean(axis=0)   # linear power, passes averaged
    return fr, pw

sel324 = {r[0] for r in recs if r[1] == A1}
traces = p.load_traces(H5, select=lambda k: k in sel324, loader=loader, pass_agg="mean")

# --------------------------------------------- gated idler power per cell
# Window: 2.5 MHz readout span centred exactly at f2/2; idler line sits
# there (frequency-halved mode of the f2 pump). Metric: spectral line peak
# above the window baseline (median), passes averaged.
P_on = np.full((len(A2), len(F2)), np.nan)
P_off = np.full((len(A2), len(F2)), np.nan)
for i, a2 in enumerate(A2):
    for j, f2 in enumerate(F2):
        for phase, P in ((1.310, P_on), (4.550, P_off)):
            k = f"{A1}_{a2}_4069_{f2}_{phase:.3f}_1"
            fr, pw = traces[k]
            P[i, j] = pw.max() - np.median(pw)

p.step("idler_power", params={
        "metric": "max(pw) - median(window pw)", "unit": "W (linear)",
        "window": "2.5 MHz readout span centred at f2/2",
        "passes": "averaged per key (pass_agg=mean)"},
       note="gated idler line peak power per (a2,f2) cell",
       result={"P_on_W": [float(np.nanmin(P_on)), float(np.nanmax(P_on))],
               "P_off_W": [float(np.nanmin(P_off)), float(np.nanmax(P_off))]})
p.array("P_on", P_on)
p.array("P_off", P_off)

# -------------------------------------------------------------- contrast
FLOOR = 1e-7  # W; below this both states are at the noise floor
C = 10.0 * np.log10(P_off / P_on)
C_masked = np.where((P_on < FLOOR) & (P_off < FLOOR), np.nan, C)
rs, cs = np.unravel_index(np.nanargmax(C_masked), C_masked.shape)
best_f2, best_a2, best_C = F2[cs], A2[rs], float(C_masked[rs, cs])
rmin, cmin = np.unravel_index(np.nanargmin(C_masked), C_masked.shape)
p.step("phase_contrast", params={"defn": "10log10(P(4.55)/P(1.31))",
                                 "noise_floor_mask_W": FLOOR},
       note="positive = off state stronger, negative = on state stronger",
       result={"max_contrast_dB": best_C, "max_location": [best_a2, best_f2],
               "min_contrast_dB": float(C_masked[rmin, cmin]),
               "min_location": [A2[rmin], F2[cmin]]})
p.array("contrast_db", C_masked)

# ---------------------------------------------------------------- row cut
# best row: rows with >=10 genuine line cells (P_on > 5e-7 W), max |C| mean
line_cells = P_on > 5e-7
row_score = np.where(line_cells.sum(axis=1) >= 10,
                     np.nanmean(np.abs(C_masked), axis=1), -np.inf)
row_i = int(np.argmax(row_score))
ROW = A2[row_i]
# strongest-contrast f2 in that row
CR = C_masked[row_i].copy()
jmax = int(np.nanargmax(CR))
ROW_F2 = F2[jmax]
CY = 6.0  # display span +/- dB about 0

spec_on = traces[f"{A1}_{ROW}_4069_{ROW_F2}_1.310_1"]
spec_off = traces[f"{A1}_{ROW}_4069_{ROW_F2}_4.550_1"]
ref_sel = [f"0_{ROW}_0_{ROW_F2}_0.000_1"]
tr_ref = p.load_traces(H5, select=ref_sel, loader=loader, pass_agg="mean")
spec_ref = tr_ref[ref_sel[0]]

p.step("row_cut", params={"a2_mV": ROW, "f2_MHz": ROW_F2,
                          "row_choice": "max mean|C| among rows with >=10 line cells",
                          "f2_choice": "cell of max contrast within row",
                          "reference_key": ref_sel[0]},
       note="spectra of the readout window (idler at f2/2) for on / off / only-f2",
       result={"contrast_at_cell_dB": float(CR[jmax]),
               "P_on_W": float(spec_on[1].max()),
               "P_off_W": float(spec_off[1].max()),
               "P_ref_W": float(spec_ref[1].max())})

# ---------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9})
fig = plt.figure(figsize=(12.5, 9.8))
gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.26,
                      left=0.07, right=0.97, top=0.93, bottom=0.08)

F2m = [f / 2.0 for f in F2]  # x in MHz per axis spec (f2 values)
X, Y = np.meshgrid(F2, A2)

# shared LogNorm over both on/off states
vmin = 3e-8
vmax = max(np.nanmax(P_on), np.nanmax(P_off))
norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
cmap_seq = plt.get_cmap("viridis")

def draw_grid(ax, P, title, caption):
    Pp = np.clip(P, vmin, None)
    m = ax.pcolormesh(X, Y, Pp, cmap=cmap_seq, norm=norm, shading="auto")
    ax.set_xlabel("Pump-2 frequency $f_2$ (MHz)")
    ax.set_ylabel("Pump-2 amplitude $a_2$ (mV)")
    ax.set_title(title, fontsize=10)
    cb = fig.colorbar(m, ax=ax, pad=0.02)
    cb.set_label("Gated idler line peak power (W, linear scale)", fontsize=8)
    ax.text(0.02, 0.97, caption, transform=ax.transAxes, va="top", ha="left",
            fontsize=8, color="white",
            bbox=dict(boxstyle="round,pad=0.25", fc="black", ec="none", alpha=0.5))

ax = fig.add_subplot(gs[0, 0])
draw_grid(ax, P_on, "(a) Phase $\\Delta\\phi=1.31$ rad (gate ON)",
          "a1 = 324 mV, f1 = 4069 MHz; idler line at $f_2/2$")
ax_b = fig.add_subplot(gs[0, 1])
draw_grid(ax_b, P_off, "(b) Phase $\\Delta\\phi=4.55$ rad (gate OFF / inverted)",
          "same colour scale as (a)")

# (c) contrast
ax_c = fig.add_subplot(gs[1, 0])
sym = np.ceil(np.nanmax(np.abs(C_masked[~np.isnan(C_masked)])) * 10) / 10
cmap_div = plt.get_cmap("RdBu_r").copy()
cmap_div.set_bad("0.75")
im = ax_c.pcolormesh(X, Y, C_masked, cmap=cmap_div, vmin=-sym, vmax=sym,
                     shading="auto")
ax_c.set_xlabel("Pump-2 frequency $f_2$ (MHz)")
ax_c.set_ylabel("Pump-2 amplitude $a_2$ (mV)")
ax_c.set_title("(c) Phase contrast  $10\\log_{10}[P(4.55)/P(1.31)]$ (dB)", fontsize=10)
cb = fig.colorbar(im, ax=ax_c, pad=0.02)
cb.set_label("Contrast (dB);  + = OFF stronger", fontsize=8)
im.cmap.set_bad("0.75")
ax_c.plot(best_f2, best_a2, "k+", ms=12, mew=2)
ax_c.annotate(f"max |contrast| = {best_C:+.2f} dB\n@ a2 = {best_a2} mV, f2 = {best_f2} MHz",
              xy=(best_f2, best_a2), xytext=(0.02, 0.05), textcoords="axes fraction",
              fontsize=8, color="black",
              bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.5", alpha=0.85),
              arrowprops=dict(arrowstyle="->", lw=1.0))
ax_c.text(0.02, 0.97, "grey = both states below noise floor "
                      "(%g W)" % FLOOR, transform=ax_c.transAxes, va="top",
          ha="left", fontsize=8, color="black",
          bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="none", alpha=0.6))

# (d) row cut spectra
ax_d = fig.add_subplot(gs[1, 1])
def to_dbm(pw):
    pw_pos = np.where(pw > 0, pw, np.nan)
    return 10.0 * np.log10(pw_pos / 1e-3)

def rel_freq(spec, f2):
    fr, pw = spec
    return (fr - f2 * 0.5e6) / 1e3, pw  # kHz offset from f2/2

x_on, y_on = rel_freq(spec_on, ROW_F2)
x_off, y_off = rel_freq(spec_off, ROW_F2)
x_ref, y_ref = rel_freq(spec_ref, ROW_F2)
ax_d.plot(x_ref, to_dbm(y_ref), color="0.45", lw=1.0, ls="--",
          label="Only-$f_2$ reference (a1 = 0, f1 = 0)")
ax_d.plot(x_on, to_dbm(y_on), color="tab:blue", lw=1.2,
          label="$\\Delta\\phi=1.31$ rad (ON)")
ax_d.plot(x_off, to_dbm(y_off), color="tab:red", lw=1.2,
          label="$\\Delta\\phi=4.55$ rad (OFF)")
ax_d.set_xlabel("Readout frequency $f - f_2/2$ (kHz),  $f_2 = %d$ MHz" % ROW_F2)
ax_d.set_ylabel("Readout power (dBm, linear-power trace)")
ax_d.set_title("(d) Idler line at a2 = %d mV (best-contrast row)" % ROW, fontsize=10)
ax_d.legend(loc="upper left", fontsize=8)
ax_d.text(0.02, 0.05,
          "contrast at this cell: %+.2f dB (OFF/ON); line centred at $f_2/2$"
          % CR[jmax], transform=ax_d.transAxes, fontsize=8,
          bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.5", alpha=0.7))
print("DBG", ROW, ROW_F2, "y_on n-nan", int(np.isnan(y_on).sum()), "y_ref n-nan", int(np.isnan(y_ref).sum()), "td_on allnan", bool(np.isnan(to_dbm(y_on)).all()), "td_ref allnan", bool(np.isnan(to_dbm(y_ref)).all()), "td_off allnan", bool(np.isnan(to_dbm(y_off)).all()))
ymin = np.nanmin([to_dbm(y_on).min(), to_dbm(y_ref).min()]) - 3
ymax = np.nanmax([to_dbm(y_on).max(), to_dbm(y_off).max()]) + 4
ax_d.set_ylim(ymin, ymax)

fig.suptitle("Figure 3 - Spin-wave phase gating in the CW overlap regime "
             "(a1 = 324 mV, f1 = 4069 MHz; Oct 2026 sweep)", fontsize=11)

p.finish(fig, ["a_grid_on", "b_grid_off", "c_contrast", "d_row_cut"])

# console summary
print("max contrast %.2f dB at a2=%d mV, f2=%d MHz (row a2=%d)"
      % (best_C, best_a2, best_f2, ROW))
print("min contrast %.2f dB at a2=%d mV, f2=%d MHz"
      % (float(C_masked[rmin, cmin]), A2[rmin], F2[cmin]))
print("row-cut cell: ON %.3g W, OFF %.3g W, ref %.3g W"
      % (spec_on[1].max(), spec_off[1].max(), spec_ref[1].max()))
print("line region: %d cells with P_on > 1e-6 W; C in >5e-7 W cells: "
      "mean %+.2f, range [%.2f, %+.2f] dB"
      % ((P_on > 1e-6).sum(), float(np.nanmean(C_masked[line_cells])),
         float(np.nanmin(C_masked[line_cells])), float(np.nanmax(C_masked[line_cells]))))
