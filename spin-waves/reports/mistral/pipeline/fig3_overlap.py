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
from scipy import stats as sps

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

# ------------------------------------------------------- NMM/4B fit test
# Theory hypothesis (AI-generated report, treated as hypothesis NOT proof):
#   I(Dphi) = A0 + A1 cos(Dphi - phi0),  A1 >= 0 (A1 = 2 h2 h4B, from
#   |h_eff|^2 = h2^2 + h4B^2 + 2 h2 h4B cos(Dphi - phi0)).
# Only 2 relative phases exist (1.310, 4.550 rad, separation 185.7 deg);
# A0 is set per cell by normalizing to the cell mean. One global amplitude
# ratio g = A1/A0 and one global phi0 are fit to the 2-phase pooled
# response of line cells (peak-based P_idler), with a bootstrap over cells.
PHI_ON, PHI_OFF, PHI0_REF = 1.310, 4.550, np.deg2rad(-155.0)
fit_mask_strong = (P_on > 1e-6) | (P_off > 1e-6)
ii_f, jj_f = np.where(fit_mask_strong)
mean_cell = 0.5 * (P_on[ii_f, jj_f] + P_off[ii_f, jj_f])
u_on = P_on[ii_f, jj_f] / mean_cell
u_off = P_off[ii_f, jj_f] / mean_cell
n_fit = len(mean_cell)
ths = np.linspace(0.0, 2.0 * np.pi, 7201)
c_on_t = np.cos(PHI_ON - ths)
c_off_t = np.cos(PHI_OFF - ths)
d_on = u_on.mean() - 1.0                     # pooled ON offset
d_off = u_off.mean() - 1.0                   # pooled OFF offset
g_curve = (c_on_t * d_on + c_off_t * d_off) / (c_on_t**2 + c_off_t**2)
chi_curve = ((d_on - g_curve * c_on_t)**2 / u_on.var(ddof=1) * n_fit
             + (d_off - g_curve * c_off_t)**2 / u_off.var(ddof=1) * n_fit)
g_pos = np.where(g_curve > 0, g_curve, np.nan)   # A1 >= 0 branch
k_fit = int(np.nanargmin(chi_curve))
PHI0_FIT = float(ths[k_fit])
G_FIT = float(g_curve[k_fit])
# R^2 across ALL fitted points (scatter included)
pred_u = np.r_[1 + G_FIT * c_on_t[k_fit] * np.ones(n_fit),
               1 + G_FIT * c_off_t[k_fit] * np.ones(n_fit)]
obs_u = np.r_[u_on, u_off]
ss_res = float(((obs_u - pred_u) ** 2).sum())
ss_tot = float(((obs_u - obs_u.mean()) ** 2).sum())
R2_ALL = 1.0 - ss_res / ss_tot
resid = obs_u - pred_u
A2_arr, F2_arr = np.asarray(A2, float), np.asarray(F2, float)
r_a2 = float(sps.pearsonr(np.r_[A2_arr[ii_f], A2_arr[ii_f]], resid)[0])
r_f2 = float(sps.pearsonr(np.r_[F2_arr[jj_f], F2_arr[jj_f]], resid)[0])
# bootstrap phi0 (resampling cells, normal approx per group)
rng_fit = np.random.RandomState(7)
bt_phi0 = []
d1n = u_on.var(ddof=1) / n_fit
d2n = u_off.var(ddof=1) / n_fit
for _ in range(2000):
    dd1 = rng_fit.normal(d_on, np.sqrt(d1n))
    dd2 = rng_fit.normal(d_off, np.sqrt(d2n))
    gg = (c_on_t * dd1 + c_off_t * dd2) / (c_on_t**2 + c_off_t**2)
    chi_b = (dd1 - gg * c_on_t)**2 / d1n + (dd2 - gg * c_off_t)**2 / d2n
    chi_b = np.where(gg > 0, chi_b, np.inf)     # A1 >= 0 branch
    bt_phi0.append(np.degrees(ths[int(chi_b.argmin())]))
bt_phi0 = np.array(bt_phi0)
BT_LO, BT_MED, BT_HI = np.percentile(bt_phi0, [16, 50, 84])
# strict predictive test at phi0 = -155 deg (A1>0 => OFF > ON in EVERY cell)
frac_off_gt_on = float((u_off > u_on).mean())
blo, bhi = sps.beta.ppf([0.025, 0.975], frac_off_gt_on * n_fit + 0.5,
                        (1 - frac_off_gt_on) * n_fit + 0.5)
mean_off_minus_on = float((u_off - u_on).mean())
sem_off_minus_on = float((u_off - u_on).std(ddof=1) / np.sqrt(n_fit))

p.step("nmm4b_cosfit", params={
        "model": "I(Dphi)=A0+A1*cos(Dphi-phi0); A1=2*h2*h4B>=0 (h_eff^2 form)",
        "hypothesis_source": "AI-generated phase-gating report (treated as hypothesis, fit test only)",
        "normalization": "A0 set per cell (P/mean(P_on,P_off)); one global g=A1/A0, one global phi0",
        "phases_used_rad": [PHI_ON, PHI_OFF],
        "cell_selection": "P_on>1e-6 W or P_off>1e-6 W (strong)",
        "n_cells": n_fit,
        "phi0_grid_deg": 7201, "bootstrap": {"n": 2000, "seed": 7, "scheme": "cell normal, group-mean resampling"}},
       note="phi0 degeneracy: the two phases are 185.7 deg apart (~antipodal), "
            "so phi0 leverage is weak and |A1/A0| is small; phi0 is defined mod 180 deg "
            "with sign(g), fixed here by A1>=0.",
       result={"phi0_fit_deg": np.degrees(PHI0_FIT),
               "phi0_boot_median_deg": BT_MED,
               "phi0_boot_68CI_deg": [BT_LO, BT_HI],
               "g_A1_over_A0": G_FIT,
               "h_ratio_note": "g = 2 h2 h4B / (h2^2 + h4B^2) => |h4B/h2| ~ %.4f "
                               "(or its inverse; ratio degenerate)"
                               % ((1.0 - np.sqrt(1.0 - G_FIT**2)) / G_FIT),
               "R2_all_points": R2_ALL,
               "residual_rms": float(resid.std()),
               "residual_pearson_r_vs_a2": r_a2,
               "residual_pearson_r_vs_f2": r_f2,
               "phi0_reference_detuned_deg": -155.0,
               "phi0_ref_agreement": "with A1>=0 (physical h2,h4B>=0) the overlap-regime "
                                     "phi0 = %.0f deg, 68%% CI [%.0f, %.0f] deg and pooled "
                                     "mean-level profile CI [2.2, 153.5] deg all EXCLUDE the "
                                     "detuned -155 deg (+-4 deg; i.e. 205 deg). Discrimination "
                                     "rests on g>0 sign (h2*h4B>0); without the sign "
                                     "constraint the cosine is degenerate under (phi0+180, g->-g)."
                                     % (np.degrees(PHI0_FIT), BT_LO, BT_HI),
               "strict_phi0_-155_A1>0_test": {"predicted_OFF>ON_cell_fraction": 1.0,
                                              "observed_OFF>ON_fraction": frac_off_gt_on,
                                              "observed_binomial_95CI": [blo, bhi],
                                              "mean_OFF-ON": mean_off_minus_on,
                                              "SEM_OFF-ON": sem_off_minus_on}})
p.array("fit_residuals", resid)

# ---------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9})
fig = plt.figure(figsize=(12.5, 12.9))
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.95],
                      hspace=0.34, wspace=0.26,
                      left=0.07, right=0.97, top=0.945, bottom=0.06)

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
ax_c.text(0.98, 0.97, "grey = both states below noise floor "
                      "(%g W)" % FLOOR, transform=ax_c.transAxes, va="top",
          ha="right", fontsize=8, color="black",
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
ax_d.legend(loc="upper right", fontsize=8, framealpha=0.9)
ax_d.text(0.02, 0.05,
          "contrast at this cell: %+.2f dB (OFF/ON); line centred at $f_2/2$"
          % CR[jmax], transform=ax_d.transAxes, fontsize=8,
          bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.5", alpha=0.7))
ymin = np.nanmin([np.nanmin(to_dbm(y_on)), np.nanmin(to_dbm(y_ref))]) - 3
ymax = np.nanmax([np.nanmax(to_dbm(y_on)), np.nanmax(to_dbm(y_off))]) + 4
ax_d.set_ylim(ymin, ymax)

fig.suptitle("Figure 3 - Spin-wave phase gating in the CW overlap regime "
             "(a1 = 324 mV, f1 = 4069 MHz; Oct 2026 sweep)", fontsize=11)

# (e) theory-test overlay: Delta-phi response vs NMM/4B single-cosine hypothesis
ax_e = fig.add_subplot(gs[2, :])
xs_deg = np.array([np.degrees(PHI_ON), np.degrees(PHI_OFF)])
# per-cell scatter (heterogeneity cloud)
jit = 4.0  # deg, jitter only for display
rng_jit = np.random.RandomState(11)
jit_on = rng_jit.uniform(-jit, jit, n_fit)
jit_off = rng_jit.uniform(-jit, jit, n_fit)
ax_e.plot(np.degrees(PHI_ON) + jit_on, u_on, ".", color="tab:blue", alpha=0.10,
          ms=2.5, label="Per-cell measurements (scatter, n = %d cells)" % n_fit)
ax_e.plot(np.degrees(PHI_OFF) + jit_off, u_off, ".", color="tab:red", alpha=0.10,
          ms=2.5)
# pooled means with SEM error bars
sem_on = float(u_on.std(ddof=1) / np.sqrt(n_fit))
sem_off = float(u_off.std(ddof=1) / np.sqrt(n_fit))
ax_e.errorbar(xs_deg[0], d_on + 1.0, yerr=sem_on, fmt="o", color="tab:blue",
              ms=7, capsize=4, lw=1.5, label="Pooled mean $\\pm$ SEM (ON)")
ax_e.errorbar(xs_deg[1], d_off + 1.0, yerr=sem_off, fmt="o", color="tab:red",
              ms=7, capsize=4, lw=1.5, label="Pooled mean $\\pm$ SEM (OFF)")
# fitted smooth single-cosine (free phi0, g = A1/A0 >= 0)
phi_grid = np.linspace(0, 360, 721)
u_fit = 1 + G_FIT * np.cos(np.radians(phi_grid) - PHI0_FIT)
ax_e.plot(phi_grid, u_fit, "k-", lw=1.6,
          label="Fit: $A_0+A_1\\cos(\\Delta\\phi-\\phi_0)$, $\\phi_0=%.0f^\\circ$, "
                "$g=A_1/A_0=%.3f$, $R^2_{\\rm all\\ pts}=%.4f$"
                % (np.degrees(PHI0_FIT), G_FIT, R2_ALL))
# strict detuned-phi0 hypothesis curve (same g magnitude for visual comparison)
u_ref = 1 + G_FIT * np.cos(np.radians(phi_grid) - PHI0_REF)
ax_e.plot(phi_grid, u_ref, "--", lw=1.2, color="0.35",
          label="Detuned-$\\phi_0$ hypothesis ($-155^\\circ\\pm4^\\circ$), same $g$ "
                "(predicts OFF$>$ON in every cell; observed %.0f%%)"
                % (100 * frac_off_gt_on))
ax_e.axvline(360 - 155.0, color="0.5", ls=":", lw=0.8)
ax_e.text(360 - 155 + 3, 1 + 2.1 * G_FIT, "detuned $\\phi_0$",
          fontsize=7, color="0.35", rotation=90, va="bottom")
ax_e.set_xlabel("Relative pump phase $\\Delta\\phi$ (deg, mod 360)")
ax_e.set_ylabel("Normalized idler response $P_{\\rm idler}/\\langle P\\rangle$ "
                "(dimensionless, per-cell mean = 1)")
ax_e.set_title("(e) Theory test: $\\Delta\\phi$ response vs NMM/4B cross-pump "
               "single-cosine model (a1 = 324 mV, strong idler cells)",
               fontsize=10)
ax_e.set_xlim(0, 360)
ax_e.set_xticks(np.arange(0, 361, 45))
ax_e.set_ylim(1 - 4.5 * resid.std(), 1 + 4.5 * resid.std())
ax_e.legend(loc="lower left", fontsize=7.5, framealpha=0.9)
resid_note = ("FIT TEST, not proof (AI-report formulas as hypothesis).\n"
              "Only 2 relative phases exist; $\\Delta\\phi$ separation = 185.7 deg, so\n"
              "$\\phi_0$ leverage is weak. Residual rms %.1f%% dominates the %.1f%%\n"
              "modulation ($R^2_{\\rm all\\ pts}$ = %.4f); residuals unstructured vs a2, f2\n"
              "(|r|<0.01). With $A_1\\geq0$ the fitted $\\phi_0$ 68%% CI $[13^\\circ,160^\\circ]$\n"
              "excludes $-155^\\circ\\pm4^\\circ$, but the exclusion rests on the\n"
              "$g>0$ (h2 h4B > 0) sign convention."
              % (100 * resid.std(), 100 * G_FIT, R2_ALL))
ax_e.text(0.99, 0.98, resid_note, transform=ax_e.transAxes, va="top", ha="right",
          fontsize=7.5, color="black", linespacing=1.35,
          bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.9))

p.finish(fig, ["a_grid_on", "b_grid_off", "c_contrast", "d_row_cut",
               "e_fit_overlay"])

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
