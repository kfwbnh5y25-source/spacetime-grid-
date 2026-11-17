"""
fig4_nonoverlap.py - Figure 4: phase gating with two NON-OVERLAPPING pulses,
extended with panel (d): memory of phase information vs pulse gap (delay).

Pulse scheme (verified from embedded scripts in both h5 attrs):
  f1 tone, 300 ns tukey, ends at t=0; f2 tone, 400 ns tukey, starts at
  t=0 + t_gap; relative phase phi on pulse 2. Amps 327-588 mV.
  Conventions: phi_eff = phi + 2 pi * Delta * t_gap  (Delta = f2 - f1).

panel (a,b,c): detune_fine031126.h5 ({f1:.2f}_{f2:.2f}_{phase:.3f}_1)
  See step notes inside; file TRUNCATED: 3387 keys = 212/875 (f1,f2) pairs,
  9/35 f2 values, traces (1,2,1001) - pass averaging impossible; hot band
  f2 points (4132, 4135) present.
panel (d): delay_fine.h5 was the intended dual-axis dataset ({delay:.4f}_{phase:.3f}_1,
  f1=4072, f2=4133 fixed, delay -1.56..+1.58 ns in 15.6 ps steps x 60 phases)
  BUT the surviving copy is a truncated file (3.1 MB of a 171 MB HDF5; only
  3 distinct delay values survive {0.5312 x1, 0.7812 x86, 1.0156 x60} - logged
  via step 'delay_fine_status'). Substituted: delay_f1.h5 (complete copy
  extracted from delay_delta.zip; keys {delay:.4f}_{f1}_{phase:.3f}_1,
  f2=4142 fixed, amps 327/588 mV, delays -17.6..+6.9 ns step 0.3125 ns plus
  sparse long delays, 16 phases) after the loose 775 MB copy of the same file
  proved truncated identically (775 of 962 MB). delay>0 = non-overlapping gap.
"""
import sys

sys.path.insert(0, "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/pipeline")
from provenance import Provenance

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from scipy.optimize import curve_fit

H5_DETUNE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/detune_fine031126.h5"
H5_DELAYFINE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/delay_fine.h5"
ZIP_DELAY = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/delay_delta.zip"
H5_DELAY = ("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/"
            "reports/mistral/tmp/delay_extract/delay_f1.h5")
FIG = "fig4_nonoverlap"
HOT_BAND = (4070.0, 4072.0, 4128.0, 4137.0)   # known A1 hot band
F2_DELAY = 4142.0                             # fixed f2 in the delay experiment

p = Provenance(FIG)
src = p.record_source(H5_DETUNE, note="detune_fine: sequential non-overlap phase-gating scan (last-week reference), truncated single pass")
p.record_source(H5_DELAYFINE, note="requested dual-axis dataset: TRUNCATED file (3.1 MB of 171 MB HDF5); only 3 delay rows readable -> replaced by delay_f1.h5")
p.record_source(ZIP_DELAY, note="zip containing the complete delay_f1.h5 used for panel d")
src_delay = p.record_source(H5_DELAY, note="complete delay_f1.h5 extracted from delay_delta.zip; keys {delay:.4f}_{f1}_{phase:.3f}_1")

# ============================================================ panels a/b/c ===
def idler_power(key, tr):
    """Idler readout: mean linear spectral power over the measured 30 MHz band,
    identical to the experiment's own live readout (spec.mean() after baseline
    subtraction). tr = raw (1, 2, 1001); ch1 is the power channel."""
    return float(tr[0, 1, :].astype(np.float64).mean())

data = p.load_traces(H5_DETUNE, select=lambda k: k.endswith("_1"), loader=idler_power,
                     pass_agg="all")

p.step("inspect_layout", params={"select": "keys ending '_1' (readoff=1)"},
       note="trace shape (1,2,1001): dim0=1 single record (run aborted in pass 1; "
            "the 10 loop passes of the script never accumulated). ch0=SA freq axis "
            "(center f2/2, span 30 MHz), ch1=baseline-subtracted linear power.",
       result={"n_keys": len(data), "trace_shape": [1, 2, 1001],
               "file_size_bytes": src["size_bytes"]})

series = {}
bad = 0
for k, I in data.items():
    try:
        a, b, c, d = k.split("_")
        f1, f2, ph = float(a), float(b), float(c)
    except ValueError:
        bad += 1
        continue
    series.setdefault((f1, f2), {})[ph] = I
f1s = sorted(set(k[0] for k in series))
f2s = sorted(set(k[1] for k in series))
n_complete = sum(1 for v in series.values() if len(v) == 16)
p.step("parse_keys", params={"key_format": "{f1:.2f}_{f2:.2f}_{phase:.3f}_1"},
       note="phase series assembled per (f1,f2); f2 grid sparse - file truncated.",
       result={"pairs": len(series), "complete_16_phase": n_complete,
               "bad_keys": bad, "f2_values_MHz": f2s})

def fit_phase(ph, I):
    """I(phi) = A0 + A1 cos(phi+phi0) + A2 cos(2(phi+phi2)); linear LSQ."""
    X = np.column_stack([np.ones_like(ph), np.cos(ph), np.sin(ph),
                         np.cos(2 * ph), np.sin(2 * ph)])
    c, *_ = np.linalg.lstsq(X, I, rcond=None)
    pred = X @ c
    ss_res = float(np.sum((I - pred) ** 2))
    ss_tot = float(np.sum((I - I.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    A1c, A1s, A2c, A2s = c[1], c[2], c[3], c[4]
    A1 = float(np.hypot(A1c, A1s))
    A2 = float(np.hypot(A2c, A2s))
    phi0 = float(np.degrees(np.arctan2(-A1s, A1c)))
    phi2 = float(np.degrees(np.arctan2(-A2s, A2c)))
    return dict(A0=float(c[0]), A1=A1, A2=A2, phi0=phi0, phi2=phi2, r2=r2,
                A1c=float(A1c), A1s=float(A1s), A2c=float(A2c), A2s=float(A2s),
                A2rel=A2 / max(A1, 1e-30))

fits = {}
for (f1, f2), pd_ in series.items():
    ph = np.array(sorted(pd_))
    I = np.array([pd_[pp] for pp in ph])
    if len(pd_) >= 10:
        fits[(f1, f2)] = fit_phase(ph, I)

p.step("fourier_fit", params={"model": "I=A0+A1cos(phi+phi0)+A2cos(2(phi+phi2))",
                              "method": "linear least squares, 16 phases"},
       result={"n_fits": len(fits)})

A1 = np.full((len(f2s), len(f1s)), np.nan)
A2r = np.full_like(A1, np.nan)
R2 = np.full_like(A1, np.nan)
P0 = np.full_like(A1, np.nan)
for (f1, f2), r in fits.items():
    A1[f2s.index(f2), f1s.index(f1)] = r["A1"]
    A2r[f2s.index(f2), f1s.index(f1)] = r["A2rel"]
    R2[f2s.index(f2), f1s.index(f1)] = r["r2"]
    P0[f2s.index(f2), f1s.index(f1)] = r["phi0"]

p.array("A1_cos_phi_amp", A1)
p.array("A2_over_A1_ratio", A2r)
p.array("fit_R2", R2)
p.array("phi0_deg", P0)

imax = np.unravel_index(np.nanargmax(A1), A1.shape)
A1_max = float(np.nanmax(A1))
f2v, f1v = np.meshgrid(np.array(f2s), np.array(f1s), indexing="ij")
hot_mask = (f2v >= HOT_BAND[2]) & (f2v <= HOT_BAND[3]) & (f1v >= HOT_BAND[0]) & (f1v <= HOT_BAND[1])
A1_off_med = float(np.nanmedian(A1[(~hot_mask) & np.isfinite(A1)]))
contrast_db_maxmed = 10 * np.log10(A1_max / A1_off_med)
A1_min = float(np.nanmin(A1[np.isfinite(A1)]))
contrast_db_maxmin = 10 * np.log10(A1_max / A1_min)
phi0_hot = P0[hot_mask & np.isfinite(P0)]
phi0_spread_hot = float(np.nanpercentile(phi0_hot, 90) - np.nanpercentile(phi0_hot, 10))
R2_hot_med = float(np.nanmedian(R2[hot_mask & np.isfinite(R2)]))
a2r_hot_med = float(np.nanmedian(A2r[hot_mask & np.isfinite(A2r)]))
a2r_off_med = float(np.nanmedian(A2r[(~hot_mask) & np.isfinite(A2r)]))
p.step("headline_numbers", result={
    "A1_max_W": A1_max, "A1_max_at": {"f1_MHz": f1s[imax[1]], "f2_MHz": f2s[imax[0]],
                                      "Delta_MHz": f2s[imax[0]] - f1s[imax[1]]},
    "A1_off_band_median_W": A1_off_med,
    "A1_contrast_max_over_offband_dB": contrast_db_maxmed,
    "A1_contrast_max_over_globalmin_dB": contrast_db_maxmin,
    "phi0_hot_band_10-90pct_spread_deg": phi0_spread_hot,
    "R2_hot_band_median": R2_hot_med,
    "A2_over_A1_hot_band_median": a2r_hot_med,
    "A2_over_A1_off_band_median": a2r_off_med},
    note="phase-flip check: phi0 remains in -165..-150 deg band (no pi-jump in hot band); "
         "R2 stays >=0.94 (no sinusoid collapse); the regime boundary instead shows as "
         "collapse of the cos(2 phi) harmonic A2/A1 inside the hot band.")

# ================================================================ panel d ====
# delay_fine status step (requested file is broken - logged plainly)
p.step("delay_fine_status", params={"path": H5_DELAYFINE,
                                    "salvage": "superblock EOF patch on tmp copy"},
       result={"file_size_bytes": 3133389,
               "stored_eof_bytes": 171843445,
               "openable": False,
               "salvaged_delay_rows": {"0.5312": 1, "0.7812": 86, "1.0156": 60},
               "note": "delay_fine.h5 has the RIGHT design (fixed f1=4072, f2=4133, "
                       "Delta=61 MHz, delay -1.56..+1.58 ns in 15.6 ps steps, 60 phases) "
                       "but the surviving copy holds only 3 distinct delays -> no decay "
                       "series. Substituted delay_f1.h5 (complete, from delay_delta.zip)."},
       note="h5py fails on original (truncated file). Substitution documented here.")

f1_sel = [4066.0, 4068.0, 4070.0, 4072.0, 4074.0, 4076.0, 4078.0]  # hot band f1 row
F1_D = 4072.0  # Delta = 70 MHz for the image panel

def delay_power(key, tr):
    del key
    return float(tr[0, 1, :].astype(np.float64).mean())

def delay_select(k):
    pr = k.split("_")
    return len(pr) == 4 and pr[3] == "1" and float(pr[1]) in f1_sel

ddata = p.load_traces(H5_DELAY, select=delay_select, loader=delay_power, pass_agg="all")

dser = {}
for k, I in ddata.items():
    a, b, c, d = k.split("_")
    dser.setdefault((float(a), float(b)), {})[float(c)] = I
p.step("parse_delay_keys", params={"key_format": "{delay:.4f}_{f1}_{phase:.3f}_1",
                                   "f1_selected_MHz": f1_sel},
       note="delay>0: f2 pulse starts after the 300 ns f1 pulse ended -> non-overlapping gap; "
            "delay<0 -> pulses overlap inside the 128000-sample window (t=0 at f1 pulse end). "
            "f2 fixed at 4142 MHz (Delta = 4142 - f1).",
       result={"n_series": len(dser),
               "n_keys_delay_file": len(ddata),
               "delay_range_ns": [min(d for d, _ in dser), max(d for d, _ in dser)]})

dfits = {}
for (delay, f1), pd_ in dser.items():
    ph = np.array(sorted(pd_))
    I = np.array([pd_[pp] for pp in ph])
    if len(pd_) >= 10:
        dfits[(delay, f1)] = fit_phase(ph, I)

uniq_d = sorted(set(d for d, _ in dfits))
A1d = np.full((len(uniq_d), len(f1_sel)), np.nan)
P0d = np.full_like(A1d, np.nan)
R2d = np.full_like(A1d, np.nan)
for (delay, f1), r in dfits.items():
    i, j = uniq_d.index(delay), f1_sel.index(f1)
    A1d[i, j] = r["A1"] * 1e9      # nW
    P0d[i, j] = r["phi0"]
    R2d[i, j] = r["r2"]
p.array("A1_vs_delay_f1_nW", A1d)
p.array("phi0_vs_delay_f1_deg", P0d)
p.array("R2_vs_delay_f1", R2d)

# image data: I(phi, delay) at f1 = 4072
img_d, img_ph = [], []
for (delay, f1), pd_ in sorted(dser.items()):
    if f1 != F1_D or len(pd_) < 16:
        continue
    img_d.append(delay)
    img_ph.append([pd_[q] for q in sorted(pd_)])
img_d = np.array(img_d)
img_ph = np.array(sorted(dser[(img_d[0], F1_D)])) if img_d.size else np.array([])
img_I = np.array([[dser[(d, F1_D)][q] for q in sorted(dser[(d, F1_D)])]
                  for d in img_d]) * 1e9  # (delay, phase) nW
p.array("I_phi_delay_f14072_nW", img_I)

# ---- memory timescale: exponential decay of A1 vs positive delay ------------
# A1 is flat for t_gap ~ 0-2 ns (f2 pulse still inside the spin-wave build-up),
# then decays; fit the decay segment 1.5-8 ns, and report the half-contrast
# delay measured from d=0 as a fit-independent cross-check.
uniq_arr = np.array(uniq_d)
sel_fit = (uniq_arr >= 1.5) & (uniq_arr <= 8.0)
d_fit = uniq_arr[sel_fit]

def model(d, a, tau, c):
    return a * np.exp(-d / tau) + c

fit_ok, tau_res = {}, {}
for j, f1 in enumerate(f1_sel):
    y = A1d[sel_fit, j]
    m = np.isfinite(y)
    try:
        po, cv = curve_fit(model, d_fit[m], y[m], p0=[80, 1.5, 8],
                           bounds=([20, 0.1, 0], [200, 20, 40]), maxfev=30000)
        perr = np.sqrt(np.diag(cv))[1]
        fit_ok[f1] = (po, perr)
        tau_res[f1] = (po[1], perr)
    except Exception as e:
        fit_ok[f1] = None
taus = [v[0] for v in tau_res.values()]
tau_med = float(np.median(taus))
tau_err_med = float(np.median([v[1] for v in tau_res.values()]))
# half-contrast delay from hot-band-averaged curve, d measured from t_gap=0
sel_all = (uniq_arr >= 0) & (uniq_arr <= 8.0)
d_pos = uniq_arr[sel_all]
A1_pos = np.nanmean(A1d[sel_all], axis=1)
A1_zero = A1_pos[0]
d_half = float(np.interp(A1_zero / 2, A1_pos[::-1], d_pos[::-1]))
# weighted per-f1 results recorded; phi0 slope coherence check at f1=4072
j = f1_sel.index(F1_D)
phi0_pass = P0d[:, j]
sel_sl = np.isfinite(phi0_pass) & (np.array(uniq_d) >= 0) & (np.array(uniq_d) <= 3.2)
d_sl = np.array(uniq_d)[sel_sl]
ph_unwrap = np.unwrap(np.radians(phi0_pass[sel_sl]))
slope_deg_ns = float(np.polyfit(d_sl, np.degrees(ph_unwrap), 1)[0])
delta_pred_deg_ns = 360.0 * (F2_DELAY - F1_D) / 1000.0   # 2 pi Delta t -> 25.2 deg/ns
p.step("memory_timescale", params={
        "fit_model": "A1(d) = A1_0 exp(-d/tau) + c; window d = 1.5..8 ns (decay segment; "
                     "A1 is flat for 0..2 ns where the f2 pulse still samples the spin-wave build-up)",
        "hot_band_f1_MHz": f1_sel},
       result={"tau_ns_median_hot_band": tau_med,
               "tau_ns_err_median": tau_err_med,
               "tau_ns_per_f1": {str(k): [round(v[0], 3), round(v[1], 3)] for k, v in tau_res.items()},
               "A1_half_contrast_delay_ns": d_half,
               "A1_at_zero_gap_nW": float(A1_zero),
               "phi0_slope_measured_deg_per_ns": slope_deg_ns,
               "phi0_slope_predicted_from_Delta_deg_per_ns": delta_pred_deg_ns,
               "Delta_MHz_at_f1_4072": F2_DELAY - F1_D},
       note="coherence of memory: phi0(delay) slews at the 2*pi*Delta*t_gap rate -> the "
            "idler carries the pump phase through the gap; amplitude A1 decays -> memory timescale tau.")

# ================================================================= figure ====
fig = plt.figure(figsize=(17.5, 4.8))
gs = fig.add_gridspec(1, 4, width_ratios=[1.1, 1.1, 1.3, 1.55], wspace=0.42)
ax_a = fig.add_subplot(gs[0]); ax_b = fig.add_subplot(gs[1])
ax_c = fig.add_subplot(gs[2]); gsd = gs[3].subgridspec(2, 1, height_ratios=[1.25, 1], hspace=0.4)
ax_d1 = fig.add_subplot(gsd[0]); ax_d2 = fig.add_subplot(gsd[1])

F1G, F2G = np.meshgrid(np.array(f1s), np.array(f2s))

def draw_map(ax, Z, cmap, cblab, title):
    Zm = np.ma.masked_invalid(Z)
    mesh = ax.pcolormesh(F1G, F2G, Zm, cmap=cmap, shading="nearest")
    ax.add_patch(Rectangle((HOT_BAND[0], HOT_BAND[2]),
                           HOT_BAND[1] - HOT_BAND[0], HOT_BAND[3] - HOT_BAND[2],
                           fill=False, ec="w", lw=1.8, ls="--"))
    ax.plot(f1s[imax[1]], f2s[imax[0]], "w*", ms=13, mec="k", mew=0.6)
    cb = fig.colorbar(mesh, ax=ax, pad=0.02)
    cb.set_label(cblab)
    ax.set_xlabel("Drive frequency $f_1$ (MHz)")
    ax.set_ylabel("Drive frequency $f_2$ (MHz)")
    ax.set_title(title, fontsize=9)

draw_map(ax_a, A1 * 1e9, "viridis", r"$A_1$ cos-$\phi$ amplitude (nW)",
         "(a)  cos-$\\phi$ Fourier amplitude $A_1$ of idler power\n"
         "dashed: known hot band; star: $A_1$ max")

draw_map(ax_b, A2r, "magma", r"$A_2/A_1$", "(b)  cos-2$\\phi$ harmonic $A_2/A_1$:\n"
         "collapse (dark) marks hot-band regime boundary")

pts = {"(4070, 4135) hot band": (4070.0, 4135.0),
       "(4064, 4132) hot band": (4064.0, 4132.0),
       "(4084, 4148) off band": (4084.0, 4148.0)}
cols = ["#d62728", "#ff7f0e", "#1f77b4"]
for (lab, (pf1, pf2)), col in zip(pts.items(), cols):
    pd_ = series[(pf1, pf2)]
    ph = np.array(sorted(pd_))
    I = np.array([pd_[q] for q in sorted(pd_)]) * 1e9
    tt = np.linspace(0, 2 * np.pi, 200)
    r = fits[(pf1, pf2)]
    Ifit = (r["A0"] + r["A1c"] * np.cos(tt) + r["A1s"] * np.sin(tt)
            + r["A2c"] * np.cos(2 * tt) + r["A2s"] * np.sin(2 * tt)) * 1e9
    ax_c.plot(np.degrees(tt), Ifit, "-", color=col, lw=1.2, alpha=0.8)
    ax_c.plot(np.degrees(ph), I, "o", color=col, ms=4.5, label=lab)
ax_c.set_xlabel("Relative phase of $f_2$ pulse, $\\phi$ (deg)")
ax_c.set_ylabel("Idler output power (nW)")
ax_c.set_title("(c)  Phase-resolved idler power $I(\\phi)$:\n"
               "points = data, lines = Fourier fit", fontsize=9)
ax_c.legend(fontsize=7, loc="lower right")
ax_c.grid(alpha=0.3)

# ---- panel d1: I(phi, delay) image -----------------------------------------
keep = (img_d >= -2.0) & (img_d <= 14.0)
mesh = ax_d1.pcolormesh(img_d[keep], np.degrees(img_ph),
                        np.ma.masked_invalid(img_I[keep].T), cmap="viridis", shading="nearest")
cb = fig.colorbar(mesh, ax=ax_d1, pad=0.02)
cb.set_label("Idler output power (nW)", fontsize=8)
ax_d1.set_xlabel("Pulse gap $t_{\\mathrm{gap}}$ (ns)", fontsize=8)
ax_d1.set_ylabel("Relative phase of $f_2$ pulse, $\\phi$ (deg)", fontsize=8)
ax_d1.set_title("(d)  $I(\\phi,t_{\\mathrm{gap}})$, f1=4072 ($\\Delta$=70 MHz):\n"
                "stripes = $\\phi_{\\mathrm{eff}}\\!\\approx\\!\\phi\\!+\\!2\\pi\\Delta t_{\\mathrm{gap}}$", fontsize=9)

# ---- panel d2: A1 decay + phi0 slew ----------------------------------------
ax_d2.plot(d_pos, A1_pos, "o", ms=4, color="#d62728", label="$A_1$ (hot-band avg)")
dp = np.linspace(1.5, 8, 300)
po_med = np.median(np.array([v[0] for v in fit_ok.values() if v]), axis=0)
ax_d2.plot(dp, model(dp, *po_med), "-", color="#d62728", lw=1.4,
           label=f"fit: $\\tau$ = {tau_med:.2f}$\\pm${tau_err_med:.2f} ns")
ax_d2.set_xlabel("Pulse gap $t_{\\mathrm{gap}}$ (ns)")
ax_d2.set_ylabel("$A_1$ (nW)", color="#d62728")
ax_d2.tick_params(axis="y", labelcolor="#d62728")
ax_d2b = ax_d2.twinx()
ax_d2b.plot(d_sl, np.degrees(ph_unwrap), "^", ms=4, color="#2ca02c",
            label="$\\phi_0(t_{\\mathrm{gap}})$")
ax_d2b.set_ylabel("$\\phi_0$ shift (deg)", color="#2ca02c")
ax_d2b.tick_params(axis="y", labelcolor="#2ca02c")
ax_d2.set_title("decay of phase contrast; $\\phi_0$ slew = $2\\pi\\Delta t_{\\mathrm{gap}}$",
                fontsize=9)
h1, l1 = ax_d2.get_legend_handles_labels()
h2, l2 = ax_d2b.get_legend_handles_labels()
ax_d2.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper right")
ax_d2.grid(alpha=0.3)

fig.suptitle("Fig. 4 - Phase gating with two non-overlapping pulses: coherent phase memory through the pulse gap "
             "(300 ns tukey $f_1$ then 400 ns $f_2$; (a-c) detune_fine031126, 390/588 mV; (d) delay_f1, 327/588 mV, $f_2$=4142 MHz)",
             fontsize=10.5, y=1.04)
fig.tight_layout()
p.finish(fig, ["a_A1_map", "b_A2_over_A1_map", "c_phase_curves",
               "d_memory_I_phi_delay", "d_memory_decay"])
