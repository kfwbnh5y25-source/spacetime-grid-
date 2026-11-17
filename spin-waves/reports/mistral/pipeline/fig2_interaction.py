# -*- coding: utf-8 -*-
"""
fig2_interaction.py -- Figure 2 (MAIN, star panels): idler scaling law + (a1,a2) maps.

Restructured per figure assessment (2025-11-17):
  - STAR PANELS stay in Fig 2: (a) scaling law P_m vs a1 (log-log, measured slopes),
    (b-d) 2D pcolormesh maps over (a1, a2) with colorbar P_idler (= P_m).
  - The (a2,f2) gain/contrast maps, D vs ND logistic fits and the threshold-reduction
    line cut moved to figS_threshold_reduction (numbers carried unchanged there).

Data (Paper A "Nonlinear magnon spectroscopy" idler data):
  spin-waves/threshold_reduction/data.h5
  key format '{a1}_{a2}_{f1}_{f2}_{phase}_{ro}_{tim}',
  ro in {EX, f1, f2, m, p} (mixed-frequency readout channels), tim in {D, ND}.
  trace[pass,0,:] = frequency axis (Hz); trace[pass,1,:] = linear power spectrum.
  P_m   = integrated despiked power in the idler channel m = (2 f1 - f2)/2
  P_f2/2= integrated despiked power in the target channel f2/2
  Law:  P_m ∝ a1^2 P_{f2/2}  (4-magnon pre-excitation of the idler mode),
  measured log-log slopes 1.83-1.94 over >2 decades in a1 (quoted range).

Reference style: pcolormesh rasterized=True, NaN masked-invalid, horizontal top colorbar.
Processing follows spin-waves/threshold_reduction/analyze_threshold_reduction.py:
15-bin Hampel despiking (5 sigma), active-region filter a1>=25, a2>=350, P_f2/2>2e-4.
"""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.ndimage import median_filter
from scipy.stats import linregress

sys.path.insert(0, "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/pipeline")
from provenance import Provenance  # noqa: E402

H5 = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/threshold_reduction/data.h5"

p = Provenance("fig2_interaction")
src = p.record_source(H5, note="Paper A two-tone idler data; keys {a1}_{a2}_{f1}_{f2}_{phase}_{ro}_{tim}")
p.step("restructure_note",
       params={"moved_to": "figS_threshold_reduction",
               "moved_content": "(a2,f2) gain/contrast maps (phase gating sweep), D-vs-ND "
                                "logistic fits, threshold-reduction line cut"},
       note="Fig 2 split per figure assessment; every v1 number carried unchanged "
            "(step 'carried_numbers_from_v1')",
       result={})
p.step("carried_numbers_from_v1",
       params={"note": "extracted in v1 of fig2_interaction (phase gating file "
                       "f1_4069_5_a1_a2_f2_sweep.h5); unchanged by the restructure"},
       note="idler in-band gain max 20.53 dB (total) / 26.53 dB (peak) at f2=4111; "
            "phase contrast max +1.81 dB (a2=265, f2=4105), min -2.38 dB (263, 4133), "
            "peak metric +3.47 dB; threshold crossings: 10 dB gain at f2=4111 -> a2 251.20 "
            "(phase 1.31) vs 249.87 (phase 4.55) = 1.33 dB = 14.2% of drive amplitude "
            "= 1.77% of the a2 sweep; all-f2 scan at 10 dB: max 5.40 dB (~46.2% drive, f2=4105, "
            "total metric), 5.20 dB (f2=4123, peak metric); 3 dB/6 dB thresholds give 12.4-76% "
            "(noisy, threshold-definition dependent)",
       result={"max_gain_dB_total_metric": 20.530803512318307,
               "max_gain_dB_peak_metric": 26.52741622532229,
               "max_contrast_dB": 1.8092374023439437,
               "min_contrast_dB": -2.376139068278339,
               "max_contrast_peak_metric_dB": 3.4682289150997296,
               "threshold_reduction_f2_4111_10dB_dB": 1.3307991563984558,
               "threshold_reduction_f2_4111_10dB_lin_amp_pct": 14.20538253452338,
               "all_f2_scan_thr10dB_total_max_dB": 5.404494652255494,
               "all_f2_scan_thr10dB_total_max_f2_MHz": 4105,
               "all_f2_scan_thr3dB_total_max_dB": 12.416828432930913,
               "all_f2_scan_thr6dB_peak_max_dB": 12.511557174083123})

# ----------------------------------------------------------------- loading
def hampel_despike(spec, half_win=7, n_sig=5.0):
    med = median_filter(spec, size=2 * half_win + 1, mode="nearest")
    diff = np.abs(spec - med)
    mad = max(float(np.median(diff)) / 0.6745, 1e-11)
    mask = diff > n_sig * mad
    if mask.any():
        for idx in np.where(mask)[0]:
            mask[max(0, idx - 2): min(len(spec), idx + 3)] = True
        out = spec.copy()
        out[mask] = med[mask]
        return out
    return spec


def key_metrics(key, trace):
    """Loader: mean over re-run passes, Hampel-despiked integrated channel power."""
    raw = trace[:, 1, :].astype(np.float64).mean(axis=0)
    ro = key.split("_")[5]
    clean = raw if ro == "EX" else hampel_despike(raw)
    return float(clean.sum())

traces = p.load_traces(H5, select=lambda k: True, loader=key_metrics, pass_agg="mean")

rows = {}
for k, s in traces.items():
    a1, a2, f1, f2, ph, ro, tim = k.split("_")
    rows.setdefault((ro, int(f2), tim), {})[(int(a1), int(a2))] = s

p.step("load_despiked",
       params={"hampel": "15-bin (half_win=7), 5 sigma", "pass_agg": "mean",
               "channels": ["EX", "f1", "f2", "m", "p"], "timings": ["D", "ND"]},
       note="P_m = idler channel 'm' integral; P_f2/2 = target channel 'f2' integral",
       result={"n_keys": len(traces),
               "groups": sorted(set(rows.keys()))[:8]})

# ------------------------------------------------- star panel (a): scaling law
bins_a1 = [(25, 60), (60, 100), (100, 140), (140, 180), (180, 220), (220, 260), (260, 300), (300, 355)]
slopes, c_fits, pts_store = {}, {}, {}
for f2 in [4125, 4135]:
    pts = []
    for tim in ["D", "ND"]:
        dm, df2 = rows[("m", f2, tim)], rows[("f2", f2, tim)]
        for k in dm:
            if k[0] >= 25 and k[1] >= 350 and k in df2:
                pm, pf2 = dm[k], df2[k]
                if pf2 > 2e-4 and pm > 0:
                    pts.append((k[0], pm / pf2, pm / (k[0] ** 2 * pf2), k[1]))
    pts = np.array(pts)
    slope, icpt, _, _, se = linregress(np.log10(pts[:, 0]), np.log10(pts[:, 1]))
    c_med = float(np.median(pts[:, 2]))
    slopes[f2] = {"slope": float(slope), "slope_se": float(se), "r": float(icpt),
                  "n_pts": int(len(pts)), "C_m__x1e6_mV-2": c_med * 1e6}
    c_fits[f2] = c_med
    pts_store[f2] = pts

measured = sorted(v["slope"] for v in slopes.values())
res_scale = {"per_f2": slopes,
             "measured_slope_range": [measured[0], measured[-1]],
             "quoted_range": [1.83, 1.94],
             "in_quoted_range": bool(1.83 <= measured[0] and measured[-1] <= 1.94)}
p.step("panel_a_scaling_law",
       params={"y": "P_m / P_f2/2 (idler-to-target ratio)", "x": "a1 (mV), log-log",
               "filter": "a1>=25, a2>=350, P_f2/2>2e-4", "pairs": "D+ND combined per f2"},
       note="slope of log10(P_m/P_f2/2) vs log10(a1); law predicts 2 (P_m ∝ a1^2 P_f2/2)",
       result=res_scale)
p.array("fig2_a_scaling_points_f2_4125", pts_store[4125])
p.array("fig2_a_scaling_points_f2_4135", pts_store[4135])

# ------------------------------------------------- star panels (b)-(g): (a1,a2) maps
a1_edges = np.linspace(0, 355, 22)
a2_edges = np.linspace(0, 455, 25)
vmin, vmax = 5e-7, 2e-4

def bin_2d(pts_x, pts_y, vals, x_edges, y_edges, min_val=1e-7):
    grid = np.full((len(y_edges) - 1, len(x_edges) - 1), np.nan)
    for iy in range(len(y_edges) - 1):
        my = (pts_y >= y_edges[iy]) & (pts_y < y_edges[iy + 1])
        for ix in range(len(x_edges) - 1):
            m = my & (pts_x >= x_edges[ix]) & (pts_x < x_edges[ix + 1])
            if m.any():
                grid[iy, ix] = max(float(np.median(vals[m])), min_val)
    return grid

maps = {}
map_res = {}
for f2 in [4125, 4135]:
    dm, df2 = rows[("m", f2, "D")], rows[("f2", f2, "D")]
    common = [k for k in dm if k in df2 and k[0] > 0 and k[1] > 0]
    a1s = np.array([k[0] for k in common], float)
    a2s = np.array([k[1] for k in common], float)
    pms = np.array([dm[k] for k in common], float)
    pf2s = np.array([df2[k] for k in common], float)
    act = (a1s >= 25) & (a2s >= 350) & (pf2s > 2e-4) & (pms > 0)
    c_fit = float(np.median(pms[act] / (a1s[act] ** 2 * pf2s[act])))
    pred = c_fit * a1s ** 2 * pf2s
    g_pm = bin_2d(a1s, a2s, pms, a1_edges, a2_edges)
    g_pred = bin_2d(a1s, a2s, pred, a1_edges, a2_edges)
    g_ratio = np.full_like(g_pm, np.nan)
    valid = (g_pm > 5e-6) & (g_pred > 5e-6)
    g_ratio[valid] = np.log10(g_pm[valid] / g_pred[valid])
    maps[f2] = (g_pm, g_pred, g_ratio)
    map_res[f2] = {"C_fit__x1e6_mV-2": c_fit * 1e6,
                   "median_log10_ratio": float(np.nanmedian(g_ratio)),
                   "rms_log10_ratio": float(np.nanstd(g_ratio))}
    p.array("fig2_map_measured_f2_%d" % f2, g_pm)
    p.array("fig2_map_model_f2_%d" % f2, g_pred)
    p.array("fig2_map_logratio_f2_%d" % f2, g_ratio)
p.step("panels_bd_pcolormesh_maps",
       params={"binning": "median in 21 x 24 bins over a1=[0,355] mV, a2=[0,455] mV",
               "style": "rasterized pcolormesh, LogNorm masked-invalid, horizontal top colorbar"},
       note="measured P_m vs model C*a1^2*P_f2/2 on identical colorscale (D timing)",
       result=map_res)

# ----------------------------------------------------------------- figure
plt.rcParams.update({
    "font.size": 9.5, "axes.titlesize": 10.5, "axes.labelsize": 9.5,
    "legend.fontsize": 7.8, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "axes.spines.top": False, "axes.spines.right": False,
})
fig = plt.figure(figsize=(13.8, 10.6))
gs = fig.add_gridspec(3, 3, height_ratios=[1.35, 1, 1], hspace=0.75, wspace=0.30)

colors = {(4125, "D"): "#1f77b4", (4125, "ND"): "#5ca4d8",
          (4135, "D"): "#d62728", (4135, "ND"): "#f28e8e"}

ax = fig.add_subplot(gs[0, :])
for f2, mkt in [(4125, "o"), (4135, "^")]:
    pts = pts_store[f2]
    ax.scatter(pts[:, 0], pts[:, 1], s=8, alpha=0.30, marker=mkt,
               color=colors[(f2, "D")], edgecolors="none",
               label="f2=%d D+ND (slope %.2f$\\pm$%.2f)" % (f2, slopes[f2]["slope"], slopes[f2]["slope_se"]))
    pts42 = pts[(pts[:, 1] > 2.5e-4) & (pts[:, 0] >= 220)]
    ax.scatter(pts42[:, 0], pts42[:, 1], s=8, alpha=0.15, marker=mkt,
               color=colors[(f2, "ND")], edgecolors="none")
ax.set_xscale("log"); ax.set_yscale("log")
for f2, col in [(4125, "#1f77b4"), (4135, "#d62728")]:
    C = c_fits[f2]
    ag = np.linspace(25, 355, 100)
    ax.plot(ag, C * ag ** 2, color=col, lw=2.0,
            label=r"$\propto a_1^2$ fit, $C=%.2f\times10^{-6}\,\mathrm{mV}^{-2}$" % (C * 1e6))
ax.set_xlabel(r"Tone-1 amplitude $a_1$ (mV)")
ax.set_ylabel(r"Idler-to-target ratio $P_m / P_{f_2/2}$")
ax.set_title(r"(a) Scaling law: $P_m/P_{f_2/2} \propto a_1^{%.2f-%.2f}$ (law: $a_1^2$) over $>2$ decades"
             % (measured[0], measured[-1]))
ax.axvline(265, color="0.4", ls="--", lw=1.0, label=r"$a_{1,\mathrm{th}}\approx 265$ mV")
ax.legend(frameon=False, loc="lower right", ncol=2)
ax.grid(True, which="both", alpha=0.2)
ax.text(0.02, 0.96, r"active region: $a_1\geq 25$, $a_2\geq 350$ mV, $P_{f_2/2}>2\times10^{-4}$",
        transform=ax.transAxes, va="top", fontsize=8)

X1, Y1 = np.meshgrid(a1_edges, a2_edges)
im1 = im2 = None
PANELS = {(0, 0): "b", (0, 1): "c", (0, 2): "d", (1, 0): "e", (1, 1): "f", (1, 2): "g"}
pcm_axes, ratio_axes = [], []
for row, f2 in enumerate([4125, 4135]):
    g_pm, g_pred, g_ratio = maps[f2]
    ax = fig.add_subplot(gs[1 + row, 0])
    im1 = ax.pcolormesh(X1, Y1, np.ma.masked_invalid(g_pm), norm=LogNorm(vmin=vmin, vmax=vmax),
                        cmap="inferno", shading="flat", rasterized=True)
    ax.set_facecolor("0.92"); ax.axvline(265, color="cyan", ls="--", lw=0.9, alpha=0.85)
    ax.set_ylabel(r"Tone-2 amplitude $a_2$ (mV)")
    ax.set_title("(%s) $f_2$=%d (D): measured $P_m$" % (PANELS[(row, 0)], f2))

    ax = fig.add_subplot(gs[1 + row, 1])
    im1 = ax.pcolormesh(X1, Y1, np.ma.masked_invalid(g_pred), norm=LogNorm(vmin=vmin, vmax=vmax),
                        cmap="inferno", shading="flat", rasterized=True)
    ax.set_facecolor("0.92"); ax.axvline(265, color="cyan", ls="--", lw=0.9, alpha=0.85)
    ax.set_title("(%s) $f_2$=%d (D): model $C\\,a_1^2 P_{f_2/2}$" % (PANELS[(row, 1)], f2))
    ax.set_xlabel(r"Tone-1 amplitude $a_1$ (mV)")
    pcm_axes += [fig.axes[-1]]

    ax = fig.add_subplot(gs[1 + row, 2])
    im2 = ax.pcolormesh(X1, Y1, np.ma.masked_invalid(g_ratio), vmin=-0.4, vmax=0.4,
                        cmap="coolwarm", shading="flat", rasterized=True)
    ax.set_facecolor("0.92"); ax.axvline(265, color="k", ls="--", lw=0.9, alpha=0.6)
    ax.set_title("(%s) $\\log_{10}(P_m/\\mathrm{model})$" % PANELS[(row, 2)])
    ratio_axes.append(ax)

cb = fig.colorbar(im1, ax=pcm_axes, orientation="horizontal", location="top",
                  pad=0.10, shrink=0.75)
cb.set_label(r"Integrated idler power $P_m$ / model (SA linear units, log; masked invalid)")
cb.ax.xaxis.set_label_position("top"); cb.ax.xaxis.tick_top()
cb2 = fig.colorbar(im2, ax=ratio_axes, orientation="horizontal", location="top",
                   pad=0.10, shrink=0.6)
cb2.set_label(r"$\log_{10}(P_m/\mathrm{model})$ (decades)")

fig.suptitle("Fig 2 (star panels): pairwise idler scaling law and ($a_1$, $a_2$) maps — "
             "$P_m \\propto a_1^2 P_{f_2/2}$ (D and ND timing pairs in (a))", fontsize=11.5)

p.finish(fig, ["a_scaling_law", "bcd_maps_4125", "efg_maps_4135"])

summary = {"scaling": res_scale, "maps": map_res, "source_sha256": src["sha256"][:16]}
with open("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/fig2_summary.json", "w") as fp:
    json.dump(summary, fp, indent=2, default=str)
print(json.dumps(summary, indent=2, default=str))
