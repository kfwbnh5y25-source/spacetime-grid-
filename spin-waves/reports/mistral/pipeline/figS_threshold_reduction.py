# -*- coding: utf-8 -*-
"""
figS_threshold_reduction.py -- Supplementary: D vs ND logistic fits + threshold reduction.

Split off from fig2_interaction v1 per the figure assessment; every v1 number is
carried unchanged (steps 'carried_numbers_from_v1' and 'threshold_sensitivity').

Data : "phase gating/f1_4069_5_a1_a2_f2_sweep.h5"
       key format '{a1}_{a2}_{f1}_{f2}_{phase}_{label}',
       trace[pass,0,:]=frequency (Hz, window at f2/2 +- 1.25 MHz),
       trace[pass,1,:]=baseline-subtracted linear power, pass axis = re-runs (mean).
Timing mapping used throughout: D ~ phase 1.310, ND ~ phase 4.550
(the threshold-reducing configuration is ND = phase 4.55).

Panels:
  (a) gain(dB) vs a2 at f2 = 4111 MHz for D and ND with logistic fits
      g(a2) = c + L / (1 + exp(-k (a2 - a2_0))); midpoint shift a2_0(ND) < a2_0(D).
  (b) threshold-reduction summary: measured crossing reductions vs threshold choice,
      i.e. 1.33 dB (14% drive) at f2 = 4111 (10 dB threshold), 5.4 dB (~46%, best-case
      f2 = 4105), 12.4 dB (~76%) / 12.5 dB (~75%) at the noisy 3/6 dB thresholds.
  (c) (a2, f2) idler gain map, phase 1.31, a1 = 324 (v1 panel a, moved here).
  (d) (a2, f2) phase-contrast map 10log10 P(4.55)/P(1.31) (v1 panel b, moved here).
"""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/pipeline")
from provenance import Provenance  # noqa: E402

H5 = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/f1_4069_5_a1_a2_f2_sweep.h5"

p = Provenance("figS_threshold_reduction")
src = p.record_source(H5, note="phase-gating f1=4069.5 a1-a2-f2 sweep (v1 fig2 source)")
p.step("source_mapping", params={"D": "phase 1.310", "ND": "phase 4.550",
                                 "timing_definition": "phase of the second pump relative to the "
                                 "seeded coherence window of pump 1"},
       note="threshold-reducing configuration is ND (phase 4.55): idler turns on at lower a2",
       result={})
A2S = np.arange(200, 276)
F2S = np.arange(4075, 4151)
A1_MAIN = 324
PH_D, PH_ND = "1.310", "4.550"

def key_metrics(key, trace):
    x = trace[:, 1, :].astype(np.float64).mean(axis=0)      # average re-run passes
    xc = np.maximum(x, 0.0)
    smooth = np.convolve(xc, np.ones(5) / 5.0, mode="same")
    return (float(xc.sum()), float(smooth.max()))

def keep(k):
    parts = k.split("_")
    return parts[0] in ("0", "324")

traces = p.load_traces(H5, select=keep, loader=key_metrics, pass_agg="mean")

i_a2 = {v: i for i, v in enumerate(A2S)}
i_f2 = {v: j for j, v in enumerate(F2S)}
ref_tot = np.full((len(A2S), len(F2S)), np.nan)
ref_pk = np.full_like(ref_tot, np.nan)
dat_tot = {(A1_MAIN, ph): np.full_like(ref_tot, np.nan) for ph in (PH_D, PH_ND)}
dat_pk = {(A1_MAIN, ph): np.full_like(ref_tot, np.nan) for ph in (PH_D, PH_ND)}
for k, (t, pk) in traces.items():
    a1, a2, f1, f2, ph, lb = k.split("_")
    a2i, f2i = int(a2), int(f2)
    if a1 == "0":
        ref_tot[i_a2[a2i], i_f2[f2i]] = t
        ref_pk[i_a2[a2i], i_f2[f2i]] = pk
    elif a1 == "324" and ph in (PH_D, PH_ND):
        dat_tot[(A1_MAIN, ph)][i_a2[a2i], i_f2[f2i]] = t
        dat_pk[(A1_MAIN, ph)][i_a2[a2i], i_f2[f2i]] = pk

p.step("carried_numbers_from_v1",
       params={"note": "v1 fig2 thresholds, unchanged"},
       note="idler in-band gain max 20.53 dB (total metric) at (a2=271, f2=4111) / 26.53 dB "
            "(peak metric) at (271, 4112); phase contrast max +1.81 dB (265, 4105), min -2.38 dB "
            "(263, 4133), peak metric +3.47 dB (239, 4113); crossings at f2=4111, 10 dB gain "
            "threshold: a2 251.20 (D) -> 249.87 (ND) = 1.33 dB = 14.2% of drive amplitude "
            "= 1.77% of the a2 sweep; all-f2 scan (10 dB): max 5.40 dB (~46.2% drive, f2=4105, "
            "total metric) and 5.20 dB (f2=4123, peak metric); at 3 dB/6 dB thresholds reductions "
            "of 12.4-12.5 dB (74-76% drive) appear but are threshold-definition dependent/noisy",
       result={
           "max_gain_dB_total_metric": 20.530803512318307,
           "max_gain_dB_peak_metric": 26.52741622532229,
           "max_contrast_dB": 1.8092374023439437,
           "min_contrast_dB": -2.376139068278339,
           "max_contrast_peak_metric_dB": 3.4682289150997296,
           "f2_4111_thr10dB": {"a2_th_1.31": 251.201712335703, "a2_th_4.55": 249.87091317930455,
                               "reduction_dB": 1.3307991563984558,
                               "lin_amp_reduction_pct": 14.20538253452338,
                               "a2_sweep_frac_reduction_pct": 1.7743988751979412},
           "all_f2_scan_thr10dB_total_max": [5.404494652255494, 4105, 266.390567484728, 260.98607283247253],
           "all_f2_scan_thr10dB_peak_max": [5.203636891595096, 4123, 250.07677563838354, 244.87313874678844],
           "all_f2_scan_thr3dB_total_max": [12.416828432930913, 4107, 239.2423818290515, 226.8255533961206],
           "all_f2_scan_thr6dB_peak_max": [12.511557174083123, 4107, 239.16024233725147, 226.64868516316835],
       })

# ----------------------------------------------------------------- panel (a): logistic fits
F2_SEG = 4111
j = i_f2[F2_SEG]
g_D = 10 * np.log10(dat_tot[(A1_MAIN, PH_D)][:, j] / ref_tot[:, j])
g_ND = 10 * np.log10(dat_tot[(A1_MAIN, PH_ND)][:, j] / ref_tot[:, j])
g_Dp = 10 * np.log10(dat_pk[(A1_MAIN, PH_D)][:, j] / ref_pk[:, j])
g_NDp = 10 * np.log10(dat_pk[(A1_MAIN, PH_ND)][:, j] / ref_pk[:, j])

def logistic(x, c, L, k, x0):
    return c + L / (1.0 + np.exp(-k * (x - x0)))

# fit the onset region only (a2 <= 270): the total-metric curves double-step
# (plateau at ~10 dB then jump to ~20 dB at a2 ~ 265), so a single sigmoid is fit
# to the peak-metric curves on the onset domain.
A2_FIT = A2S[A2S <= 270]
fits = {}
for name, g, gp in [("D_phase1.31", g_D, g_Dp), ("ND_phase4.55", g_ND, g_NDp)]:
    ok = np.isfinite(gp) & (A2S <= 270)
    try:
        p0 = [0.0, 25.0, 0.25, 250.0]
        popt, _ = curve_fit(logistic, A2S[ok].astype(float), gp[ok], p0=p0, maxfev=20000)
        fits[name] = {"c": float(popt[0]), "L": float(popt[1]), "k": float(popt[2]),
                      "a2_0": float(popt[3]), "fit_data": "peak-metric gain, a2<=270"}
    except Exception as e:
        fits[name] = {"error": str(e)}

shift = None
if "a2_0" in fits["D_phase1.31"] and "a2_0" in fits["ND_phase4.55"]:
    shift = fits["D_phase1.31"]["a2_0"] - fits["ND_phase4.55"]["a2_0"]
res_fit = {"f2_MHz": F2_SEG, "metric": "total in-band gain (dB rel. only-f2 reference)",
           "logistic_fits": fits, "midpoint_shift_a2": shift}
p.step("panel_a_D_vs_ND_logistic_fits",
       params={"model": "g(a2) = c + L/(1+exp(-k(a2-a2_0)))", "phase_D": 1.31, "phase_ND": 4.55},
       note="logistic (sigmoid) onset fits to the peak-metric gain curves on the onset domain "
            "(a2 <= 270; total-metric curves double-step, so a single sigmoid is only fit to the "
            "peak metric); the threshold midpoint a2_0(ND) vs a2_0(D) tests phase-gated "
            "threshold reduction; crossings use gain over the only-f2 reference",
       result=res_fit)
p.array("figS_a_gain_curves_f2_4111", np.vstack([A2S, g_D, g_ND, g_Dp, g_NDp]))

# ----------------------------------------------------------------- panel (b): threshold reduction summary
def crossing(avals, g, thr_db):
    idx = np.where(g >= thr_db)[0]
    if len(idx) == 0:
        return None
    i = int(idx[0])
    if i == 0:
        return float(avals[0])
    return float(avals[i - 1] + (thr_db - g[i - 1]) * (avals[i] - avals[i - 1]) / (g[i] - g[i - 1]))

thr_list = [3.0, 6.0, 10.0]
bars = {}
for thr in thr_list:
    t1, t4 = crossing(A2S, g_D, thr), crossing(A2S, g_ND, thr)
    ent = {"thr_dB": thr}
    if t1 is not None and t4 is not None and t4 < t1:
        d = t1 - t4
        ent.update({"a2_th_D": t1, "a2_th_ND": t4, "reduction_dB": d,
                    "lin_amp_reduction_pct": (1 - 10 ** (-d / 20)) * 100})
    else:
        ent.update({"a2_th_D": t1, "a2_th_ND": t4, "reduction_dB": None})
    bars[str(thr)] = ent
p.step("panel_b_threshold_summary",
       params={"criterion": "first a2 where gain(D or ND) >= thr_dB above only-f2 reference",
               "thresholds_dB": thr_list},
       note="f2=4111 local values; global best-case numbers carried from v1 (5.40 dB total / "
            "5.20 dB peak at 10 dB; 12.4-12.5 dB at 3/6 dB, noisy)",
       result=bars)
p.array("figS_b_threshold_bars", np.array(
    [[bars[s].get("reduction_dB") or np.nan for s in bars]]))

# ----------------------------------------------------------------- panels (c)+(d): moved v1 maps
gain_tot_db = 10 * np.log10(dat_tot[(A1_MAIN, PH_D)] / ref_tot)
contrast_tot_db = 10 * np.log10(dat_tot[(A1_MAIN, PH_ND)] / dat_tot[(A1_MAIN, PH_D)])
p.array("figS_c_gain_map", gain_tot_db)
p.array("figS_d_contrast_map", contrast_tot_db)
p.step("panels_cd_moved_v1_maps",
       params={"panel_c": "gain map 10log10(P_D/P_ref)", "panel_d": "contrast map 10log10(P_ND/P_D)"},
       note="moved verbatim from fig2 v1; max gain 20.53 dB at (271, 4111); contrast +1.81/-2.38 dB",
       result={})

# ----------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9})
fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.2))

ax = axes[0, 0]
ax.plot(A2S, g_D, "o-", ms=3, color="tab:blue", label="D (phase 1.31), total metric")
ax.plot(A2S, g_ND, "s-", ms=3, color="tab:red", label="ND (phase 4.55), total metric")
ax.plot(A2S, g_Dp, "o", ms=2, mfc="none", color="tab:blue", alpha=0.5, label="D, peak metric")
ax.plot(A2S, g_NDp, "s", ms=2, mfc="none", color="tab:red", alpha=0.5, label="ND, peak metric")
xx = np.linspace(A2S[0], A2S[-1], 400)
for name, col in [("D_phase1.31", "tab:blue"), ("ND_phase4.55", "tab:red")]:
    f = fits[name]
    if "a2_0" in f:
        ax.plot(xx, logistic(xx, f["c"], f["L"], f["k"], f["a2_0"]), color=col,
                lw=1.8, alpha=0.8,
                label="logistic fit %s ($a_{2,0}$=%.1f, $k$=%.2f)" % (name.split("_")[0], f["a2_0"], f["k"]))
        ax.axvline(f["a2_0"], color=col, lw=0.9, alpha=0.5)
ax.set_xlabel("Pump-2 amplitude $a_2$ (key units)")
ax.set_ylabel("Idler in-band power gain (dB rel. only-$f_2$ ref)")
ax.set_title("(a) D vs ND idler onset, logistic fits ($f_2$ = %d MHz)" % F2_SEG)
ax.legend(fontsize=7, loc="upper left")
ax.text(0.30, 0.10, "sigmoid midpoints: $a_{2,0}$(D) = %.1f, $a_{2,0}$(ND) = %.1f"
        % (fits["D_phase1.31"].get("a2_0", np.nan), fits["ND_phase4.55"].get("a2_0", np.nan)),
        transform=ax.transAxes, fontsize=8.5,
        bbox=dict(fc="w", alpha=0.85, ec="0.5"))

ax = axes[0, 1]
xs = np.array([0.0, 1.0, 2.0, 3.6])
w = 0.38
vals_tot = [bars[str(t)].get("reduction_dB") if bars[str(t)].get("reduction_dB") is not None else 0.0
            for t in thr_list]
ax.bar(xs[:3] - w / 2, vals_tot, w, color="tab:purple", label="this file, $f_2$ = 4111 (total metric)")
max10_tot, max10_pk = 5.404494652255494, 5.203636891595096
ax.bar([xs[3] - w / 2], [max10_tot], w, color="tab:orange",
       label="best-$f_2$ scan max (total): 5.40 dB $\\approx$ 46% drive ($f_2$ = 4105)")
ax.hlines(max10_pk, xs[3] - w, xs[3], color="k", ls=":",
          label="best-$f_2$ scan max (peak): 5.20 dB")
ax.annotate("10 dB threshold, $f_2$=4111:\n%.2f dB $\\approx$ %.0f%% of drive\n(best case %.1f dB $\\approx$ 46%%)\n"
            "FLAG: 50%% not reproduced as stated;\n3/6 dB thresholds hit 12.4-12.5 dB\n(74-76%%) but are noisy"
            % (bars["10.0"]["reduction_dB"], bars["10.0"]["lin_amp_reduction_pct"], max10_tot),
            xy=(0.98, 0.98), xycoords="axes fraction", ha="right", va="top", fontsize=8,
            bbox=dict(fc="w", alpha=0.9, ec="0.5"))
ax.set_xticks(xs)
ax.set_xticklabels(["3 dB", "6 dB", "10 dB", "best-$f_2$\nscan max"])
ax.set_ylabel("Measured threshold reduction (dB of $a_2$)")
ax.set_title("(b) Threshold reduction vs threshold choice, $a_1$ = 324")
ax.legend(fontsize=7.5, loc="lower left")

X, Y = np.meshgrid(F2S, A2S)
ax = axes[1, 0]
im = ax.pcolormesh(X, Y, gain_tot_db, shading="auto", cmap="viridis", rasterized=True)
cb = fig.colorbar(im, ax=ax, orientation="horizontal", location="top", pad=0.10, shrink=0.8)
cb.set_label("Idler in-band power gain (dB rel. only-$f_2$ ref), phase 1.31")
cb.ax.xaxis.set_label_position("top"); cb.ax.xaxis.tick_top()
ax.set_xlabel("Pump-2 frequency $f_2$ (MHz)")
ax.set_ylabel("Pump-2 amplitude $a_2$ (key units)")
ax.set_title("(c) Idler gain map, $a_1$ = 324 (moved from Fig 2)")
ax.plot(4111, 271, "r*", ms=13, mec="w")

ax = axes[1, 1]
vmax = np.nanmax(np.abs(contrast_tot_db))
im = ax.pcolormesh(X, Y, contrast_tot_db, shading="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   rasterized=True)
cb = fig.colorbar(im, ax=ax, orientation="horizontal", location="top", pad=0.10, shrink=0.8)
cb.set_label("Contrast $10\\log_{10}[P(4.55)/P(1.31)]$ (dB)")
cb.ax.xaxis.set_label_position("top"); cb.ax.xaxis.tick_top()
ax.set_xlabel("Pump-2 frequency $f_2$ (MHz)")
ax.set_ylabel("Pump-2 amplitude $a_2$ (key units)")
ax.set_title("(d) Phase-contrast map (moved from Fig 2); +1.81 / $-$2.38 dB extremes")

fig.suptitle("Fig S (threshold reduction): D vs ND onset fits and idler threshold reduction, "
             "phase-gating sweep $f_1$ = 4069.5 MHz", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.96])

p.finish(fig, ["a_D_vs_ND_logistic", "b_threshold_reduction", "c_gain_map", "d_contrast_map"])

summary = {"fits": res_fit, "bars": bars, "source_sha256": src["sha256"][:16]}
with open("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/figS_summary.json", "w") as fp:
    json.dump(summary, fp, indent=2, default=str)
print(json.dumps(summary, indent=2, default=str))
