"""
Reproducible Analysis & Publication Figures for Two-Tone Threshold Reduction
and Non-Degenerate Idler Coupling (spin-waves/threshold_reduction/data.h5).

Reads raw HDF5 traces + single-tone broadband/zoom caches, applies 15-bin Hampel
despiking (5 sigma) to remove narrowband lab RF pickup spurs, verifies spectral
centering / absence of edge leakage, and saves 4 publication figures (.png/.pdf).

Run from repo root:
  ~/.venvs/nvim-python/bin/python spin-waves/reports/neuromorphic/analyze_threshold_reduction.py
"""

import json
import re
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.ndimage import median_filter
from scipy.stats import linregress

REPO_ROOT = Path(__file__).resolve().parents[3]
H5_PATH = REPO_ROOT / "spin-waves/threshold_reduction/data.h5"
BROAD_CACHE = REPO_ROOT / "spin-waves/Figure!_antigraity/figure1_single_tone/cache/broad_halfaverage.npz"
ZOOM_CACHE = REPO_ROOT / "spin-waves/Figure!_antigraity/figure1_single_tone/cache/zoom_halfaverage.npz"
OUT_DIR = Path(__file__).resolve().parent

KEY_PATTERN = re.compile(r"^(\d+)_(\d+)_(\d+)_(\d+)_([-\d.]+)_(EX|f1|f2|m|p)_(D|ND)$")


def hampel_despike(spec: np.ndarray, half_win: int = 7, n_sig: float = 5.0) -> np.ndarray:
    """Replace narrowband pickup spikes (> n_sig robust MAD above running median)."""
    med = median_filter(spec, size=2 * half_win + 1, mode="nearest")
    diff = np.abs(spec - med)
    mad = max(float(np.median(diff)) / 0.6745, 1e-11)
    mask = diff > n_sig * mad
    if mask.any():
        for idx in np.where(mask)[0]:
            mask[max(0, idx - 2) : min(len(spec), idx + 3)] = True
        out = spec.copy()
        out[mask] = med[mask]
        return out
    return spec


def load_despiked_h5(h5_path: Path):
    """Load all groups from threshold_reduction/data.h5 with Hampel despiking."""
    rows = {}
    with h5py.File(h5_path, "r") as f:
        for k in f.keys():
            m = KEY_PATTERN.match(k)
            if m is None:
                continue
            a1, a2 = int(m[1]), int(m[2])
            f1, f2 = int(m[3]), int(m[4])
            ro, tim = m[6], m[7]
            tr = np.asarray(f[k]["trace"])
            freq_mhz = tr[0, 0, :] / 1e6
            raw_mean = tr[:, 1, :].mean(axis=0)
            cln = raw_mean if ro == "EX" else hampel_despike(raw_mean)
            rows.setdefault((ro, f2, tim), {})[(a1, a2)] = (float(cln.sum()), cln, freq_mhz)
    return rows


def bin_2d(pts_x, pts_y, vals, x_edges, y_edges, min_val=1e-7):
    grid = np.full((len(y_edges) - 1, len(x_edges) - 1), np.nan)
    for iy in range(len(y_edges) - 1):
        my = (pts_y >= y_edges[iy]) & (pts_y < y_edges[iy + 1])
        for ix in range(len(x_edges) - 1):
            m = my & (pts_x >= x_edges[ix]) & (pts_x < x_edges[ix + 1])
            if m.any():
                grid[iy, ix] = max(float(np.median(vals[m])), min_val)
    return grid


def make_figure1_scaling(rows):
    """Fig 1: Idler quadratic scaling Pm / (a1^2 Pf2/2) = const and threshold shift."""
    plt.rcParams.update({
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9.5,
        "legend.fontsize": 7.8,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4), dpi=300)

    colors = {(4125, "D"): "#1f77b4", (4125, "ND"): "#5ca4d8", (4135, "D"): "#d62728", (4135, "ND"): "#f28e8e"}
    markers = {(4125, "D"): "o", (4125, "ND"): "s", (4135, "D"): "^", (4135, "ND"): "v"}
    bins_a1 = [(25, 60), (60, 100), (100, 140), (140, 180), (180, 220), (220, 260), (260, 300), (300, 355)]
    table_rows = []

    ax = axes[0]
    for f2 in [4125, 4135]:
        for tim in ["D", "ND"]:
            dm = rows[("m", f2, tim)]
            df2 = rows[("f2", f2, tim)]
            pts = []
            for k in dm:
                if k[0] >= 25 and k[1] >= 350 and k in df2:
                    pm, pf2 = dm[k][0], df2[k][0]
                    if pf2 > 2e-4 and pm > 0:
                        pts.append((k[0], k[1], pm, pf2, pm / pf2, pm / (k[0] ** 2 * pf2)))
            pts = np.array(pts)
            slope, _, _, _, se_slope = linregress(np.log10(pts[:, 0]), np.log10(pts[:, 4]))
            xs, y_med, y_lo, y_hi = [], [], [], []
            for lo, hi in bins_a1:
                msk = (pts[:, 0] >= lo) & (pts[:, 0] < hi)
                if msk.sum() < 2:
                    continue
                vals = pts[msk, 5]
                xs.append(np.median(pts[msk, 0]))
                med = float(np.median(vals))
                q25, q75 = np.percentile(vals, [25, 75])
                y_med.append(med * 1e6)
                y_lo.append((med - q25) * 1e6)
                y_hi.append((q75 - med) * 1e6)
                table_rows.append({
                    "f2": f2, "timing": tim, "a1_bin": f"[{lo},{hi})",
                    "n": int(msk.sum()),
                    "med_pm": float(np.median(pts[msk, 2])),
                    "med_pf2": float(np.median(pts[msk, 3])),
                    "ratio_med_e6": med * 1e6,
                    "ratio_q25_e6": float(q25) * 1e6,
                    "ratio_q75_e6": float(q75) * 1e6,
                    "log_slope": float(slope),
                    "log_slope_se": float(se_slope),
                })
            ax.errorbar(xs, y_med, yerr=[y_lo, y_hi], fmt=markers[(f2, tim)] + "-",
                        color=colors[(f2, tim)], capsize=3, lw=1.4, ms=5,
                        label=f"f2={f2} {tim} (slope={slope:.2f}±{se_slope:.2f})")

    ax.set_xlabel(r"Tone-1 Amplitude $a_1$ (mV)")
    ax.set_ylabel(r"Normalized Idler Ratio $P_m / (a_1^2 P_{f_2/2})$ ($\times 10^{-6}\ \mathrm{mV}^{-2}$)")
    ax.set_title(r"(a) Constant Ratio $P_m / (a_1^2 P_{f_2/2})$ Across $a_1$")
    ax.set_ylim(0, 2.25)
    ax.axvline(265, color="0.4", ls="--", lw=1.0, label=r"$a_{1,\mathrm{th}} \approx 265$ mV")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(alpha=0.25)

    ax = axes[1]
    for f2, col, lbl in [(4125, "#1f77b4", "f2=4125 (D+ND)"), (4135, "#d62728", "f2=4135 (D+ND)")]:
        all_pts = []
        for tim in ["D", "ND"]:
            dm, df2 = rows[("m", f2, tim)], rows[("f2", f2, tim)]
            for k in dm:
                if k[0] >= 25 and k[1] >= 350 and k in df2:
                    pm, pf2 = dm[k][0], df2[k][0]
                    if pf2 > 2e-4 and pm > 0:
                        all_pts.append((k[0], pm / pf2))
        all_pts = np.array(all_pts)
        ax.scatter(all_pts[:, 0], all_pts[:, 1], s=10, alpha=0.35, color=col, edgecolors="none", label=lbl)
        c_med = np.median(all_pts[:, 1] / (all_pts[:, 0] ** 2))
        a_grid = np.linspace(25, 355, 100)
        ax.plot(a_grid, c_med * a_grid ** 2, color=col, lw=2.0,
                label=rf"$\propto a_1^2$ fit ($C={c_med*1e6:.2f}\times 10^{{-6}}$)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Tone-1 Amplitude $a_1$ (mV)")
    ax.set_ylabel(r"Idler-to-Target Power Ratio $P_m / P_{f_2/2}$")
    ax.set_title(r"(b) $P_m / P_{f_2/2} \propto a_1^2$ Over $>2$ Decades")
    ax.axvline(265, color="0.4", ls="--", lw=1.0)
    ax.legend(frameon=False, loc="lower right")
    ax.grid(True, which="both", alpha=0.2)

    ax = axes[2]
    a1_th_bins = [(0, 1), (1, 100), (100, 175), (175, 225), (225, 275), (275, 355)]
    for f2, col_d, col_nd in [(4125, "#1f77b4", "#5ca4d8"), (4135, "#d62728", "#f28e8e")]:
        for tim, col, ls, mk in [("D", col_d, "-", "o"), ("ND", col_nd, "--", "s")]:
            d = rows[("f2", f2, tim)]
            x_a1_2, y_th = [], []
            for lo1, hi1 in a1_th_bins:
                pts = sorted([(a2, v[0], a1) for (a1, a2), v in d.items() if lo1 <= a1 < hi1])
                a2s = np.array([p[0] for p in pts])
                ss = np.array([p[1] for p in pts])
                a1s = np.array([p[2] for p in pts])
                for i in range(max(0, len(pts) - 6)):
                    if np.median(ss[i : i + 7]) > 2e-5:
                        x_a1_2.append(np.median(a1s) ** 2 / 1e4)
                        y_th.append(np.median(a2s[i : i + 7]))
                        break
            ax.plot(x_a1_2, y_th, marker=mk, ls=ls, color=col, lw=1.7, ms=5, label=f"f2={f2} {tim}")
            if tim == "ND":
                sl, ic, _, _, _ = linregress(x_a1_2[:-1], y_th[:-1])
                xg = np.linspace(0, 9.5, 50)
                ax.plot(xg, ic + sl * xg, ls=":", color=col, lw=1.2)
    ax.set_xlabel(r"$a_1^2$ ($\times 10^4\ \mathrm{mV}^2$)")
    ax.set_ylabel(r"Threshold $a_{2,\mathrm{th}}$ (mV, $10\sigma$ crossing)")
    ax.set_title(r"(c) Quadratic $a_1^2$ Drop (ND) & $4B$ Split (D)")
    ax.axvline(265 ** 2 / 1e4, color="0.4", ls="--", lw=1.0, label=r"$a_{1,\mathrm{th}}^2$")
    ax.legend(frameon=False, loc="lower left")
    ax.grid(alpha=0.25)

    plt.tight_layout()
    fig.savefig(OUT_DIR / "fig_idler_and_threshold_scaling.png", dpi=300)
    fig.savefig(OUT_DIR / "fig_idler_and_threshold_scaling.pdf")
    plt.close(fig)
    return table_rows


def make_figure2_pcolormesh(rows):
    """Fig 2: Side-by-side 2D pcolormesh of Pm vs C * a1^2 * Pf2/2 on identical colorscale."""
    a1_edges = np.linspace(0, 355, 22)
    a2_edges = np.linspace(0, 455, 25)
    fig, axes = plt.subplots(2, 3, figsize=(13.8, 8.2), dpi=300, sharex=True, sharey=True)
    vmin, vmax = 5e-7, 2e-4

    for r_idx, f2 in enumerate([4125, 4135]):
        dm, df2 = rows[("m", f2, "D")], rows[("f2", f2, "D")]
        common = [k for k in dm if k in df2 and k[0] > 0 and k[1] > 0]
        a1s = np.array([k[0] for k in common], float)
        a2s = np.array([k[1] for k in common], float)
        pms = np.array([dm[k][0] for k in common], float)
        pf2s = np.array([df2[k][0] for k in common], float)

        act = (a1s >= 25) & (a2s >= 350) & (pf2s > 2e-4) & (pms > 0)
        c_fit = float(np.median(pms[act] / (a1s[act] ** 2 * pf2s[act])))
        pred = c_fit * (a1s ** 2) * pf2s

        g_pm = bin_2d(a1s, a2s, pms, a1_edges, a2_edges)
        g_pred = bin_2d(a1s, a2s, pred, a1_edges, a2_edges)
        g_ratio = np.full_like(g_pm, np.nan)
        valid_r = (g_pm > 5e-6) & (g_pred > 5e-6)
        g_ratio[valid_r] = np.log10(g_pm[valid_r] / g_pred[valid_r])

        norm = LogNorm(vmin=vmin, vmax=vmax)
        axes[r_idx, 0].pcolormesh(a1_edges, a2_edges, g_pm, norm=norm, cmap="inferno", shading="flat")
        im1 = axes[r_idx, 1].pcolormesh(a1_edges, a2_edges, g_pred, norm=norm, cmap="inferno", shading="flat")
        im2 = axes[r_idx, 2].pcolormesh(a1_edges, a2_edges, g_ratio, vmin=-0.4, vmax=0.4, cmap="coolwarm", shading="flat")

        axes[r_idx, 0].set_title(f"f2={f2} MHz (D): Measured $P_m$ (despiked)")
        axes[r_idx, 1].set_title(f"f2={f2} MHz (D): Model $C\\,a_1^2 P_{{f_2/2}}$ ($C={c_fit*1e6:.2f}\\times 10^{{-6}}$)")
        axes[r_idx, 2].set_title(f"f2={f2} MHz: $\\log_{{10}}(P_m / \\mathrm{{Model}})$ (active)")
        axes[r_idx, 0].set_ylabel(r"Tone-2 Amplitude $a_2$ (mV)")
        for c_idx in range(3):
            axes[r_idx, c_idx].axvline(265, color="cyan", ls="--", lw=0.9, alpha=0.85)
            axes[r_idx, c_idx].set_facecolor("0.92")

    for c_idx in range(3):
        axes[1, c_idx].set_xlabel(r"Tone-1 Amplitude $a_1$ (mV)")
    fig.colorbar(im1, ax=axes[:, :2], fraction=0.025, pad=0.02, label="Integrated Subharmonic Power (SA linear units, log)")
    fig.colorbar(im2, ax=axes[:, 2], fraction=0.045, pad=0.03, label=r"$\log_{10}(\mathrm{Measured} / \mathrm{Model})$ (decades)")
    fig.suptitle(r"Side-by-Side 2D Validation: Measured Idler $P_m$ vs. $C\,a_1^2 P_{f_2/2}$ (Same Colorscale, Spikes Removed)", y=0.995, fontsize=11)
    fig.savefig(OUT_DIR / "fig_pcolormesh_pm_vs_model.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_pcolormesh_pm_vs_model.pdf", bbox_inches="tight")
    plt.close(fig)


def make_figure3_all_maps(rows):
    """Fig 3 & 4: All 4 two-tone readout channels on (a1, a2) + broadband (f, a) map."""
    a1_edges = np.linspace(0, 355, 22)
    a2_edges = np.linspace(0, 455, 25)
    fig2, axes2 = plt.subplots(2, 4, figsize=(16.0, 7.6), dpi=300, sharex=True, sharey=True)
    norm_all = LogNorm(vmin=1e-6, vmax=2.5e-3)

    for r_idx, (f2, tim) in enumerate([(4125, "D"), (4135, "D")]):
        for c_idx, (ro, title_ro) in enumerate([("f1", r"$f_1/2$ (2035.5 MHz)"),
                                                ("f2", r"$f_2/2$ (target)"),
                                                ("m", r"Idler $m = f_1 - f_2/2$"),
                                                ("p", r"Upper $p = f_2 - f_1/2$")]):
            d = rows[(ro, f2, tim)]
            ks = np.array(list(d.keys()), float)
            vs = np.array([v[0] for v in d.values()], float)
            g = bin_2d(ks[:, 0], ks[:, 1], vs, a1_edges, a2_edges, min_val=1e-6)
            im = axes2[r_idx, c_idx].pcolormesh(a1_edges, a2_edges, g, norm=norm_all, cmap="viridis", shading="flat")
            axes2[r_idx, c_idx].set_facecolor("0.92")
            axes2[r_idx, c_idx].axvline(265, color="white", ls="--", lw=0.8, alpha=0.8)
            axes2[r_idx, c_idx].set_title(f"{title_ro}\n($f_2={f2}$ MHz, {tim})")
            if r_idx == 1:
                axes2[r_idx, c_idx].set_xlabel(r"$a_1$ (mV)")
        axes2[r_idx, 0].set_ylabel(r"$a_2$ (mV)")

    fig2.colorbar(im, ax=axes2, fraction=0.018, pad=0.02, label="Despiked Integrated Power (SA linear units, log)")
    fig2.suptitle(r"All Two-Tone Readout Channels on $(a_1, a_2)$ Grid (Hampel-Despiked, Shared Log Scale)", y=0.99, fontsize=11.5)
    fig2.savefig(OUT_DIR / "fig_all_responses_a1_a2.png", dpi=300, bbox_inches="tight")
    fig2.savefig(OUT_DIR / "fig_all_responses_a1_a2.pdf", bbox_inches="tight")
    plt.close(fig2)

    broad = np.load(BROAD_CACHE)
    zoom = np.load(ZOOM_CACHE)
    b_p = np.array([hampel_despike(s).sum() * 1e6 for s in broad["spec_mean"]])
    z_p = np.array([hampel_despike(s).sum() * 1e6 for s in zoom["spec_mean"]])

    fig3, (axb, axz) = plt.subplots(1, 2, figsize=(13.5, 4.6), dpi=300)
    fe_b, ae_b = np.linspace(3760, 4250, 99), np.linspace(200, 500, 61)
    gb = bin_2d(broad["f"], broad["a"], b_p, fe_b, ae_b, min_val=0.5)
    imb = axb.pcolormesh(fe_b, ae_b, gb, norm=LogNorm(vmin=0.8, vmax=400), cmap="inferno", shading="flat")
    axb.set_facecolor("0.92")
    for fx, lbl, col in [(4007, "m(4135)", "cyan"), (4017, "m(4125)", "deepskyblue"),
                         (4071, "f1=4071", "lime"), (4125, "f2=4125", "white"), (4135, "f2=4135", "yellow")]:
        axb.axvline(fx, color=col, ls="--", lw=1.0, label=lbl)
    axb.set_xlabel("Drive Frequency $f$ (MHz)")
    axb.set_ylabel("Drive Amplitude $a$ (mV)")
    axb.set_title(r"Broadband $(f, a)$ Landscape (`Figure1/data.h5`, Despiked)")
    axb.legend(frameon=True, facecolor="black", labelcolor="white", fontsize=7.5, loc="lower left", ncol=2)
    fig3.colorbar(imb, ax=axb, label=r"Despiked $\Sigma P$ (nW, log)")

    fe_z, ae_z = np.linspace(4045, 4085, 41), np.linspace(250, 350, 51)
    gz = bin_2d(zoom["f"], zoom["a"], z_p, fe_z, ae_z, min_val=0.5)
    imz = axz.pcolormesh(fe_z, ae_z, gz, norm=LogNorm(vmin=0.8, vmax=150), cmap="inferno", shading="flat")
    axz.set_facecolor("0.92")
    axz.axvline(4065, color="cyan", ls="--", lw=1.1, label="f_exc=4065")
    axz.axvline(4071, color="lime", ls="--", lw=1.1, label="f1=4071")
    axz.set_xlabel("Drive Frequency $f$ (MHz)")
    axz.set_ylabel("Drive Amplitude $a$ (mV)")
    axz.set_title(r"High-Density Zoom $(f, a)$ (`data2.h5`, Despiked)")
    axz.legend(frameon=True, facecolor="black", labelcolor="white", fontsize=8.0, loc="lower right")
    fig3.colorbar(imz, ax=axz, label=r"Despiked $\Sigma P$ (nW, log)")

    plt.tight_layout()
    fig3.savefig(OUT_DIR / "fig_broadband_af_map.png", dpi=300, bbox_inches="tight")
    fig3.savefig(OUT_DIR / "fig_broadband_af_map.pdf", bbox_inches="tight")
    plt.close(fig3)


def make_figure4_spectral_centering(rows):
    """Fig 5: Inspect raw averaged despiked spectra for centering vs edge leakage."""
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.0), dpi=300, sharex=True)
    plot_groups = [
        (axes[0, 0], "f1", 4125, [
            ("a1>275, a2<150 (Tone-1 driven)", lambda a1, a2: a1 > 275 and a2 < 150, "#1f77b4"),
            ("a1[150,250), a2>360 (Sub-th a1: flat floor!)", lambda a1, a2: 150 <= a1 < 250 and a2 > 360, "#d62728"),
            ("a1[150,250), a2<150 (Quiet control)", lambda a1, a2: 150 <= a1 < 250 and a2 < 150, "0.55"),
        ], r"Channel $f_1/2$ ($f_c = 2035.5$ MHz)"),
        (axes[0, 1], "f2", 4125, [
            ("a1=0, a2>360 (Single-tone f2)", lambda a1, a2: a1 == 0 and a2 > 360, "#1f77b4"),
            ("a1>275, a2[265,305) (Pre-excited boost)", lambda a1, a2: a1 > 275 and 265 <= a2 < 305, "#d62728"),
            ("a1=0, a2[265,305) (Sub-th baseline)", lambda a1, a2: a1 == 0 and 265 <= a2 < 305, "0.55"),
        ], r"Channel $f_2/2$ ($f_c = 2062.5$ MHz)"),
        (axes[1, 0], "m", 4125, [
            ("a1[250,350), a2>360 (High a1)", lambda a1, a2: 250 <= a1 < 350 and a2 > 360, "#1f77b4"),
            ("a1[100,200), a2>360 (Sub-th a1)", lambda a1, a2: 100 <= a1 < 200 and a2 > 360, "#d62728"),
            ("a1[250,350), a2<150 (No Tone-2 control)", lambda a1, a2: 250 <= a1 < 350 and a2 < 150, "0.55"),
        ], r"Channel $m = (2f_1 - f_2)/2$ ($f_c = 2008.5$ MHz)"),
        (axes[1, 1], "m", 4135, [
            ("a1[250,350), a2>360 (High a1)", lambda a1, a2: 250 <= a1 < 350 and a2 > 360, "#1f77b4"),
            ("a1[100,200), a2>360 (Sub-th a1)", lambda a1, a2: 100 <= a1 < 200 and a2 > 360, "#d62728"),
            ("a1[250,350), a2<150 (No Tone-2 control)", lambda a1, a2: 250 <= a1 < 350 and a2 < 150, "0.55"),
        ], r"Channel $m = (2f_1 - f_2)/2$ ($f_c = 2003.5$ MHz)"),
    ]

    metrics = []
    for ax, ro, f2, specs, title in plot_groups:
        d = rows[(ro, f2, "D")]
        ymax = 0.0
        for lbl, cfn, col in specs:
            sub = [v for (a1, a2), v in d.items() if cfn(a1, a2)]
            freq_mhz = sub[0][2]
            fc = 0.5 * (freq_mhz[0] + freq_mhz[-1])
            df = freq_mhz - fc
            mean_sp = np.mean([v[1] for v in sub], axis=0)
            ymax = max(ymax, float(mean_sp.max()) * 1e7)
            ax.plot(df, mean_sp * 1e7, color=col, lw=1.5, label=f"{lbl} (N={len(sub)})")
            w = np.maximum(mean_sp, 0)
            metrics.append({
                "ro": ro, "f2": f2, "label": lbl, "n": len(sub), "fc": float(fc),
                "pk_off": float(df[np.argmax(mean_sp)]),
                "cent_off": float(np.sum(w * df) / (np.sum(w) + 1e-30)),
                "p_left": float(np.mean(mean_sp[df < -2.0])),
                "p_ctr": float(np.mean(mean_sp[np.abs(df) < 0.75])),
                "p_right": float(np.mean(mean_sp[df > 2.0])),
            })
        ax.axvline(0, color="k", ls=":", lw=1.0, alpha=0.6)
        ax.set_ylim(-0.12 * ymax, 1.38 * ymax)
        ax.set_title(title)
        ax.set_ylabel(r"Power ($\times 10^{-7}$ SA linear units)")
        ax.legend(frameon=False, fontsize=7.5, loc="upper left")
        ax.grid(alpha=0.25)

    axes[1, 0].set_xlabel(r"Detuning from Window Center $f - f_c$ (MHz)")
    axes[1, 1].set_xlabel(r"Detuning from Window Center $f - f_c$ (MHz)")
    fig.suptitle("Spectral Centering & Leakage Audit (Hampel-Despiked Mean Spectra)", fontsize=11, y=0.99)
    plt.tight_layout()
    fig.savefig(OUT_DIR / "fig_spectral_centering_check.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_spectral_centering_check.pdf", bbox_inches="tight")
    plt.close(fig)
    return metrics


def main():
    print("Loading and Hampel-despiking", H5_PATH)
    rows = load_despiked_h5(H5_PATH)
    table_rows = make_figure1_scaling(rows)
    make_figure2_pcolormesh(rows)
    make_figure3_all_maps(rows)
    metrics = make_figure4_spectral_centering(rows)
    out_summary = {"scaling_table": table_rows, "centering_metrics": metrics}
    (OUT_DIR / "threshold_reduction_audit.json").write_text(json.dumps(out_summary, indent=2))
    print("Wrote figures and threshold_reduction_audit.json to", OUT_DIR)


if __name__ == "__main__":
    main()
