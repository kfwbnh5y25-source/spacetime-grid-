"""
Generates unbinned, Hampel-despiked pcolormesh-interpolated 2D (a1, a2) maps
for ALL 4 (f1, f2) combinations in spin-waves/threshold_reduction/data.h5 (72,920 groups):
  1. fig_all_maps_4pairs.{png,pdf}: 4 (f1,f2) pairs x 4 channels (f1/2, f2/2, m, p) [D]
  2. fig_pcolormesh_pm_vs_model_all4.{png,pdf}: 4 (f1,f2) pairs x 3 cols (Measured P_m, Model C*a1^2*P_f2, Log10 Residual)
"""

from pathlib import Path
import h5py
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = Path(__file__).resolve().parent
TR_DIR = HERE.parent.parent / "threshold_reduction"
H5_TR = TR_DIR / "data.h5"


def hampel_batch(specs: np.ndarray, hw: int = 7, nsig: float = 5.0) -> np.ndarray:
    pad = np.pad(specs, ((0, 0), (hw, hw)), mode="edge")
    win = sliding_window_view(pad, 2 * hw + 1, axis=1)
    med = np.median(win, axis=2)
    mad = np.median(np.abs(win - med[:, :, None]), axis=2)
    sig = 1.4826 * mad
    mask = (sig > 0) & (np.abs(specs - med) > nsig * sig)
    out = specs.copy()
    out[mask] = med[mask]
    return out


def gaussian_pcolormesh_interp(pts_x: np.ndarray, pts_y: np.ndarray, vals: np.ndarray,
                               x_grid: np.ndarray, y_grid: np.ndarray,
                               sig_x: float = 14.0, sig_y: float = 14.0) -> np.ndarray:
    """
    Smooth kernel-weighted interpolation onto a dense regular (x_grid, y_grid) mesh
    for pcolormesh without dropping or binning any Sobol data points.
    """
    Xg, Yg = np.meshgrid(x_grid, y_grid, indexing="xy")  # (ny, nx)
    dx = (Xg[:, :, None] - pts_x[None, None, :]) / sig_x
    dy = (Yg[:, :, None] - pts_y[None, None, :]) / sig_y
    w = np.exp(-0.5 * (dx * dx + dy * dy))
    w_sum = np.sum(w, axis=2)
    z = np.sum(w * vals[None, None, :], axis=2) / np.maximum(w_sum, 1e-12)
    z[w_sum < 1e-3] = np.nan
    return z


def main():
    print(f"Loading {H5_TR}...")
    meta, specs = [], []
    with h5py.File(H5_TR, "r") as f:
        for k in f.keys():
            p = k.split("_")
            if len(p) < 7 or p[5] == "EX":
                continue
            tr = f[k]["trace"][:]
            s = tr[:, 1, :].mean(0) if tr.ndim == 3 else tr[1, :]
            meta.append((int(p[0]), int(p[1]), int(p[2]), int(p[3]), p[5], p[6]))
            specs.append(s)

    clean = hampel_batch(np.array(specs, dtype=np.float64))
    raw_sum = np.sum(clean, axis=1)

    data = {}
    for idx, (a1, a2, f1, f2, ro, tim) in enumerate(meta):
        data.setdefault((f1, f2, ro, tim), {})[(a1, a2)] = float(raw_sum[idx])

    pairs = [(4071, 4125), (4071, 4135), (4081, 4125), (4081, 4135)]
    xg = np.linspace(0, 355, 90)
    yg = np.linspace(0, 455, 90)
    norm_all = LogNorm(vmin=1e-6, vmax=5e-3)
    norm_pm = LogNorm(vmin=1e-6, vmax=5e-4)
    channels = [("f1", "$f_1/2$"), ("f2", "$f_2/2$"), ("m", "Idler $m$"), ("p", "Upper $p$")]

    for tim in ["D", "ND"]:
        # --- Figure A: All 4 (f1, f2) pairs x 4 channels (f1/2, f2/2, m, p) ---
        fig, axes = plt.subplots(4, 4, figsize=(16, 14), dpi=180, sharex=True, sharey=True)
        fig.suptitle(f"Complete Despiked 2D ($a_1, a_2$) Landscapes Across All 4 $(f_1, f_2)$ Pairs ({tim}, Unbinned pcolormesh)",
                     fontsize=13, fontweight="bold", y=0.995)

        for row, (f1, f2) in enumerate(pairs):
            for col, (ro, ch_lbl) in enumerate(channels):
                ax = axes[row, col]
                pts = data[(f1, f2, ro, tim)]
                px = np.array([k[0] for k in pts.keys()], dtype=float)
                py = np.array([k[1] for k in pts.keys()], dtype=float)
                pv = np.maximum(np.array(list(pts.values()), dtype=float), 1e-6)
                Z = gaussian_pcolormesh_interp(px, py, pv, xg, yg)
                im = ax.pcolormesh(xg, yg, Z, shading="gouraud", norm=norm_all, cmap="inferno")
                if row == 0:
                    ax.set_title(ch_lbl, fontsize=11, fontweight="bold")
                if col == 0:
                    ax.set_ylabel(f"({f1}, {f2}) MHz\n$a_2$ (mV)", fontsize=9.5)
                if row == 3:
                    ax.set_xlabel("$a_1$ (mV)")

        cbar = fig.colorbar(im, ax=axes, shrink=0.6, pad=0.02)
        cbar.set_label("Hampel-Despiked Integrated Power $\\Sigma P$ (a.u., log scale)")
        fig.savefig(TR_DIR / f"fig_all_maps_4pairs_{tim}.png", bbox_inches="tight")
        fig.savefig(TR_DIR / f"fig_all_maps_4pairs_{tim}.pdf", bbox_inches="tight")
        plt.close(fig)

        # --- Figure B: P_m vs Model C*a1^2*P_{f2/2} vs Residual across all 4 (f1, f2) pairs ---
        fig, axes = plt.subplots(4, 3, figsize=(13.5, 14.5), dpi=180, sharex=True, sharey=True)
        fig.suptitle(f"2D $(a_1, a_2)$ Validation of $P_m = C\\,a_1^2 P_{{f_2/2}}$ Across All 4 $(f_1, f_2)$ Pairs ({tim}, pcolormesh)",
                     fontsize=13, fontweight="bold", y=0.995)

        for row, (f1, f2) in enumerate(pairs):
            m_d = data[(f1, f2, "m", tim)]
            f2_d = data[(f1, f2, "f2", tim)]
            common = sorted(set(m_d.keys()) & set(f2_d.keys()))
            px = np.array([k[0] for k in common], dtype=float)
            py = np.array([k[1] for k in common], dtype=float)
            pm = np.array([m_d[k] for k in common], dtype=float)
            pf2 = np.array([f2_d[k] for k in common], dtype=float)

            mask_c = (px >= 60) & (py >= 350) & (pf2 > 1e-4) & (pm > 2e-6)
            C_med = float(np.median(pm[mask_c] / (px[mask_c]**2 * pf2[mask_c])))
            p_mod = C_med * (px**2) * pf2

            Zm = gaussian_pcolormesh_interp(px, py, np.maximum(pm, 1e-6), xg, yg)
            Zmod = gaussian_pcolormesh_interp(px, py, np.maximum(p_mod, 1e-6), xg, yg)
            with np.errstate(divide="ignore", invalid="ignore"):
                Zres = np.log10(Zm / Zmod)
            Zres[Zm < 4e-6] = np.nan

            im0 = axes[row, 0].pcolormesh(xg, yg, Zm, shading="gouraud", norm=norm_pm, cmap="viridis")
            im1 = axes[row, 1].pcolormesh(xg, yg, Zmod, shading="gouraud", norm=norm_pm, cmap="viridis")
            im2 = axes[row, 2].pcolormesh(xg, yg, Zres, shading="gouraud", vmin=-0.5, vmax=0.5, cmap="coolwarm")

            axes[row, 0].set_ylabel(f"({f1}, {f2}) MHz\n$C={C_med:.2e}$\n$a_2$ (mV)", fontsize=9)
            if row == 0:
                axes[0, 0].set_title("Measured Idler $P_m$", fontweight="bold")
                axes[0, 1].set_title("Model $C\\,a_1^2 P_{f_2/2}$", fontweight="bold")
                axes[0, 2].set_title("Log-Ratio $\\log_{10}(P_m / \\mathrm{Model})$", fontweight="bold")
            if row == 3:
                for c in range(3):
                    axes[3, c].set_xlabel("$a_1$ (mV)")

        cb0 = fig.colorbar(im1, ax=axes[:, :2], shrink=0.5, pad=0.02)
        cb0.set_label("Idler Power $P_m$ (log scale)")
        cb2 = fig.colorbar(im2, ax=axes[:, 2], shrink=0.5, pad=0.03)
        cb2.set_label("$\\log_{10}(P_m / \\mathrm{Model})$ (decades)")

        fig.savefig(TR_DIR / f"fig_pcolormesh_pm_vs_model_all4_{tim}.png", bbox_inches="tight")
        fig.savefig(TR_DIR / f"fig_pcolormesh_pm_vs_model_all4_{tim}.pdf", bbox_inches="tight")
        plt.close(fig)

    print("Saved D and ND maps to", TR_DIR)


if __name__ == "__main__":
    main()
