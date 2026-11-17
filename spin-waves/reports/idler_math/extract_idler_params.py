"""
Extracts physical coupling, linewidth, and memory parameters from:
  - spin-waves/threshold_reduction/data.h5 (72,920 groups: f1 in {4071, 4081}, f2 in {4125, 4135})
  - spin-waves/neuromorphic/data_followup.h5 (E3FIX delay sweep 0..300 ns)

Outputs:
  - idler_extracted_params.json
  - fig_idler_theory_vs_data.png / .pdf
"""

import json
from pathlib import Path
import h5py
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SW_DIR = HERE.parent.parent
H5_TR = SW_DIR / "threshold_reduction" / "data.h5"
H5_FU = SW_DIR / "neuromorphic" / "data_followup.h5"


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


def lorentzian(df, amp, gamma, offset, f0):
    return offset + amp * (gamma**2) / ((df - f0)**2 + gamma**2)


def fit_linewidth(det_mhz: np.ndarray, mean_spec: np.ndarray):
    p0 = [float(mean_spec.max() - mean_spec.min()), 0.8, float(mean_spec.min()), 0.0]
    try:
        popt, pcov = curve_fit(lorentzian, det_mhz, mean_spec, p0=p0,
                               bounds=([0, 0.05, -1, -2], [1, 5.0, 1, 2]), maxfev=5000)
        gamma = float(popt[1])
        perr = float(np.sqrt(np.diag(pcov))[1])
        return gamma, perr, float(popt[3])
    except Exception:
        return np.nan, np.nan, 0.0


def estimate_threshold_hinge(a2_arr: np.ndarray, rs_arr: np.ndarray, thr: float = 2.0e-5, n_boot: int = 200, seed: int = 42):
    """
    Robustly estimate a2_th in an a1 bin by fitting a smooth boundary between
    sub-threshold (rs <= thr) and supra-threshold (rs > thr) points:
    finds the midpoint between highest sub-threshold and lowest supra-threshold
    (via logistic regression P(rs > thr | a2)) with bootstrap SE.
    """
    y = (rs_arr > thr).astype(int)
    if np.sum(y == 1) < 3 or np.sum(y == 0) < 3:
        return np.nan, np.nan
    # Fit 1D logistic decision boundary a2_th where P=0.5 using grid search over a2 for max accuracy / min cross-entropy
    def boundary(a2v, yv):
        # Grid search over candidate thresholds minimizing logistic loss / misclassification with tie-break by midpoint
        cands = np.linspace(float(a2v.min()), float(a2v.max()), 400)
        best_c, best_loss = np.nan, 1e18
        for c in cands:
            # soft hinge loss on transition width 4 mV
            logits = (a2v - c) / 3.0
            loss = np.sum(np.logaddexp(0.0, np.where(yv == 1, -logits, logits)))
            if loss < best_loss:
                best_loss = loss
                best_c = c
        return float(best_c)

    est = boundary(a2_arr, y)
    rng = np.random.default_rng(seed)
    boots = []
    n = len(a2_arr)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yb = y[idx]
        if np.sum(yb == 1) >= 2 and np.sum(yb == 0) >= 2:
            boots.append(boundary(a2_arr[idx], yb))
    se = float(np.std(boots)) if len(boots) > 5 else 2.0
    return est, max(se, 0.5)


def main():
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

    specs = np.array(specs, dtype=np.float64)
    clean = hampel_batch(specs)
    raw_sum = np.sum(clean, axis=1)
    det_mhz = np.linspace(-2.5, 2.5, 1001)

    sums = {}
    spec_buckets = {}
    p_buckets = {}
    for idx, (a1, a2, f1, f2, ro, tim) in enumerate(meta):
        sums.setdefault((f1, f2, ro, tim), {})[(a1, a2)] = raw_sum[idx]
        if ro == "m" and tim == "D" and a1 >= 220 and a2 >= 360:
            spec_buckets.setdefault((f1, f2, "m"), []).append(clean[idx])
        elif ro == "f1" and tim == "D" and a1 >= 310 and a2 <= 100:
            spec_buckets.setdefault((f1, f2, "f1"), []).append(clean[idx])
        elif ro == "f2" and tim == "D" and a1 <= 100 and a2 >= 360:
            spec_buckets.setdefault((f1, f2, "f2"), []).append(clean[idx])
        elif ro == "p" and a2 >= 350:
            if a1 >= 275:
                p_buckets.setdefault((f1, f2, "supra"), []).append(clean[idx])
            elif 50 <= a1 < 220:
                p_buckets.setdefault((f1, f2, "sub"), []).append(clean[idx])

    results = {"pairs": {}, "memory_e3fix": {}}

    # 8 fine bins across full [0, 355) range for ND (and D)
    a1_bins = [(0, 20), (20, 80), (80, 140), (140, 190), (190, 235), (235, 275), (275, 315), (315, 360)]

    for f1 in [4071, 4081]:
        for f2 in [4125, 4135]:
            lw = {}
            for ro in ["f1", "f2", "m"]:
                arr = np.array(spec_buckets.get((f1, f2, ro), []))
                if len(arr) > 0:
                    mu = np.median(arr, axis=0)
                    g_mhz, g_err, f0 = fit_linewidth(det_mhz, mu)
                    tau_ns = 1e3 / (2.0 * np.pi * g_mhz) if np.isfinite(g_mhz) else np.nan
                    lw[ro] = {"gamma_hwhm_mhz": g_mhz, "err_mhz": g_err, "f0_offset_mhz": f0, "tau_amp_ns": tau_ns}

            # Idler m scaling
            m_d = sums[(f1, f2, "m", "D")]
            f2_d = sums[(f1, f2, "f2", "D")]
            c_vals = [pm / (a1**2 * f2_d[(a1, a2)]) for (a1, a2), pm in m_d.items()
                      if 60 <= a1 <= 240 and a2 >= 350 and (a1, a2) in f2_d and f2_d[(a1, a2)] > 1e-4 and pm > 2e-6]
            c_med = float(np.median(c_vals))
            c_iqr = [float(np.percentile(c_vals, 25)), float(np.percentile(c_vals, 75))]

            # Weak upper side peak p (supra-th a1 vs sub-th a1)
            arr_p_sup = np.array(p_buckets[(f1, f2, "supra")])
            arr_p_sub = np.array(p_buckets[(f1, f2, "sub")])
            diff_sup = arr_p_sup[:, 400:600].mean(1) - 0.5 * (arr_p_sup[:, :200].mean(1) + arr_p_sup[:, 800:].mean(1))
            t_p_sup = float(diff_sup.mean() / (diff_sup.std(ddof=1) / np.sqrt(len(diff_sup)) + 1e-18))
            diff_sub = arr_p_sub[:, 400:600].mean(1) - 0.5 * (arr_p_sub[:, :200].mean(1) + arr_p_sub[:, 800:].mean(1))
            t_p_sub = float(diff_sub.mean() / (diff_sub.std(ddof=1) / np.sqrt(len(diff_sub)) + 1e-18))

            # Channel p scaling C_p = P_p / (a2^2 * P_{f1/2}) for a1 >= 275, a2 >= 350
            cp_vals = []
            for tim in ["D", "ND"]:
                p_t = sums[(f1, f2, "p", tim)]
                f1_t = sums[(f1, f2, "f1", tim)]
                for (a1, a2), pp in p_t.items():
                    if a1 >= 275 and a2 >= 350 and (a1, a2) in f1_t and f1_t[(a1, a2)] > 2e-6 and pp > 1e-6:
                        cp_vals.append(pp / (a2**2 * f1_t[(a1, a2)]))
            cp_med = float(np.median(cp_vals)) if cp_vals else np.nan

            # Threshold curves vs a1^2 for ND and D using logistic boundary + bootstrap SE
            th_curves = {}
            for tim in ["ND", "D"]:
                f2_t = sums[(f1, f2, "f2", tim)]
                bx, by, bse = [], [], []
                for lo, hi in a1_bins:
                    pts = [(a1, a2, rs) for (a1, a2), rs in f2_t.items() if lo <= a1 < hi]
                    if len(pts) < 15:
                        continue
                    a1_arr = np.array([p[0] for p in pts], dtype=float)
                    a2_arr = np.array([p[1] for p in pts], dtype=float)
                    rs_arr = np.array([p[2] for p in pts], dtype=float)
                    est, se = estimate_threshold_hinge(a2_arr, rs_arr)
                    if np.isfinite(est):
                        bx.append(float(np.mean(a1_arr) ** 2))
                        by.append(est)
                        bse.append(se)
                th_curves[tim] = {"a1sq": bx, "a2th": by, "se": bse}

            # Weighted linear regression on ND across FULL sub-threshold + ND curve
            bx_nd = np.array(th_curves["ND"]["a1sq"])
            by_nd = np.array(th_curves["ND"]["a2th"])
            w_nd = 1.0 / (np.array(th_curves["ND"]["se"]) ** 2)
            # Weighted polyfit degree 1
            coef = np.polyfit(bx_nd, by_nd, 1, w=np.sqrt(w_nd))
            slope_a2, a2_0 = float(coef[0]), float(coef[1])
            y_hat = a2_0 + slope_a2 * bx_nd
            ss_res = np.sum(w_nd * (by_nd - y_hat) ** 2)
            ss_tot = np.sum(w_nd * (by_nd - np.average(by_nd, weights=w_nd)) ** 2)
            r2_nd = float(1.0 - ss_res / max(ss_tot, 1e-12))
            k_2b_mv2 = float(-slope_a2 / a2_0)

            results["pairs"][f"{f1}_{f2}"] = {
                "linewidths": lw,
                "idler_C_med_mv2": c_med,
                "idler_C_iqr_mv2": c_iqr,
                "channel_p_t_supra": t_p_sup,
                "channel_p_t_sub": t_p_sub,
                "channel_p_C_med_mv2": cp_med,
                "a2_th_0_mv": a2_0,
                "a2_th_slope_mv_per_mv2": slope_a2,
                "K_2B_mv_minus2": k_2b_mv2,
                "nd_fit_r2": r2_nd,
                "th_curves": th_curves,
            }

    # E3FIX relaxation fit
    if H5_FU.exists():
        with h5py.File(H5_FU, "r") as f:
            for amp in [350, 380]:
                ds = [0, 10, 20, 40, 80, 150, 300]
                js = []
                for d in ds:
                    ab, ba = [], []
                    for sk in f.keys():
                        if not sk.endswith("__SIG"):
                            continue
                        nm = sk[:-5]
                        if nm.startswith(f"E3FIX_a{amp}_AB_d{d}_"):
                            bl = f[nm + "__BL"]["trace"][:, 1, :].mean(0)
                            ab.append(f[sk]["trace"][:, 1, :].mean(0) - bl)
                        elif nm.startswith(f"E3FIX_a{amp}_BA_d{d}_"):
                            bl = f[nm + "__BL"]["trace"][:, 1, :].mean(0)
                            ba.append(f[sk]["trace"][:, 1, :].mean(0) - bl)
                    A = hampel_batch(np.array(ab))
                    B = hampel_batch(np.array(ba))
                    j = float(np.sum((A.mean(0) - B.mean(0))**2) / (np.sum(A.var(0) + B.var(0)) + 1e-12))
                    js.append(j)
                popt, _ = curve_fit(lambda d, j0, tau, jinf: j0 * np.exp(-d / tau) + jinf,
                                    np.array(ds, float), np.array(js, float), p0=[js[0], 10.0, 0.06],
                                    bounds=([0, 1.0, 0], [50, 200.0, 1.0]))
                results["memory_e3fix"][str(amp)] = {"delays_ns": ds, "fisher_j": js, "tau_J_ns": float(popt[1])}

    (HERE / "idler_extracted_params.json").write_text(json.dumps(results, indent=2))

    # --- 4-Panel Publication Figure ---
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.5), dpi=170)
    fig.suptitle("Experimental Parameter Extraction vs. Two-Mode Parametric Self-Energy Theory", fontsize=13, fontweight="bold")

    # (a) Modal Linewidths (f1/2, m, f2/2)
    ax = axes[0, 0]
    for (f1, f2, ro, col, lbl) in [
        (4071, 4125, "f1", "#9467bd", "f1/2 (2035.5 MHz)"),
        (4081, 4125, "f1", "#8c564b", "f1/2 (2040.5 MHz)"),
        (4081, 4125, "m", "#d62728", "Idler m (2018.5 MHz)"),
        (4081, 4125, "f2", "#2ca02c", "f2/2 (2062.5 MHz)"),
    ]:
        mu = np.median(spec_buckets[(f1, f2, ro)], axis=0)
        mu_n = (mu - mu.min()) / (mu.max() - mu.min() + 1e-12)
        g_mhz = results["pairs"][f"{f1}_{f2}"]["linewidths"][ro]["gamma_hwhm_mhz"]
        ax.plot(det_mhz, mu_n, lw=1.3, alpha=0.85, color=col, label=f"{lbl}: HWHM={g_mhz:.2f} MHz")
    ax.set_xlabel("Detuning $f - f_c$ (MHz)")
    ax.set_ylabel("Normalized PSD")
    ax.set_title("(a) Modal Linewidth Extraction ($\\Gamma_1, \\Gamma_k, \\Gamma_m$)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (b) Proper Threshold Reduction Fit across FULL a1^2 range (ND vs D)
    ax = axes[0, 1]
    for f1, col in [(4071, "#1f77b4"), (4081, "#d62728")]:
        pr = results["pairs"][f"{f1}_4125"]
        tc_nd = pr["th_curves"]["ND"]
        tc_d = pr["th_curves"]["D"]
        bx_nd = np.array(tc_nd["a1sq"]) / 1e4
        ax.errorbar(bx_nd, tc_nd["a2th"], yerr=tc_nd["se"], fmt="s", ms=5, color=col, capsize=3,
                    label=f"f1={f1} ND (R²={pr['nd_fit_r2']:.2f}, $K_{{2B}}={pr['K_2B_mv_minus2']:.2e}$)")
        xg = np.linspace(0, bx_nd.max(), 50)
        ax.plot(xg, pr["a2_th_0_mv"] + pr["a2_th_slope_mv_per_mv2"] * xg * 1e4, "--", lw=1.6, color=col)
        bx_d = np.array(tc_d["a1sq"]) / 1e4
        ax.errorbar(bx_d, tc_d["a2th"], yerr=tc_d["se"], fmt="o-", ms=5, lw=1.2, alpha=0.75, color=col, mfc="white",
                    label=f"f1={f1} D (4B pre-exc. split)")
    ax.set_xlabel("Tone-1 Drive Power Proxy $a_1^2$ ($\\times 10^4\\,\\mathrm{mV}^2$)")
    ax.set_ylabel("Tone-2 Threshold $a_{2,\\mathrm{th}}$ (mV)")
    ax.set_title("(b) $f_2=4125\\,\\mathrm{MHz}$: Logistic Threshold Fit (`ND` Linear in $a_1^2$, `D` Split)")
    ax.legend(fontsize=7.5, loc="lower left")
    ax.grid(True, alpha=0.3)

    # (c) Weak Upper Side Peak in Channel p (supra-th a1 vs sub-th a1)
    ax = axes[1, 0]
    for (f1, f2, col, fc) in [(4071, 4125, "#1f77b4", 2089.5), (4081, 4125, "#d62728", 2084.5)]:
        arr_sup = np.array(p_buckets[(f1, f2, "supra")])
        arr_sub = np.array(p_buckets[(f1, f2, "sub")])
        # 50-pt moving boxcar smoothing to expose weak ~5e-9 peak cleanly
        k_sm = np.ones(40) / 40.0
        mu_sup = np.convolve(arr_sup.mean(0), k_sm, mode="same") * 1e9
        mu_sub = np.convolve(arr_sub.mean(0), k_sm, mode="same") * 1e9
        t_s = results["pairs"][f"{f1}_{f2}"]["channel_p_t_supra"]
        ax.plot(det_mhz[30:-30], mu_sup[30:-30], lw=1.8, color=col, label=f"({f1},{f2}) p@{fc} supra-th a1 (t={t_s:+.1f})")
        ax.plot(det_mhz[30:-30], mu_sub[30:-30], lw=1.2, ls="--", color=col, alpha=0.6, label=f"({f1},{f2}) sub-th a1 (chance)")
    ax.set_xlabel("Detuning $f - f_p$ (MHz)")
    ax.set_ylabel("Smoothed Excess PSD (pW)")
    ax.set_title("(c) Weak Upper Side Peak (Channel `p`, Turns On Only for $a_1 > a_{1,\\mathrm{th}}$)")
    ax.legend(fontsize=7.5)
    ax.grid(True, alpha=0.3)

    # (d) E3FIX Relaxation Curve
    ax = axes[1, 1]
    for amp, col in [("350", "#1f77b4"), ("380", "#2ca02c")]:
        if amp in results["memory_e3fix"]:
            me = results["memory_e3fix"][amp]
            ax.semilogy(me["delays_ns"], me["fisher_j"], "o-", lw=1.6, color=col,
                        label=f"{amp} mV ($\\tau_J={me['tau_J_ns']:.1f}$ ns)")
    ax.set_xlabel("Inter-Symbol Delay $d$ (ns)")
    ax.set_ylabel("Order Fisher Separability $J$")
    ax.set_title("(d) Relaxation Memory Horizon Extraction (E3FIX)")
    ax.legend(fontsize=8.5)
    ax.grid(True, which="both", alpha=0.3)

    plt.tight_layout()
    fig.savefig(HERE / "fig_idler_theory_vs_data.png", bbox_inches="tight")
    fig.savefig(HERE / "fig_idler_theory_vs_data.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Done! Saved idler_extracted_params.json and fig_idler_theory_vs_data.png")


if __name__ == "__main__":
    main()
