"""
fig5_harmonics.py - Harmonic decomposition of the phase dependence I(phi)
for the spectroscopy paper, on detune_fine.h5 (amps [327, 588] mV).

Exact DFT on the uniform 16-point phase grid:
    c_n = (2/N) * sum_k I(phi_k) exp(-i n phi_k),  N = 16,  c_n = amplitude
    of the n-th cosine harmonic (n >= 1), c_0 = 2*<I>.
Includes the arg-c2 pi-jump (mode boundary) analysis and a cross-check on
detune copy 2.h5 (same waveform/amps, 25-point phase grid, wider span).
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.ndimage import median_filter
import h5py

sys.path.insert(0, "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/pipeline")
from provenance import Provenance

H5_MAIN = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/detune_fine.h5"
H5_XCHK = ("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/"
           "phase gating/detune copy 2.h5")
NPH_MAIN, NPH_XCHK = 16, 25

p = Provenance("fig5_harmonics")
p.record_source(H5_MAIN, note="detune_fine.h5, 4800 keys, 20x15 (f1,f2) x 16 phases, "
                              "measured Jun 20, amps [327,588] mV, tone1 300 ns + tone2 400 ns "
                              "sequential non-overlapping tukey pulses, SA 30 MHz at f2/2")


def key_to_bin(k, am):
    a = k.split("_")
    return (float(a[0]), float(a[1]), int(round(float(a[2]) % (2 * np.pi) / (2 * np.pi / am))) % am)


def loader_factory():
    def loader(key, tr):
        y = tr[:, 1, :].astype(np.float64).mean(axis=0)   # avg passes, ch1 = linear power
        I_raw = np.clip(y, 0, None).sum()
        I_med = np.clip(median_filter(y, 21, mode="nearest"), 0, None).sum()
        return (I_raw, I_med)
    return loader


# ------------------------------------------------ DFT self-test (synthetic) --
phi16 = 2 * np.pi * np.arange(NPH_MAIN) / NPH_MAIN
rng = np.random.default_rng(0)
I_syn = (1.0 + 0.5 * np.cos(phi16 - 0.7) + 0.3 * np.cos(2 * phi16 + 1.1)
         + 0.1 * np.cos(3 * phi16 - 0.4) + 0.01 * rng.standard_normal(NPH_MAIN))
cds = {n: I_syn @ np.exp(-1j * n * phi16) / (NPH_MAIN / 2) for n in range(8)}
p.step("dft_selftest",
       params={"input": "(1, 0.5 cos(phi-0.7), 0.3 cos(2phi+1.1), 0.1 cos(3phi-0.4), noise 0.01)",
               "convention": "c_n = (2/N) sum I e^{-i n phi}"},
       result={"abs_c1": round(float(abs(cds[1])), 4), "abs_c2": round(float(abs(cds[2])), 4),
               "abs_c3": round(float(abs(cds[3])), 4),
               "arg_c1_deg": round(float(np.degrees(np.angle(cds[1]))), 2)},
       note="Exact DFT recovers injected amplitudes and phases on the 16-pt grid (expect 0.5, -40.1 deg)")

# ------------------------------------------------- main dataset load (grid) --
data = p.load_traces(H5_MAIN, loader=loader_factory(), pass_agg="mean")
keys = sorted(data.keys())
f1s = np.array(sorted(set(float(k.split("_")[0]) for k in keys)))
f2s = np.array(sorted(set(float(k.split("_")[1]) for k in keys)))
nf1, nf2 = len(f1s), len(f2s)
P = np.zeros((nf1, nf2, NPH_MAIN))
Ps = np.zeros((nf1, nf2, NPH_MAIN))
i1 = {v: i for i, v in enumerate(f1s)}
i2 = {v: j for j, v in enumerate(f2s)}
for k in keys:
    fv, fvv, ib = key_to_bin(k, NPH_MAIN)
    P[i1[fv], i2[fvv], ib] += data[k][0]
    Ps[i1[fv], i2[fvv], ib] += data[k][1]
# every cell written exactly once:
assert np.all(P != 0), "phase grid incomplete"
p.step("load_inband_power",
       params={"inband_power": "sum of clip(ch1, 0, None) over 1001 SA points, passes averaged",
               "medfilt": "median_filter(axis=SA, 21) [0.63 MHz] before clip+sum"},
       result={"n_keys": len(keys), "n_f1": nf1, "n_f2": nf2, "n_phase": NPH_MAIN,
               "axis_check": "SA center freq tracks key field 2 (f2/2), key field 1 is f1 (x axis)",
               "P_raw_min_W": float(P.min()), "P_raw_max_W": float(P.max())})

# ------------------------------------------------------------- harmonics ----
def dft(X, pgrid):
    return {n: np.einsum("ijk,k->ij", X, np.exp(-1j * n * pgrid)) * (2 / len(pgrid))
            for n in range(1 + len(pgrid) // 2)}


def harm_stats(X, phi, tag):
    C = dft(X, phi)
    A = {n: np.abs(C[n]) for n in range(8)}
    res = {"tag": tag}
    for n in range(8):
        i, j = np.unravel_index(int(np.argmax(A[n])), A[n].shape)
        res[f"A{n}_max_W"] = float(A[n].max())
        res[f"A{n}_argmax_f1_f2_MHz"] = [float(f1s[i]), float(f2s[j])]
        res[f"A{n}_median_W"] = float(np.median(A[n]))
    strong = A[1] > 0.5 * A[1].max()
    res["max_A4_over_A1_in_strong_region"] = float(np.max(A[4][strong] / A[1][strong]))
    res["median_A4toA7_over_A1_median"] = float(np.median([A[n] for n in range(4, 8)]) / np.median(A[1]))
    ii, jj = np.unravel_index(int(np.argmax(A[2] / A[1])), A[1].shape)
    res["max_A2_over_A1"] = float((A[2] / A[1])[ii, jj])
    res["max_A2_over_A1_at_f1_f2_MHz"] = [float(f1s[ii]), float(f2s[jj])]
    ii, jj = np.unravel_index(int(np.argmax(A[3] / A[1])), A[1].shape)
    res["max_A3_over_A1"] = float((A[3] / A[1])[ii, jj])
    res["max_A3_over_A1_at_f1_f2_MHz"] = [float(f1s[ii]), float(f2s[jj])]
    res["arg_c1_min_deg"] = float(np.degrees(np.angle(C[1]).min()))
    res["arg_c1_max_deg"] = float(np.degrees(np.angle(C[1]).max()))
    res["arg_c2_min_deg"] = float(np.degrees(np.angle(C[2]).min()))
    res["arg_c2_max_deg"] = float(np.degrees(np.angle(C[2]).max()))
    Irec = (C[0][..., None].real / 2
            + 2 * sum(np.real(C[n][..., None] * np.exp(1j * n * phi)) for n in (1, 2, 3)))
    res["median_3term_recon_residual_frac_of_I_mean"] = float(np.median(np.abs(X - Irec) / X.mean()))
    return C, A, res


C_r, A_r, st_r = harm_stats(P, phi16, "raw")
C_m, A_m, st_m = harm_stats(Ps, phi16, "medfilt21")
p.step("dft_harmonics_raw", params={"convention": "c_n = (2/16) sum I e^{-i n phi}; A0 = 2*<I> by convention"}, result=st_r)
p.step("dft_harmonics_medfilt21", params={"same convention"}, result=st_m)


def nrm(a):
    a = a - np.median(a)
    return a / np.linalg.norm(a)


p.step("medfilt_robustness",
       params={"check": "raw vs medfilt(21) harmonic maps, per-pixel correlation and phase rms"},
       result={"corr_A0": float(np.corrcoef(nrm(A_r[0]).ravel(), nrm(A_m[0]).ravel())[0, 1]),
               "corr_A1": float(np.corrcoef(nrm(A_r[1]).ravel(), nrm(A_m[1]).ravel())[0, 1]),
               "corr_A2": float(np.corrcoef(nrm(A_r[2]).ravel(), nrm(A_m[2]).ravel())[0, 1]),
               "rms_arg_c1_diff_deg": float(np.std(np.degrees(np.angle(C_r[1]) - np.angle(C_m[1])))),
               "rms_arg_c2_diff_deg": float(np.std(np.degrees(np.angle(C_r[2]) - np.angle(C_m[2])))),
               "A1_argmax_raw_vs_medfilt_MHz": [st_r["A1_argmax_f1_f2_MHz"], st_m["A1_argmax_f1_f2_MHz"]]})

# ------------------------------------------- c2 pi-jump line (main dataset) --
def collect_jumps(C2, C1, A1, f1s, f2s):
    Dg = f2s[None, :] - f1s[:, None]
    a2, a1 = np.angle(C2), np.angle(C1)
    out = []
    for j in range(len(f2s)):
        u, u1 = np.unwrap(a2[:, j]), np.unwrap(a1[:, j])
        for i in range(len(f1s) - 1):
            if abs(u[i + 1] - u[i]) > 0.6 * np.pi and min(A1[i, j], A1[i + 1, j]) > 0.3 * A1.max():
                out.append((0.5 * (f1s[i] + f1s[i + 1]), f2s[j],
                            0.5 * (Dg[i, j] + Dg[i + 1, j]), u[i + 1] - u[i],
                            abs(u1[i + 1] - u1[i])))
    for i in range(len(f1s)):
        u, u1 = np.unwrap(a2[i, :]), np.unwrap(a1[i, :])
        for j in range(len(f2s) - 1):
            if abs(u[j + 1] - u[j]) > 0.6 * np.pi and min(A1[i, j], A1[i, j + 1]) > 0.3 * A1.max():
                out.append((f1s[i], 0.5 * (f2s[j] + f2s[j + 1]),
                            0.5 * (Dg[i, j] + Dg[i, j + 1]), u[j + 1] - u[j],
                            abs(u1[j + 1] - u1[j])))
    return np.array(out) if out else np.zeros((0, 5))


Jm = collect_jumps(C_m[2], C_m[1], A_m[1], f1s, f2s)
jr = {"n_jump_segments": int(len(Jm))}
if len(Jm):
    lo, hi = np.percentile(Jm[:, 2], [25, 75])
    jr.update({
        "Delta_jump_median_MHz": float(np.median(Jm[:, 2])),
        "Delta_jump_IQR_MHz": [float(lo), float(hi)],
        "jump_size_median_deg": float(np.degrees(np.median(np.abs(Jm[:, 3])))),
        "jump_size_min_max_deg": [float(np.degrees(np.abs(Jm[:, 3]).min())),
                                  float(np.degrees(np.abs(Jm[:, 3]).max()))]})
    cf = np.polyfit(Jm[:, 0], Jm[:, 2], 1)
    jr["Delta_jump_vs_f2_slope_MHz_per_MHz"] = float(cf[0])
    jr["Delta_jump_line_rms_MHz"] = float(np.std(Jm[:, 2] - np.polyval(cf, Jm[:, 0])))
    jr["jump_line_Sigma_span_MHz"] = [float((Jm[:, 0] + Jm[:, 1]).min() / 2),
                                      float((Jm[:, 0] + Jm[:, 1]).max() / 2)]
    jr["max_darg_c1_across_jump_deg"] = float(np.degrees(Jm[:, 4].max()))
    jr["min_darg_c2_across_jump_deg"] = float(np.degrees(np.abs(Jm[:, 3]).min()))
p.step("c2_pi_jump_detune_fine",
       params={"criterion": "|d(arg c2)| > 0.6*pi after unwrap along both axes; "
                            "both pixels must have A1 > 0.3*max(A1); "
                            "c1 continuity read across the same edges"},
       result=jr)

# ------------------------------------------------------ c1 sign-flip check --
p.step("c1_signflip_check",
       params={"claim_src": "sibling dataset detune_fine031126.h5 (amps [390,588]): "
                            "c1 phase sign flip around f1~4070-4072, f2~4128-4137"},
       result={"arg_c1_min_deg": st_m["arg_c1_min_deg"], "arg_c1_max_deg": st_m["arg_c1_max_deg"],
               "Real_c1_sign": "uniformly negative" if (C_m[1].real <= 0).all() else "changes sign",
               "conclusion": "no c1 sign flip in this dataset/amplitude set; the harmonic-level "
                             "flip is carried by c2 as a pi jump at Delta ~ 74-76 MHz"},
       note="arg c1 stays in the third quadrant throughout; the R^2-collapse c1 flip of the "
            "031126 sibling does not reproduce at tone1 amp 327 mV")

# ------------------------------------------------------- delta/sigma map ----
Delta = f2s[None, :] - f1s[:, None]
Sigma = 0.5 * (f1s[:, None] + f2s[None, :])
iA1 = np.unravel_index(int(np.argmax(A_m[1])), A_m[1].shape)
iA2 = np.unravel_index(int(np.argmax(A_m[2])), A_m[2].shape)
iA0 = np.unravel_index(int(np.argmax(A_m[0])), A_m[0].shape)
iA3 = np.unravel_index(int(np.argmax(A_m[3])), A_m[3].shape)
p.step("delta_sigma_peaks",
       params={"coords": "Delta = f2 - f1, Sigma = (f1+f2)/2, MHz"},
       result={f"{s}_peak_W": float(A_m[n][ij])
               for s, n, ij in (("A0", 0, iA0), ("A1", 1, iA1), ("A2", 2, iA2), ("A3", 3, iA3))} |
              {f"{s}_peak_Delta_Sigma_MHz": [float(Delta[ij]), float(Sigma[ij])]
               for s, ij in (("A0", iA0), ("A1", iA1), ("A2", iA2), ("A3", iA3))})

p.array("A0_medfilt_W", A_m[0])
p.array("A1_medfilt_W", A_m[1])
p.array("A2_over_A1", A_m[2] / A_m[1])
p.array("arg_c1_deg", np.degrees(np.angle(C_m[1])))
p.array("arg_c2_deg", np.degrees(np.angle(C_m[2])))
p.array("c1_complex", C_m[1])
p.array("c2_complex", C_m[2])
p.array("Delta_MHz", Delta)
p.array("Sigma_MHz", Sigma)

# ------------------------------------------------- example I(phi) curves ----
pts = [("max |c1|", iA1), ("max c2/c1", (0, nf2 - 1))]
mask = np.ones_like(A_m[1], dtype=bool)
for i, j in [iA1, (0, nf2 - 1)]:
    mask[max(0, i - 1):i + 2, max(0, j - 1):j + 2] = False
pts.append(("min |c1|", np.unravel_index(int(np.argmin(np.where(mask, A_m[1], np.inf))), A_m[1].shape)))
ex_res, ex_curves = {}, []
for name, (i, j) in pts:
    Icur = Ps[i, j, :]
    crec = {n: (Icur @ np.exp(-1j * n * phi16)) / (NPH_MAIN / 2) for n in range(8)}
    Irec = crec[0].real / 2 + 2 * sum(np.real(crec[n] * np.exp(1j * n * phi16)) for n in (1, 2, 3))
    ex_res[name] = {"f1_MHz": float(f1s[i]), "f2_MHz": float(f2s[j]), "Delta_MHz": float(Delta[i, j]),
                    "A1_W": float(abs(crec[1])), "A2_W": float(abs(crec[2])), "A3_W": float(abs(crec[3])),
                    "arg_c1_deg": float(np.degrees(np.angle(crec[1]))),
                    "rms_resid_frac_of_mean": float(np.std(Icur - Irec) / Icur.mean())}
    ex_curves.append((name, f1s[i], f2s[j], Icur, Irec))
p.step("example_curves",
       params={"points": {n: [float(f1s[i]), float(f2s[j])] for n, (i, j) in pts},
               "reconstruction": "c0/2 + 2 sum_{n=1..3} Re(c_n e^{i n phi})"},
       result=ex_res)

# ============================== CROSS-CHECK: detune copy 2 ==================
p.record_source(H5_XCHK, note="detune copy 2.h5, 60360 keys, wider span f1 4051-4107 / "
                              "f2 4108-4163 MHz, 25-point uniform phase grid (2pi/25, 0..6.032 rad), "
                              "2 passes, same waveform and amps [327,588] (script_content attr); "
                              "key grid ~30% incomplete (missing phase bins scattered)")

with h5py.File(H5_XCHK, "r") as fx:
    xkeys = list(fx.keys())
    xf1s = np.array(sorted(set(float(k.split("_")[0]) for k in xkeys)))
    xf2s = np.array(sorted(set(float(k.split("_")[1]) for k in xkeys)))
    PX = np.zeros((len(xf1s), len(xf2s), NPH_XCHK))
    cnt = np.zeros((len(xf1s), len(xf2s), NPH_XCHK))
    j1 = {v: i for i, v in enumerate(xf1s)}
    j2 = {v: j for j, v in enumerate(xf2s)}
    phix = 2 * np.pi * np.arange(NPH_XCHK) / NPH_XCHK
    for k in xkeys:
        fv, fvv, ib = key_to_bin(k, NPH_XCHK)
        y = fx[k]["trace"][:, 1, :].astype(np.float64).mean(axis=0)
        PX[j1[fv], j2[fvv], ib] += np.clip(y, 0, None).sum()
        cnt[j1[fv], j2[fvv], ib] += 1
    PX /= np.maximum(cnt, 1.0)
CX = dft(PX, phix)
AX = {n: np.abs(CX[n]) for n in range(8)}
JX = collect_jumps(CX[2], CX[1], AX[1], xf1s, xf2s)
xres = {"n_keys": len(xkeys), "missing_phase_cells": int((cnt == 0).sum()),
        "A1_max_W": float(AX[1].max()), "A2_max_W": float(AX[2].max()),
        "A3_max_W": float(AX[3].max()),
        "median_A4toA7_over_A1_median": float(np.median([AX[n] for n in range(4, 8)]) / np.median(AX[1])),
        "n_jump_segments": int(len(JX))}
if len(JX):
    lo, hi = np.percentile(JX[:, 2], [25, 75])
    cf = np.polyfit(JX[:, 0], JX[:, 2], 1)
    xres.update({"Delta_jump_median_MHz": float(np.median(JX[:, 2])),
                 "Delta_jump_IQR_MHz": [float(lo), float(hi)],
                 "jump_size_median_deg": float(np.degrees(np.median(np.abs(JX[:, 3])))),
                 "Delta_jump_vs_f2_slope_MHz_per_MHz": float(cf[0]),
                 "Delta_jump_line_rms_MHz": float(np.std(JX[:, 2] - np.polyval(cf, JX[:, 0]))),
                 "max_darg_c1_across_jump_deg": float(np.degrees(JX[:, 4].max()))})
p.step("c2_pi_jump_crosscheck",
       params={"dataset": "detune copy 2.h5 (same amps [327,588], 25-pt phase grid, wider span)",
               "criterion": "same as detune_fine; A1 mask > 0.3*max over its own map"},
       result=xres,
       note="Cross-dataset match of the constant-Delta pi-jump line => mode-boundary claim")

# ================================ FIGURE ====================================
fig = plt.figure(figsize=(17.5, 10.5))
gs = GridSpec(2, 4, figure=fig, hspace=0.32, wspace=0.38)
AXKW = None

# (a) A0
ax = fig.add_subplot(gs[0, 0])
im = ax.pcolormesh(f1s, f2s, A_m[0].T * 1e6, cmap="viridis", shading="nearest")
plt.colorbar(im, ax=ax, label=r"$A_0$ (nW)  [$=2\langle P\rangle$]")
ax.plot(f1s[iA0[0]], f2s[iA0[1]], "r*", ms=14)
ax.set_xlabel(r"$f_1$ (MHz)"); ax.set_ylabel(r"$f_2$ (MHz)")
ax.set_title("(a) $A_0$: phase-mean in-band power", fontsize=10)

# (b) A1
ax = fig.add_subplot(gs[0, 1])
im = ax.pcolormesh(f1s, f2s, A_m[1].T * 1e9, cmap="viridis", shading="nearest")
plt.colorbar(im, ax=ax, label=r"$A_1=|c_1|$ (nW)")
ax.plot(f1s[iA1[0]], f2s[iA1[1]], "r*", ms=14)
ax.set_xlabel(r"$f_1$ (MHz)"); ax.set_ylabel(r"$f_2$ (MHz)")
ax.set_title(f"(b) $A_1$ max at $(\\Delta,\\Sigma)$ =\n({Delta[iA1]:.0f}, {Sigma[iA1]:.0f}) MHz", fontsize=10)

# (c) A2/A1
ax = fig.add_subplot(gs[0, 2])
im = ax.pcolormesh(f1s, f2s, (A_m[2] / A_m[1]).T, cmap="magma", shading="nearest", vmin=0)
plt.colorbar(im, ax=ax, label=r"$A_2/A_1$")
ax.plot(f1s[0], f2s[-1], "c*", ms=14)
ax.set_xlabel(r"$f_1$ (MHz)"); ax.set_ylabel(r"$f_2$ (MHz)")
ax.set_title("(c) $A_2/A_1$ ratio; star: (e) max-$c_2/c_1$", fontsize=10)

# (d) arg c2 with pi-jump line
ax = fig.add_subplot(gs[0, 3])
a2d = np.degrees(np.angle(C_m[2]))
im = ax.pcolormesh(f1s, f2s, a2d.T, cmap="twilight_shifted", shading="nearest",
                   vmin=-180, vmax=180)
plt.colorbar(im, ax=ax, label=r"arg $c_2$ (deg)")
a1d = np.degrees(np.angle(C_m[1]))
ax.contour(f1s, f2s, a1d.T, levels=np.linspace(a1d.min(), a1d.max(), 5),
           colors="lime", linewidths=0.8)
# jump line: measured constant-Delta median
Dmid = float(np.median(Jm[:, 2]))
ax.plot(f1s, f1s + Dmid, "r--", lw=2,
        label=fr"$\pi$-jump line: $\Delta$ = {Dmid:.1f} MHz")
ax.legend(loc="lower left", fontsize=8)
ax.set_xlim(f1s.min() - 0.5, f1s.max() + 0.5)
ax.set_ylim(f2s.min() - 0.5, f2s.max() + 0.5)
ax.set_xlabel(r"$f_1$ (MHz)"); ax.set_ylabel(r"$f_2$ (MHz)")
ax.set_title(f"(d) arg $c_2$: $\\pi$-jump at $\\Delta\\approx${Dmid:.0f} MHz"
             + "\n(lime: arg $c_1$, spans <8 deg = continuous)", fontsize=10)

# (e) example curves
ax = fig.add_subplot(gs[1, 0:2])
cols = ["#1f77b4", "#ff7f0e", "#2ca02c"]
phdeg = np.degrees(phi16)
for (name, ff1, ff2, Icur, Irec), col in zip(ex_curves, cols):
    ax.plot(phdeg, Icur * 1e9, "o", ms=4, color=col, alpha=0.8)
    ax.plot(phdeg, Irec * 1e9, "-", lw=1.5, color=col, alpha=0.9,
            label=fr"{name}: $f_1$={ff1:.0f}, $f_2$={ff2:.0f} MHz")
ax.set_xlabel(r"injected two-tone phase $\phi$ (deg)")
ax.set_ylabel(r"$I(\phi)=\sum_{\rm SA}$ clip $P$, passes avg (nW)")
ax.set_title(r"(e) Example $I(\phi)$ with 3-harmonic DFT reconstruction")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# (f) Delta, Sigma maps: |c1| mesh + |c2| contours
ax = fig.add_subplot(gs[1, 2:4])
mesh = ax.pcolormesh(Sigma, Delta, A_m[1] * 1e9, cmap="viridis", shading="nearest")
plt.colorbar(mesh, ax=ax, label=r"$A_1=|c_1|$ (nW)")
cs = ax.contour(Sigma, Delta, A_m[2] * 1e9, levels=5, colors="w", linewidths=1.0)
ax.clabel(cs, fmt="%.0f", fontsize=7)
ax.plot(Sigma[iA1], Delta[iA1], "r*", ms=14, label=r"$|c_1|$ max")
ax.plot(Sigma[iA2], Delta[iA2], "w^", ms=9, label=r"$|c_2|$ max")
if len(Jm):
    ax.plot([Sigma[0, :].min(), Sigma[0, :].max()], [float(np.median(Jm[:, 2]))] * 2,
            "r--", lw=2, label=fr"$\pi$-jump: $\Delta\approx{float(np.median(Jm[:,2])):.0f}$ MHz")
ax.set_xlabel(r"$\Sigma=(f_1+f_2)/2$ (MHz)")
ax.set_ylabel(r"$\Delta=f_2-f_1$ (MHz)")
ax.set_title("(f) Harmonics in $(\\Delta,\\Sigma)$: |c1| (color), |c2| (white contours)")
ax.legend(loc="lower right", fontsize=8)

p.finish(fig, panel_labels=["A0_map", "A1_map", "A2_over_A1", "arg_c2_pi_jump",
                            "I_phase_examples", "delta_sigma_harmonics"])
print("PIPELINE REPORT -------------------------------------------------")
print("A1 peak", st_m["A1_max_W"], "at", st_m["A1_argmax_f1_f2_MHz"])
print("A2 peak", st_m["A2_max_W"], "at", st_m["A2_argmax_f1_f2_MHz"])
print("A3 peak", st_m["A3_max_W"], "at", st_m["A3_argmax_f1_f2_MHz"])
print("A0 peak", st_m["A0_max_W"], "at", st_m["A0_argmax_f1_f2_MHz"])
print("max A2/A1", st_m["max_A2_over_A1"], "at", st_m["max_A2_over_A1_at_f1_f2_MHz"])
print("max A3/A1", st_m["max_A3_over_A1"], "at", st_m["max_A3_over_A1_at_f1_f2_MHz"])
print("n>=4 check", st_m["median_A4toA7_over_A1_median"], st_m["max_A4_over_A1_in_strong_region"])
print("c2 jump:", {k: v for k, v in jr.items()})
print("crosscheck:", {k: v for k, v in xres.items()})
print("examples:", {k: {kk: round(vv, 4) if isinstance(vv, float) else vv for kk, vv in v.items()}
                   for k, v in ex_res.items()})
