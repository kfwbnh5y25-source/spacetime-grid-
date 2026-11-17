# Synthesis figure for Paper A (GAMMA_EFF_MATH.md section 5, last row).
# Overlay of predicted threshold/memory structure, Eqs (2),(6),(7),(8),(8a),(9),(10)
# [W and Psi via (9)+(10)], with the MEASURED points/tag numbers from the ledger.
# BOSS RULE: no mode numbers anywhere. Modes: k (seed), m (target/idler).
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

OUT = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral"

# ---------------------------------------------------------------- scales
# Gamma tension (GAMMA_EFF_MATH.md section 1, flagged 1117):
G_NMM = (8.2e6, 9.8e6)          # MHz*1e6, NMM report scale for alpha=0.002 (UNVERIFIED sim)
G_MEAS = (1.50e6, 1.61e6)       # measured HWHM linewidths, threshold_reduction/data.h5
G_DYN = (185e6, 233e6)          # dynamic decay 1/(2 tau), delay + "3" datasets (VERIFIED fits)
TAU_B = 2.15                    # ns, verified build-up fit ("3" dataset)
TAU_B_ERR = 0.08

fig = plt.figure(figsize=(13.0, 8.2))
fig.suptitle("Synthesis: two-tone threshold and phase-memory predictions vs. measured points\n"
             r"modes $k$ (seed) and $m$ (target/idler); $V$ convention: pump element $= iVh$",
             fontsize=11)

# =====================================================================
# (a) Arnold threshold vs detuning: Eq (2) with the two Gamma scales
# =====================================================================
ax = fig.add_subplot(2, 3, 1)
d = np.linspace(0.0, 25.0, 400)                        # MHz
fill = ax.fill_between(d, np.sqrt(G_MEAS[0]**2 + (d*1e6)**2)/1e6,
                       np.sqrt(G_MEAS[1]**2 + (d*1e6)**2)/1e6,
                       color="#d62728", alpha=0.25, lw=0)
l1, = ax.plot(d, np.sqrt(1.55**2 + d**2), color="#d62728", lw=1.6,
              label=r"measured linewidths $\Gamma=1.5$--$1.6$ MHz (VERIFIED)")
fill2 = ax.fill_between(d, np.sqrt(G_NMM[0]**2 + (d*1e6)**2)/1e6,
                        np.sqrt(G_NMM[1]**2 + (d*1e6)**2)/1e6,
                        color="#1f77b4", alpha=0.25, lw=0)
ax.plot(d, np.sqrt(9.0**2 + d**2), color="#1f77b4", lw=1.6)
l2 = plt.Rectangle((0, 0), 1, 1, fc="#1f77b4", alpha=0.25)
ax.set_xlabel(r"pump detuning from mode half-frequency, $\delta = 2\pi(f_m - f_p/2)$ (MHz)")
ax.set_ylabel(r"$V\,h_{\rm th}(\delta) = \sqrt{\Gamma^2+\delta^2}$   (rate scale, MHz)")
ax.set_title(r"(a) single-tone threshold, Eq. (2): the $\Gamma$ scale tension")
ax.legend(handles=[l1, l2], loc="lower right", fontsize=8)
# measured phase-INSENSITIVE reduction, Eq (8) first term: 1 - K2B a1^2 = 10-15% at a1=327 mV
y0a = 1.55
ax.annotate("", xy=(2.6, y0a*0.86), xytext=(2.6, y0a),
            arrowprops=dict(arrowstyle="-|>", color="k", lw=1.4))
ax.plot([2.6], [y0a*0.875], "k.", ms=5)
ax.text(3.4, 1.62, "measured, Eq. (8):\n" r"$1-K_{2B}a_1^2 = -10$--$15\%$"
        "\n" r"at $a_1=327$ mV (VERIFIED,\n" "relative: scale-free)",
        fontsize=7.5)
ax.grid(alpha=0.3, lw=0.5)
ax.text(0.97, 0.05,
        "blue band = NMM sim. scale $\Gamma=8.2$--$9.8$ MHz (UNVERIFIED);\n"
        "detuning $\delta\\neq\\Delta$: $\delta_m \\sim$ few MHz\n"
        "while the tone-tone detuning $\\Delta = 50$--$110$ MHz",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7)

# =====================================================================
# (b) Phase-resolved threshold: Eq (6) exact + Eq (7) linear; inset Eq (8a) depth
# =====================================================================
axb = fig.add_subplot(2, 3, 2)
psi = np.linspace(0, 2*np.pi, 721)
for w, col, lbl in [(0.04, "#1f77b4", "NMM-scale $w=W/\\Gamma$ (sim., UNVERIFIED)"),
                    (0.005, "#d62728", r"measured bound $|\!h_{4B}/h_2\!|\approx 0.005$ (fig3, VERIFIED bound)")]:
    hexact = np.sqrt(1 - (w*np.sin(psi))**2) - w*np.cos(psi)
    hlin = 1 - w*np.cos(psi)
    axb.plot(np.rad2deg(psi), hexact, color=col, lw=1.6, label=lbl + r", Eq. (6)")
    axb.plot(np.rad2deg(psi), hlin, color=col, lw=1.0, ls="--", alpha=0.8,
             label=r"$\qquad$linearized, Eq. (7)" if w == 0.04 else None)
axb.set_xlim(0, 360)
axb.set_ylim(0.90, 1.10)
axb.set_xlabel(r"relative phase at turn-on, $\Psi = \phi - \phi_W$ (deg)")
axb.set_ylabel(r"$h_{\rm th}(\Psi)\,V/\sqrt{\Gamma^2+\delta^2}$")
axb.set_title(r"(b) phase-gated threshold, Eqs. (6)/(7) [$\delta=0$]")
axb.legend(fontsize=7, loc="upper right")
axb.grid(alpha=0.3, lw=0.5)
axb.text(180, 1.075, r"peak-to-valley $= 2W/V =: 2h_W$ (measures $W$ in field units)",
         ha="center", fontsize=8)
axb.text(0.02, 0.03,
         r"dashed $\equiv$ linearized (7) $\equiv \Gamma_{\rm eff}/\Gamma = 1 - w\cos\Psi$ [(8)];"
         "\n" r"period $2\pi$ (not $\pi$) - falsifiable",
         transform=axb.transAxes, fontsize=7)
# inset: sub-threshold depth, Eq (8a)
axins = axb.inset_axes([0.52, 0.10, 0.46, 0.36])
x = np.linspace(0.0, 0.999, 900)
for w, col in [(0.04, "#1f77b4"), (0.005, "#d62728")]:
    xm = x[x < 1 - w]
    depth = 1 - ((1 - xm - w)/(1 - xm + w))**2
    axins.plot(xm, depth, color=col, lw=1.2, label=None)
axins.axhspan(0.49, 0.60, color="green", alpha=0.18, lw=0)
axins.axhspan(0.32, 0.43, color="gray", alpha=0.15, lw=0)
axins.plot([0.970, 0.978], [0.49, 0.60], "r_", ms=8, mew=1.6)
axins.plot([0.758, 0.822], [0.49, 0.60], "b_", ms=8, mew=1.6)
axins.set_xlim(0.70, 1.0)
axins.set_ylim(0, 1)
axins.set_xlabel(r"$x=Vh_2/\Gamma$", fontsize=7, labelpad=1)
axins.set_ylabel(r"$1-P_{\min}/P_{\max}$", fontsize=7, labelpad=1)
axins.tick_params(labelsize=6)
axins.set_title(r"Eq. (8a) depth", fontsize=7, pad=2)
axins.text(0.98, 0.16, "measured\ndepth\n49-60%\n(green)", fontsize=5.5,
           ha="right", color="darkgreen",
           transform=axins.transAxes)

# =====================================================================
# (c) Memory law (10): tau_mem(Delta) - theory bands vs measured points
# =====================================================================
axc = fig.add_subplot(2, 3, 3)
axc.set_yscale("log")
axc.set_ylim(0.05, 700)
axc.set_xlim(45, 95)
axc.axhspan(51, 61, color="#1f77b4", alpha=0.25, lw=0, zorder=1)
axc.axhspan(313, 331, color="#d62728", alpha=0.18, lw=0, zorder=1)
axc.axhspan(2.1, 2.7, color="#ff7f0e", alpha=0.22, lw=0, zorder=1)
axc.axvspan(74, 88, color="crimson", alpha=0.10, lw=0, zorder=0)
axc.axvline(74.5, color="k", ls="--", lw=1.0, zorder=2)
axc.text(74.9, 320, "identity-switch line\n" r"$\Delta=74.5$ MHz (VERIFIED)", fontsize=7)
# measured points: filled = exp-fit tau; open = lower bound; diamond = half-contrast time
axc.errorbar([76], [1.8], yerr=[[0.7], [0.7]], fmt="o", color="k", ms=5,
             capsize=3, label=r"16-phase ranges (VERIFIED)", zorder=5)
axc.errorbar([66], [2.75], yerr=[[1.25], [1.25]], fmt="o", color="k", ms=5, capsize=3, zorder=5)
axc.errorbar([84], [0.45], yerr=[[0.35], [0.35]], fmt="o", color="k", ms=5, capsize=3, zorder=5)
axc.errorbar([74], [2.13], yerr=[[0.08], [0.08]], fmt="s", color="crimson", ms=6,
             capsize=3, zorder=6)
axc.plot([54], [7.0], marker="v", mfc="none", mec="k", ms=7, ls="none", zorder=5)
axc.annotate(r"$\tau>7$ ns at $\Delta\leq 58$", xy=(54, 7.2), xytext=(46.5, 18),
             fontsize=7.5, arrowprops=dict(arrowstyle="-|>", lw=0.8))
axc.plot([61], [7.4], "D", mfc="w", mec="k", ms=6, zorder=5)
axc.annotate(r"$t_{1/2}=7.4$ ns, no decay" "\n" r"to 6.5 ns ($\Delta=61$)", xy=(61.4, 7.4),
             xytext=(63.5, 60), fontsize=7.5,
             arrowprops=dict(arrowstyle="-|>", lw=0.8))
axc.annotate(r"$\Delta=74$: $\tau=2.13\pm0.08$ ns" "\n" "(not a clean exponential)",
             xy=(74, 2.13), xytext=(76.5, 11), fontsize=7.5,
             arrowprops=dict(arrowstyle="-|>", lw=0.8))
axc.set_xlabel(r"tone-tone detuning $\Delta = f_2 - f_1$ (MHz)")
axc.set_ylabel(r"phase-contrast memory time (ns)")
axc.set_title(r"(c) memory, Eq. (10): $\tau_{\rm mem}=1/(2\Gamma_k)$ vs. $W(t)=W_0e^{-t/\tau}$")
axc.grid(alpha=0.3, lw=0.5, which="both")
axc.text(45.8, 95, r"Eq. (10) with NMM $\Gamma$ (UNVERIFIED): 51--61 ns", fontsize=7.5,
         color="#1f4d80")
axc.text(45.8, 480, r"Eq. (10) with measured $\Gamma_k=1.5$--$1.6$ MHz: 310--330 ns",
         fontsize=7.5, color="#8b1a1a")
axc.text(45.8, 3.4, r"dynamic decay $1/(2\tau_b)=185$--233 MHz (VERIFIED fit)", fontsize=7.5,
         color="#a34d00")
axc.text(58, 0.085,
         r"PREDICTION (10) FAILS for $\Delta \gtrsim 74$ MHz"
         "\n(near the identity-switch line): measured collapse\n"
         "1-3 orders below both theory bands.\n"
         "Honest reading: $\tau_{\rm mem}(\Delta)$ is detuning-\n"
         "dependent; (10) is the LOW-$\\Delta$ limit only\n"
         "(nonlinear seed ring-down / hybrid-mode\n"
         "damping candidates).", fontsize=6.8, color="crimson")

# =====================================================================
# (d) Phase slew, Eq (9): dphi0/dt_gap = 2 pi Delta vs measured
# =====================================================================
axd = fig.add_subplot(2, 3, 4)
dd = np.linspace(40, 100, 200)
axd.plot(dd, 0.36*dd, "k-", lw=1.4, label=r"Eq. (9): $2\pi\Delta$")
axd.fill_between(dd, 0.74*0.36*dd, 0.79*0.36*dd, color="gray", alpha=0.25, lw=0,
                 label=r"measured shortfall band $74$--$79\%$")
axd.errorbar([61], [15.32], yerr=[0.04], fmt="o", color="crimson", ms=5, capsize=3,
             label=r"60-phase data (VERIFIED)")
axd.errorbar([74], [20.70], yerr=[0.07], fmt="o", color="crimson", ms=5, capsize=3)
axd.plot([70], [20.9], "^", color="crimson", ms=6, mfc="w", label=r"fig4, $\Delta=70$ MHz")
axd.text(59.5, 13.6, "69.7%\n(79% in\ncollapse-free\nwindow)", fontsize=7)
axd.text(74.6, 22.0, "77.7%", fontsize=7.5)
axd.text(71.0, 18.1, "82.9%", fontsize=7.5)
axd.set_xlabel(r"$\Delta$ (MHz)")
axd.set_ylabel(r"$d\varphi_0/dt_{\rm gap}$ (deg/ns)")
axd.set_title(r"(d) phase $\equiv$ delay, Eq. (9): $\varphi_{\rm eff} = \varphi + 2\pi\Delta\,t_{\rm gap}$")
axd.legend(fontsize=7.5, loc="upper left")
axd.grid(alpha=0.3, lw=0.5)
axd.text(0.97, 0.05,
         r"form of (9) VERIFIED; effective detuning reduced by"
         "\n" r"$\delta_{\rm eff}\approx(0.26$--$0.37)\Delta$ (mode-pulling candidate, OPEN)",
         transform=axd.transAxes, ha="right", fontsize=7)

# =====================================================================
# (e) Turn-on edge smearing (meas5): A1/A0 vs edge width
# =====================================================================
axe = fig.add_subplot(2, 3, 5)
te = np.array([0.5, 1, 2, 4, 7, 10])
r = np.array([0.38, 0.35, 0.19, 0.045, 0.073, 0.042])
axe.semilogx(te, r, "o-", color="#2ca02c", ms=6, label=r"meas5, $\Delta=74$ op. point (VERIFIED)")
axe.semilogx([1], [0.34], "s", mfc="none", mec="#2ca02c", ms=8,
             label=r"sanity GO-run, 1 ns edge, 0.34 (VERIFIED)")
axe.semilogx([32], [0.005], "v", mfc="none", mec="crimson", ms=7,
             label=r"32 ns tukey: $\leq 0.005$ (sanity, VERIFIED)")
axe.axvspan(2.07, 2.23, color="orange", alpha=0.25, lw=0)
axe.set_ylim(0, 0.46)
axe.set_xlim(0.4, 45)
axe.set_xlabel(r"tone-2 turn-on edge width $t_{\rm edge}$ (ns)")
axe.set_ylabel(r"$A_1/A_0$  (phase-gating amplitude)")
axe.set_title("(e) gate survives sharp edges only (meas5: Outcome A)")
axe.legend(fontsize=7, loc="upper right")
axe.grid(alpha=0.3, lw=0.5, which="both")
axe.text(0.02, 0.72,
         r"orange band: measured smearing scales $\tau_{\rm mem}(\Delta{=}74)=2.13$ ns"
         "\n" r"$\approx\tau_b=2.15$ ns (ring-up, VERIFIED)"
         "\n" "end-control 1 ns re-run: 0.35 (no drift)",
         fontsize=7)
axe.text(0.02, 0.30,
         "NO parameter-free theory curve here: (9)/(10)\n"
         "assume a delta-function $t_{\\rm on}$.\n"
         "The nonzero $\\approx5$--$7\\%$ floor for $t_{\\rm edge}\\geq4$ ns is OPEN.",
         fontsize=7)

# =====================================================================
# (f) Threshold build-up ("3" dataset): a_th(tau_p) vs Eq with tau_b
# =====================================================================
axf = fig.add_subplot(2, 3, 6)
tp, ath = zip(*[
    (0.3125, 590), (0.4688, 590), (0.6250, 580), (0.7812, 560), (0.9375, 530),
    (1.0938, 520), (1.2500, 480), (1.4062, 480), (1.5625, 460), (1.7188, 460),
    (1.8750, 440), (2.0312, 440), (2.1875, 420), (2.3438, 420), (2.500, 420),
    (2.656, 420), (2.812, 410), (2.969, 390), (3.125, 400), (3.281, 390),
    (3.438, 390), (3.594, 390), (3.750, 390), (3.906, 370), (4.062, 370),
    (4.219, 360), (4.375, 370), (4.531, 360), (4.688, 360), (4.844, 370),
    (5.000, 350), (5.156, 360), (5.312, 340), (5.469, 340), (5.625, 350),
    (5.781, 350), (5.938, 350), (6.094, 350), (6.250, 350), (6.406, 350),
    (6.562, 340), (6.719, 330), (6.875, 340), (7.031, 330), (7.188, 330),
    (7.344, 340), (7.500, 330), (7.656, 330), (7.812, 330), (7.969, 330),
    (8.125, 320), (8.281, 320), (8.438, 340), (8.594, 330), (8.750, 330),
    (8.906, 330), (9.062, 330), (9.219, 320)])
tp = np.array(tp); ath = np.array(ath)
tb = np.linspace(0.2, 10, 400)
axf.plot(tb, 320/np.sqrt(1 - np.exp(-tb/TAU_B)), color="#9467bd", lw=1.6,
         label=r"$a_\infty/\sqrt{1-e^{-\tau_p/\tau_b}}$, $\tau_b=2.15\pm0.08$ ns (FIT, form asserted)"
         )
axf.fill_between(tb, 320/np.sqrt(1 - np.exp(-tb/(TAU_B-TAU_B_ERR))),
                 320/np.sqrt(1 - np.exp(-tb/(TAU_B+TAU_B_ERR))),
                 color="#9467bd", alpha=0.22, lw=0)
axf.plot(tp, ath, ".", color="k", ms=5,
         label=r"single-tone threshold $a_{\rm th}(\tau_p)$ (VERIFIED, 50%-duty train)")
axf.axhline(320, color="gray", ls=":", lw=1.0)
axf.text(6.0, 328, r"$a_\infty=320$ mV: saturation for $\tau_p\gtrsim 8$ ns", fontsize=7.5)
axf.set_xlim(0, 10)
axf.set_ylim(250, 650)
axf.set_xlabel(r"pulse width $\tau_p$ (ns)")
axf.set_ylabel(r"single-tone threshold $a_{\rm th}$ (mV)")
axf.set_title('(f) population build-up ("3" dataset): rates out of the\n'
              r"small-signal $\Gamma$ window")
axf.legend(fontsize=7, loc="upper right")
axf.grid(alpha=0.3, lw=0.5)
axf.text(0.03, 0.06,
         r"This curve forces $\Gamma_{\rm dyn}=1/(2\tau_b)\approx233$ MHz:"
         "\n" r"$\tau_b=8$ ns saturation EXCLUDES the NMM 8--10 MHz"
         "\n" r"scale at these drives (panel (d),(c) tension)",
         transform=axf.transAxes, fontsize=7)

fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUT + "/fig_synthesis.png", dpi=200)
fig.savefig(OUT + "/fig_synthesis.pdf")
print("saved", OUT + "/fig_synthesis.png and .pdf")

# quick numeric footprints for the manifest
w = 0.005
for dd_ in (0.49, 0.60):
    print("w=0.005 depth %.2f -> 1-x = %.4f" % (dd_, w*(1+np.sqrt(1-dd_))/(1-np.sqrt(1-dd_))))
w = 0.04
for dd_ in (0.49, 0.60):
    print("w=0.04 depth %.2f -> 1-x = %.4f" % (dd_, w*(1+np.sqrt(1-dd_))/(1-np.sqrt(1-dd_))))
