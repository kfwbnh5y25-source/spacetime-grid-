# Core messages — Paper A (two-tone phase gating / parametric idler spectroscopy)

Status key: VERIFIED (number extracted by us from data), CLAIMED (user's number, verification queued), THEORY (model prediction, not yet fitted).

## M1 — The f/2 idler is a coherent, phase-sensitive parametric response, not heating or leakage.
- Evidence: idler line at f/2 tracks the drive; sa_mask zeroing aside, gain up to +20.5 dB over the f2-only reference (f1_4069_5, a2=271 mV, f2=4111 MHz) — far above any pickup.
- I(phi) = A0 + A1 cos(phi_eff + phi0) with period 2pi (not cos 2phi): 4-phase Fourier c1 dominant in the pump2-scan map (|c1|/|c2| = 1.39), SVD rank-1 ~ 95%.
- Status: VERIFIED. Needs: absolute idler power calibration for the paper (SA ref level chain).

## M2 — Phase gating is strong and switch-like in the overlap (CW) regime.
- Evidence: +20.5 dB max gain at a2 = 271 mV/f2 = 4111 MHz; contrast between phi = 4.55 and phi = 1.31 states visible but single-point contrast in this a2 <= 275 mV dataset is modest (max |C| ~ 2.4 dB full-band) — the strong contrast lives at a2 approaching threshold (~ 288+ mV class and the A_old_params a2 200-588 data).
- ptp phase dependence in A_old_params (a2 up to 588 mV): up to ~45-47 units at f2 ~ 4135/430 mV; CSI2 dynamic range ~ 36-38 dB.
- Status: PARTIALLY VERIFIED. Needs: the money number (M4) pinned at one operating point.

## M3 — The interaction is NOT reciprocal: which pump you scan matters, because the two pumps play different roles.
- Evidence: fixed-pump1 4069/324 (scan pump2): two hot bands, cos(phi) gating, ±36-38 dB range. Fixed-pump 4130/266 (scan pump1): cos(2phi) dominant (c2/c1 median 1.9), saturation at ~347 mV, no suppression lobe.
- Physics reading: with pump1 strong and near resonance, pump2 couples into the pump1-driven mode; the interaction channel (four-wave mixing term ~ a1^2 a2* + c.c. vs a2^2 a1*) is asymmetric when amplitudes/roles differ — consistent with a driven-mode + coupler picture rather than a symmetric pair.
- Status: OBSERVED. Needs: tie to Gamma_eff model so Fig 2 caption can state WHY (this is the deeper-than-threshold result; candidate highlight message).

## M4 — Money number: phase gating depth (CORRECTED by boss).
- BOSS CORRECTION: the "50 percent" target is about PHASE GATING DEPTH — relative phase should
  modulate/idle the idler response by ~50 percent or MORE between the favorable and unfavorable
  phase ("the more the better"). It is NOT a threshold-reduction claim (the earlier
  threshold-a2-shift reading was a misinterpretation of the boss's prompt).
- Verification path: phase-contrast across datasets. RESULT (delay recheck): with depth
  defined as $1 - P_{\min}/P_{\max}$, overlap ($t_{\rm gap}=0$) data give 49–60% for
  $\Delta = 50$–$70$ MHz (delay_delta.h5, clean blocks); with the symmetrized depth
  $(P_{\max}-P_{\min})/(P_{\max}+P_{\min})$: 32–43%. fig3 overlay fit (two-phase, antipodal):
  gain ratio $A_1/A_0 = 0.011$ at $\Delta\phi = 185.7^\circ$ — the two-phase design cannot
  reach the cosine extremum by construction; 16-phase remeasure at
  $\phi_0/\phi_0+\pi$/quadratures recommended. Recorded phases 1.31/4.55 rad remain LOWER
  BOUNDS ($\phi_0 \approx -155^\circ$).
- Status: VERIFIED at the 49–60% level (delay_delta, one definition) / 32–43% (symmetrized);
  the >=50% claim depends on the depth definition — state the definition in the paper.

## M5 — Phase control equals time-of-flight control: phase accumulates in the delay between pulses.
- phi_eff = phi + 2 pi * Delta * t_gap (Delta = f2 - f1). Verified by the delay_vs_delta fit on delay_delta.h5 (fit convention matches). This converts a ns-scale delay knob into a well-defined phase knob, which is the gating mechanism of Fig 4.
- Status: VERIFIED (fit). Needs: one clean re-measure for Fig 4 quality.

## M6 — In the non-overlap fine-detune map, the regime boundary is a HARMONIC structure change, not a phase flip (CORRECTED).
- Older plots suggested R^2 collapse + c1 phase sign flip in the hot band; the fresh detune_fine031126.h5 analysis (truncated file:
  212/875 pairs, single pass) does NOT show a c1 sign flip: phi0 = -155 deg +- 4 deg stable, R^2 >= 0.94 everywhere.
- What changes at the A1 island instead: the cos-2phi harmonic COLLAPSES (A2/A1: 0.33 off-band -> 0.01-0.07 in the hot-band core);
  max A1 contrast 2.7 dB (hot-band core vs off-band median) at f1=4070, f2=4135, Delta = 65 MHz, A1 = 127.5 nW.
- BOSS PREDICTION CONFIRMED (fig5_harmonics, detune_fine.h5 complete grid + cross-check): **arg c2 shows a PI JUMP (~pi, median 132 deg)
  along a straight constant-Delta line at Delta = 74.5 MHz** (IQR 73.5-74.5, line rms 0.22 MHz in detune_fine.h5; Delta = 74.5 MHz,
  IQR 73.5-75.5, rms 0.72 MHz in detune copy 2.h5 — agrees across datasets) while **arg c1 varies <= 1.4 deg across the same line**.
  Interpretation: mode-pair identity switch at constant detuning = the paper's mode-mapping claim. STRONG (two datasets).
- Full-harmonic content (uniform 16-phase DFT, no fitting): c1 DOMINATES the strong-response zone (Delta ~ 65-85 MHz, Sigma ~ 4106-4116 MHz);
  A2/A1 <= 0.3 and A3/A1 <= 0.25 on the strong ridge; A2 competes only in the low-A1 high-Delta corner (Delta >= 85, where A1 is noise-limited,
  A2/A1 up to 0.67). Higher harmonics negligible: median |c4..c7| = 6.4 percent of median A1; max A4/A1 = 0.20 in the strong region.
  3-harmonic truncation validated; I(phi) reconstruction residual <= 6 percent rms except at the single max-|c1| point (36 percent,
  real higher-order content there).
- Peak values (medfilt maps, f1/f2 MHz): A0 = 635 nW at (4065,4154); A1 = 117 nW at (4075,4150); A2 = 35.6 nW at (4077,4142);
  A3 = 27.6 nW at (4080,4150).
- c1 sign flip: ABSENT in detune_fine.h5 (arg c1 confined to -156..-148 deg, Re c1 uniformly negative). The 031126 flip region sits at
  f2 4128-4137, outside this dataset's f2 window (>= 4140); the flip is carried by c2, not c1.
- Status: VERIFIED (M6 core: c2 pi-jump mode boundary, cross-dataset; harmonic dominance map; c1 continuous). Remaining: tie jump line
  to eigenmode spectrum (decisive test listed below); explain the max-|c1| point's higher-order content.

## M7 — Sequential vs overlapping regimes separate the physics: overlap = quasi-static four-wave interaction; non-overlap = time-of-flight phase imprint.
- Fig 3 vs Fig 4 is a controlled comparison (same tones, a1 327 -> 390 mV only + waveform change).
- STATUS UPDATE (delay recheck, independent DFT extraction from raw 16-phase sweeps):
  - fig4's tau = 2.4 +- 0.7 ns CONFIRMED for delay_f1.h5 conditions (f2 = 4142, hot band
    f1 4066-4078, i.e. Delta = 64-76): my per-f1 fits give hot-band median tau = 1.55 ns
    [IQR 1.08-1.84] (fit from 2 ns) / 1.80 ns [1.28-2.09] (fit from 1.5 ns = fig4's window).
    A1 plateau ~95-100 nW over t_gap 0-2 ns, collapsed to <=10 nW by ~5 ns; phi0 slew
    74-79% of 2 pi Delta in both datasets.
  - NEW: tau is DETUNING-DEPENDENT (delay_delta.h5, Sigma = 4106, Delta 50-88, t_gap <= 7 ns):
    at Delta 50-58 the ~30-36% contrast PERSISTS over the whole 7 ns window (no decay);
    1.1-2.5 ns near the c2 pi-jump line (Delta 74-78); < 1 ns at Delta >= 80.
  - So "the memory time" must be quoted as tau_mem(Delta); at Delta ~ 74-76 it is 2.4 ns
    (fig4 value), far below the naive 1/(2 Gamma_k) = 51-61 ns, but at Delta <= 58 the data
    give a LOWER BOUND tau > 7 ns COMPATIBLE with 51-61 ns. Candidate readings:
    nonlinear ringdown of the over-threshold driven seed (likely near the pi-line) and/or
    mode-pair hybridization changing the damping across the identity-switch line.
  - DATA HYGIENE: delay_f1.h5 and delay_delta.h5 also contain 4,745 "1EX" sweeps
    (SA span 60 MHz, center 4142 MHz, sweep 10 ms) — readoff field must be filtered to "1"
    in every analysis; ~9% of phase blocks carry a baseline-recalibration step
    (fake |c1|). All numbers above are from clean blocks.
- Status: FRAME + MEASURED (tau(Delta) curve now in GAMMA_EFF_MATH.md section 4.2).

## Math backbone (Gamma_eff doc) — now anchored on TWO AI-drafted reports + our data (user caution: reports unverified, treat as hypothesis generators)
Sources (read, cross-checked to each other):
- spin-waves/reports/idler_math/IDLER_MATH_REPORT.md (sub-threshold idler, self-energy, K_2B)
- spin-waves/reports/phase-gating/REPORT.md (NMM 4ABCD dissection, h_eff(Delta phi) formula, far-detuned turn-on kick)

Verified-in-data so far:
- Adiabatic idler law P_m ∝ a1^2 P_f2/2 (log-log 1.83-1.94, >2 decades); C_m ≈ K_2B within 10-20 percent = self-consistency headline.
- phi_eff = phi + 2 pi Delta t_gap fit (delay_delta) = EQUIVALENT to the phase-gating report's "relative phase at tone-2 turn-on" —
  delay and phase control are one knob. State this equivalence prominently in the paper.
- Candidate unified threshold law (from phase-gating report's h_eff(Delta phi) + idler K_2B):
      h2,th(phi) = h2,th(0) * [1 - K_2B h1^2 - K_phi h1 h4B cos(phi_eff)]  (form to be re-derived, not trusted)
  4B internal cross-pump = whole-Arnold-tongue vertical collapse (0.4212 -> 0.1093 mT claim, -74 pct; VERIFY analytically before printing);
  4A Kerr = pure horizontal shift, tip unmoved.

OPEN TENSION (candidate highlight result): phase-gating report predicts 4B averages OUT for Delta f >= 50 MHz in steady-state CW,
yet our CW-overlap data (f1_4069_5, Delta ~ 58-64 MHz) shows clean rank-1 cos(phi) with strong contrast. Hypothesis: burst-wise
(1,15,1) re-selection of attractor basins + SA integration over bursts = phase-dependent probability-weighted average. Test with
fig3/fig4 agent numbers; if confirmed, the paper contributes the resolution itself (steady-state vs burst-wise regime map).
Bistability band (R^2 collapse + phase sign flip) = separatrix signature, consistent with attractor-selection picture.

Trust rules before anything enters the paper: (1) re-derive h_eff(Delta phi) and the 4A/4B dissection closed-form ourselves;
(2) NMM tensor magnitudes (simulation-specific) = quote only as theory support;
(3) dB/x-fold claims in the reports are UNVERIFIED; (4) check the d'Aquino 2026 arXiv citation exists (2601.07798).
(5) BOSS RULE (standing): NO MODE NUMBERS anywhere in the paper, figures or captions until the eigenmode/bias-field fit is done —
    our simulated eigenmodes are on the boss's PC, the exact bias field is unknown, and the fit is slow. Refer to modes generically
    ("a parametric pair", "the driven mode", "mode-pair identity switch") only.

