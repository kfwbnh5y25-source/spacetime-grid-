# PAPER DRAFT — Phase-gated parametric spin-wave memory in a two-tone pumped magnetic microdisk

**Status: SKELETON.** Structure, figure specs, and all already-verified numbers are in place;
`[PENDING]` marks numbers that tonight's campaign (meas1–meas3) or later work fills in.
Math backbone: `GAMMA_EFF_MATH.md` (equations referenced below as GEM-x, NOT copied —
the .md files link, the paper text states its own derivations compactly).
`IDLER_MATH_REPORT.md` is hypothesis-generator only until verified (todo: read/verify).

**BOSS RULE (standing): NO mode numbers.** Modes are $k$ (seed) and $m$ (target/idler), pumps
"tone 1" ($f_1 = 2f_k$) and "tone 2" ($f_2 = 2f_m$), detuning $\Delta = f_2 - f_1$.

---

## Title options

1. "Phase-gated two-tone parametric excitation: relative-phase control of idler ignition,
   threshold, and memory in a BiYIG microdisk" 
2. "A 200 ns phase-coherent memory window between two parametric pumps" (if meas1 goes well)

Decision: after meas1. Draft uses option 1.

## Abstract (skeleton)

- [Setup 1 sentence] Two-tone parallel pumping of a confined magnets Modes pair ($k$ seeded by
  tone 1's half-frequency response, $m$ idler under tone 2).
- [Finding 1] The idler response is PHASE-GATED: controlled by the relative phase
  $\phi$ between tone 2 and the internal 4-wave-mixing cross-pump born in mode $k$
  (GEM Eq. 7: threshold modulation $\pm W/V$, depth grows toward threshold).
- [Finding 2 = money] The gate survives a time gap between the pumps: phase information
  persists $\tau_{\rm mem}(\Delta)$, measured from $>7$ ns at $\Delta \le 58$ MHz out to
  **[PENDING meas1] 0 → 200 ns**, testing the $e^{-2\Gamma_k t_{\rm gap}}$ law (GEM Eq. 10).
- [Finding 3] Threshold shift $h_{2,\rm th}(\Psi)$ follows
  $[\sqrt{\Gamma^2+\delta^2 - W^2\sin^2\Psi} - W\cos\Psi]/V$ **[PENDING meas2]**.
- [Claim] Burst-wise re-selection extends phase sensitivity beyond the CW-coherence bound
  $|\Delta| \lesssim \Gamma$ (GEM §4.3) — the memory is read out, not stored, by each burst.

## 1. Introduction

- Parametric excitation in confined geometries; two-tone/non-commutative steady states
  (cite d'Aquino 2026 arXiv:2601.07798v2; Soares 2026 arXiv:2601.05775 — VERIFIED citations).
- Gap: mode-interaction theory predicts phase-sensitive thresholds (Eqs. 55–60, Routh-Hurwitz
  App. C.3 of d'Aquino) but experiments report phase-averaged idler power only.
- Contributions: (i) phase-resolved threshold and gain at a fixed two-tone operating point;
  (ii) temporal sequence control (non-overlap gap) → phase memory measurement and its
  detuning dependence; (iii) closed-form threshold/gating model connecting all of it (GEM).

### 1b. Literature map (from boss's Zotero, lib 22169132, + web; verification level: V = full text verified, A = abstract-level, S = metadata incomplete — chase full ref before citing)

**Parametric pumping foundations.**
- Suhl 1957, J. Phys. Chem. Solids 1, 209 (A): first theory of high-power FMR
  instabilities; the Suhl processes are the parent of all parametric magnon work. Cite via
  Gurevich & Melkov; check the exact page before print.
- Gurevich & Melkov 1996, *Magnetization Oscillations and Waves* (book, V): standard
  monograph — parallel-pumping threshold, 3- and 4-magnon processes; our Eq. (2)
  convention-check source for laboratories.
- Stancil & Prabhakar 2009, *Spin Waves: Theory and Applications* (book, V): general
  textbook anchor; substrate/delay-line background.
- Brächer, Pirro, Hillebrands 2017, Phys. Rep. 699 (A): modern parallel-pumping review —
  amplification and manipulation of magnon spin currents; frames the field our operating
  point lives in. Full text skim for the intro paragraph still TODO.

**Mode-resolved spectroscopy and micromagnetics (our sample family).**
- d'Aquino, Serpico, Miano, Forestiere 2009, J. Comput. Phys. 228, 5134 (A/bib): numerical
  computation of magnetization modes; origin of the eigenmode-solver line we use.
- Srivastava et al. 2023, Phys. Rev. Applied 19, 064078 (A): parallel pumping triggers
  ~all fundamental eigenmodes of an in-plane-magnetized YIG microdisk; MRFM + micro-BLS
  label >40 modes. THIS is the sample/geometry lineage; the eigenmode-vs-bias fit (todo,
  boss PC) closes on their labeling.
- Perna et al. 2022, JMMM 546, 168683 (V, equation-level): NMM reduction — modal ODE
  hierarchy with higher-order nonlinear terms; our derivation backbone (GEM cites its
  Eqs. directly).
- Massouras et al. 2024, PRB 110, 064435 (A): mode-resolved micromagnetics of parallel
  pumping in disks — frequency pulling, mutual phase locking, higher-order scattering.
  Language support for the slew-shortfall (mode pulling) caveat T2.
- Abert 2019, EPJ B 92, 57 (A): micromagnetics models/numerics review — SM-only cite.
- Yemeli et al. 2025, arXiv:2503.07159 (A): high-power FMR of BiYIG nanodisks (OUR
  material system): saturation + self-modulation instability; mode count explodes above
  the critical threshold. Supports A1 weakening above threshold.

**Two-tone nonlinear interaction theory (direct prior).**
- d'Aquino, Perna, Merbouche, de Loubens 2026, arXiv:2601.07798 (V, page+equation
  verified): amplitude equations (42)-(43), steady states z/u/c±± (55)-(60), Routh-Hurwitz
  App. C.3; threshold (20)/(24). THEORY WE TEST.
- Soares et al. 2026, arXiv:2601.05775 (A): experimental two-tone parallel pumping in a
  YIG microdisk: steady state depends on mode pair, detuning, and TEMPORAL SEQUENCE of
  tones (non-commutativity). The companion experiment — we must differentiate (see below).

**Differentiation vs arXiv:2601.05775 (Perna/Soares 2026) — the "must cite and
differentiate" paragraph, to be placed in the Discussion:**
Both works study two-tone parallel pumping of a YIG microdisk and both observe
sequence-dependent (non-commutative) responses. The distinctions are sharp and go in
our favor on the temporal/phase axis: (i) 2601.05775 resolves STEADY STATES (which
attractor $z/u_h/u_n/c_{\pm\pm}$ the pair settles into, classified with the d'Aquino
2026 theory) and the ORDER of threshold crossings — i.e. memory of the excitation
SEQUENCE in a coarse, discrete sense. We resolve memory of the excitation PHASE: a
continuous relative phase $\phi$ and a variable gap $t_{\rm gap}$ between non-overlapping
pulses, with the response following the winding $\phi_{\rm eff} = \phi + 2\pi\Delta\,t_{\rm
gap}$ (GEM Eq. 9, VERIFIED) and a memory factor decaying as $e^{-2\Gamma_k t_{\rm gap}}$.
(ii) Their tones OVERLAP in steady state; ours are separated by a programmable gap —
the gap scan is the regime they do not access. (iii) We quantify a CONTINUOUS gating
observable ($P_m(\phi)$ contrast $\ge 50\%$ depth, threshold shift $h_{2,\rm th}(\Psi)$
vs phase) rather than a discrete attractor classification, which gives an extraction of
the internal 4B cross-pump amplitude $W$ in field units. (iv) They are free-standing
YIG microdisk, CW-detection regime; we operate pulsed with idler readout on the sample
path (BiYIG disks — the Yemeli 2025 material family), so their multi-attractor hopping
maps onto OUR far-detuned burst-wise regime (GAMMA §4.3) rather than onto the coherent
low-$\Delta$ regime where we measure the depth. Net: cite as the
sequence-non-commutativity companion; claim phase-continuum memory + gap-control as ours.

**Phase states of parametric oscillators (the ON/OFF physics).**
- Elyasi, Saitoh, Bauer 2022, PRB 105, 054403 (arXiv:2109.09117) (A): magnon parametron —
  bistable phase states above threshold form an Ising-like spin; predicts telegraph
  switching from thermal fluctuations (their "thermally activated switching" = our
  multi-attractor noise floor, GAMMA §4.3). NEW ADDITION — cite for phase-state framing.

**Coherence, mixing, combs (adjacent nonlinear magnon physics).**
- Hula et al. 2022, APL 121, 082401 (A): spin-wave frequency combs from cascaded
  four-magnon scattering — the comb mechanism our meas4 duty-control discriminates
  against (line scaling with duty × on-time).
- Dreyer et al. 2022, Nat. Comms 13, 4924 (A): imaging + phase-locking of nonlinear spin
  waves; coherence of half-integer processes. Supports §4.3 coherence question.
- An et al. 2024, Nat. Comms 15, 7238 (A): coherent nonlinear magnons in waveguides at
  ultrahigh frequency resolution — detection-methodology only.
- Li et al. 2026, Nat. Comms (DOI 10.1038/s41467-026-71916-9) (A): parametric pumping →
  spontaneous + idler propagating modes via four-wave mixing in a YIG delay line — the
  traveling-wave cousin of our 4B idler mechanism.
- Kim & Merbouche 2024, "Stimulated magnon scattering by non-degenerate parametric
  excitation", arXiv:2406.09032 (v1 Jun 2024, v2 Sep 2024); HAL hal-04625827 (A,
  resolved 1118): micromagnetics of two-frequency pumping generating NEW magnon modes
  via stimulated scattering; outcome depends on film geometry AND SEQUENCE of the pumping
  fields. Closest prior art to our tone-1 3-wave $k{+}m$ idler channel. STILL CHASE the
  published venue before print (arXiv may be the final version — it is a simulation-only
  companion, whereas our 3-wave channel is measured analytically extracted).
- Zautkin & Orel (with Cherepanov per Zotero note) ~1980, "Parametric excitation of spin
  waves by noise modulation of their frequencies", best match: Zh. Eksp. Teor. Fiz. 79,
  281–287 (1980) [Sov. Phys. JETP 52, 142–147] (A, resolved 1118 via citation lists in
  the L'vov S-theory school's collected reprints — NOT the article itself). S-theory-era
  precursor for noise-modulated / partially-coherent parametric excitation: supports our
  A3 (thermal seeding homogeneous) and the burst-wise statistics picture (GAMMA §4.3).
  CHASE the article (Google Scholar/library) to confirm title, author list
  (is Cherepanov on this one?), and page numbers before print.

**Neuromorphic framing (ONE sentence max in the paper — no hype).**
- Körber et al. 2023, Nat. Comms 14, 3723 (A): magnon-scattering reservoir — sequence-
  dependent mode scattering = the physical basis of our burst-wise re-selection picture.
- Papp, Porod, Csaba 2021, Nat. Comms 12, 6422 (A): spin-wave neural network hardware.
- Wright et al. 2022 Nature 601 & Momeni et al. 2025 Nature 646 (A): physical neural
  network training context; cite only if the "gating primitive" sentence survives.
- Brown 1963, Phys. Rev. 130, 1677 (V, in GAMMA doc): FDT thermal bath (SM).

## 2. Experimental setup

- Sample, antenna angles, bias field; SA chain; AWG 64 GS/s.
- Waveform grammar (Fig. 2 / Methods): `build_loop` pulse list $[P(a_i, f_i, t_0, T_{\rm on}, \phi)]$;
  June geometry for non-overlap: tone1 300 ns → gap $t_{\rm gap}$ → tone2 400 ns, 2 µs loop,
  32 ns tukey edges, duty 35%; CW-overlap geometry for Fig 3/meas3.
- **Data hygiene boxes**: exclude 1EX sweeps (different SA settings) — fig4/delay recheck
  learned this the hard way; baseline recalibration ~90 s vs ~8 s blocks → repeat-block
  design of meas3 kills the step artifacts.
- Safety envelope: freqs 3.5–4.5 GHz; |waveform| ≤ 1; peak ≤ [amplifier compression ceiling
  — value TBD from probe]; `MAX_PEAK_MV` scoping.

## 3. Model (compact derivation, full version in SM)

From the linearized NMM amplitude equations (d'Aquino Eqs. 42–43), with our $V$ convention
(off-diagonal pump element $iVh$, threshold $h_{\rm th} = \Gamma/V$, factor-2
versions elsewhere note their $V/2$):

1. Single-mode threshold: $h_{\rm th}^{(0)}(\delta) = \sqrt{\Gamma^2+\delta^2}/V$ (GEM 2);
   maps to d'Aquino Eq. (24) with $V \leftrightarrow \omega_h|(\varphi_h,\varphi_h^*)|/2$.
2. Internal 4B cross-pump $W e^{i\phi_W} \propto P_{mk} b_k^2$, stationary iff
   $|\Delta| \lesssim \Gamma$ (CW) or re-selected per burst (pulsed, GEM §4.3).
3. Threshold rule (THE formula):
   $$h_{2,\rm th}(\Psi) = \frac{\sqrt{\Gamma^2+\delta^2 - W^2\sin^2\Psi} - W\cos\Psi}{V}
     \xrightarrow{W\,\ll\,\Gamma} h_{\rm th}^{(0)} - \frac{W}{V}\cos\Psi \quad \text{(GEM 6, 7)}$$
4. Effective damping split (phase-insensitive $C h_1^2$ vs phase-sensitive $W\cos\Psi$,
   GEM 8); self-consistency $C_m \sim K_{2B}$ at 10–20% (idler report, VERIFIED scale).
5. Gap phase accumulation: $\Psi(t_{\rm on}) = \phi + 2\pi\Delta\, t_{\rm gap} + \rm const$
   (GEM 9, measured slew 74–80% of $2\pi\Delta$ — mode pulling, honest caveat).
6. Memory: $W(t_{\rm gap}) = W_0 e^{-2\Gamma_k t_{\rm gap}}$ (GEM 10) — low-$\Delta$ limit;
   measured collapse toward $\Delta \approx 74.5$ MHz (identity-switch line) = open.
   **[PENDING] meas1 at $\Delta = 50, 58$ MHz tests the limit case directly.**

## 4. Results — figure plan

### Fig 1 — SINGLE-TONE overview (FIXED per boss: Fig 1 is the one-tone figure)
Port of `figure1_single_tone.py` → `figures/fig1_overview.{png,pdf}` (DONE, per ledger).
Panels (a)–(i); (h) = pulse-timing diagram; (i) = RESERVED for sample characterization
(fallback eigen_27mT_*.npz) — HELD behind the BOSS RULE until the eigenmode-vs-bias-field
fit is done. Data: `Figure1_arnold_plot/data.h5` + `arnold_tongue_zoom/data2.h5`,
halfaverage obs, full drive calibration. Physics content: single-tone Arnold tongue,
threshold line $h_{\rm th}(\delta) = \sqrt{\Gamma^2+\delta^2}/V$ (GEM Eq. 2), zoom + hinge.
Prediction overlay (GEM Eq. 2 contours): synthesis step (todo: fig-synthesis).

### Fig 2 — Scaling & $\Gamma_{\rm eff}$ shift
- a: $P_m \propto a_1^2 P_{f_2/2}$ log-log, slopes 1.83–1.94 (VERIFIED).
- b: $a_1 \times f_1$ map A1 peak at $\Delta = 70$–80 MHz (july map) vs the 50–58 gap scan.
- c: threshold shift $C h_1^2$ (phase-insensitive term) — links idler K_{2B}.
Status: HAVE data, skeleton only.

### Fig 3 — Phase gating at fixed OP (depth claim) ← meas3 (+ old 16-phase overlap)
- a: repeat-averaged $P_m(\phi)$, 60 phases × ≥5 repeats, error bars; fresh baseline per block.
- b: cos fit, $A_1/A_0$, $\phi_0$ with CI.
- c: (optional panel) two-phase overlay fit from old data at $\phi_0 \approx -155°$
  (excluded by 031126) — or drop, superseded.
- **Money number**: $1 - P_{\min}/P_{\max} \ge 50\%$ with CI. Old single-shot estimate:
  49–60% (definition-sensitive; symmetric-depth version 32–43% — report BOTH definitions,
  primary = the symmetric amplitude ratio to be chosen by meas3 statistics). [PENDING]
- **NEW (meas5 edge sweep, 1117):** the pulse edge width is the phase-smearing knob —
  $A_1/A_0 = 0.38/0.35/0.19/0.045/0.073/0.042$ at edge 0.5/1/2/4/7/10 ns, end-control
  1 ns = 0.35 (no drift), A0 constant while A1 dies. Falloff scale ~2 ns matches
  $\tau_{\rm mem}(\Delta = 74) = \tau_b$. Methodological trap worth an SM paragraph:
  over-smoothed programming edges DESTROY the gating signal without moving the line
  power — all data must be taken with ~1 ns edges (June geometry). Candidate SM figure:
  A1/A0 vs edge width at fixed OP.

### Fig 4 — Memory / long-gap delay ← meas1 (the campaign centerpiece)
- a: $P_m(\phi)$ heat maps vs $t_{\rm gap}$ (0 → 200 ns), $\Delta = 50$ and $58$ MHz.
- b: extracted $A_1(t_{\rm gap})$ decay + exponential fit → $\tau_{\rm mem}(\Delta = 50, 58)$.
- c: $\phi_{\rm eff}$ slew check $d\phi_0/dt_{\rm gap}$ vs $2\pi\Delta$ (expect ~75–100%).
- Existing inset data: old delay.h5 stop at +15 ns; delay-recheck2: Δ=74 τ=2.13±0.08 ns,
  Δ=61 no decay to 7.4 ns. meas1 extends to 200 ns. [PENDING — the open-tension killer]

### Fig 5 — Threshold vs relative phase ← meas2 
- a: $h_{2,\rm th}$ vs $\Psi$ (2 mV steps, 4–8 phases over $2\pi$), fit to GEM (6) full form
  → extract $W/V$ independent of gain fits.
- b: regime map table as figure (GEM §5).
Status: script ready, data PENDING.

### Fig 6 (maybe, or 5b) — Identity switch / fine detune map
- $c_2$ $\pi$-jump at $\Delta = 74.5$ MHz (VERIFIED cross-dataset), $c_1$ continuous;
  July bracketing 76±2 / 103±2 at different $\Sigma$.
- Memory-collapsevs detuning: $\tau_{\rm mem}(\Delta)$ from old data + meas1 low-Δ points.

## 5. Discussion

- $\Gamma_{\rm eff}$ as instantaneous, programmable gate: device implication (isolation-free
  phase control of parametric onset; potential neuromorphic gating primitive — careful,
  one sentence max, no hype).
- Factor-2 / convention hygiene paragraph (defensive, short).
- Honest limitations: burst-wise re-selection interpretation is still a hypothesis
  (pulse-width tomography = proposed follow-up); $\delta_{\rm eff} \approx (0.26-0.37)\Delta$
  slew shortfall unexplained (mode pulling?); memory law fails near identity line.

## 6. Methods / SM

- Full derivation appendix from d'Aquino (42)–(43) → GEM 6/7/8 (todo: derive idler ODEs
  ourselves where IDLER_MATH_REPORT is wrong — pending verification).
- Numerical self-tests (tmp/test_gamma_eqs.py v2, ALL PASS: exact eigenvalue bisection etc.).
- Data pipeline: extraction scripts per fig, all (t_gap, f1) blocks listed, exclusions logged.

### SM-S1 — Sub-threshold self-energy, quadratic idler, and the measured coupling constants

Two-tone scattering picture (derived by us; sympy-checked, `tmp/check_idler2.py` +
`tmp/check_idler3.py`; NMM cross-check reproduced live, `idler_math/verify_idler_nmm.mjs`):

- Lower idler $m$ ($f_m = f_1 - f_2/2$): tone 1 scatters the tone-2-pumped seed into
  $P_m \propto a_1^2\,P_{f_2/2}$ (Fig 2a). Adiabatic elimination of $c_m^*$ gives the
  phase-INSENSITIVE self-energy $\Sigma_k = \kappa_{km}\kappa_{mk}^* h_1^2/(\Gamma_m + j\delta_m)$
  and the relative threshold shift $h_{2,\rm th}(h_1) = h_{2,\rm th}(0)(1 - K_{2B}h_1^2)$.
- Upper side $p$ ($f_p = f_2 - f_1/2$): EXACT MIRROR under $1\leftrightarrow2$,
  $k\leftrightarrow h$ (sympy-verified resonance condition); predicts
  $P_p \propto a_2^2\,P_{f_1/2}$ — VERIFIED in data ($t = +11.0$ to $+5.9$ across 4 pairs),
  weak only because $P_{f_1/2} \approx 50\times$ smaller than $P_{f_2/2}$ at the working
  point ($C_p \approx C_m$).
- Self-consistency: in the aligned limit ($\kappa_{km}=\kappa_{mk}$, $\beta_{mk}=1$,
  $\Gamma_k=\Gamma_m$) $K_{2B} \equiv C_m$ EXACTLY (sympy); measured 10–20% agreement
  (table below) = two independent measurements of the same tensor element.
- Absolute threshold reduction at the working point $a_1 = 327$ mV:
  $K_{2B} h_1^2 \approx 10\text{–}15\%$. NOTE: this is NOT the $\ge 50\%$ gating-depth
  claim (main text Fig 3); different observable, keep separate.

Measured constants (`threshold_reduction/data.h5`, 72,920 groups; $k$-mode linewidth
$\Gamma_k$, idler linewidth $\Gamma_m$; both constants in mV$^{-2}$):

| $(f_1, f_2)$ MHz | $\Gamma_k$ MHz | $\Gamma_m$ MHz | $C_m$ | $K_{2B}$ (ND) |
|---|---|---|---|---|
| 4071, 4125 | 1.51 | 1.54 | $1.16\times10^{-6}$ | $0.94\times10^{-6}$ ($R^2=0.92$) |
| 4071, 4135 | 1.61 | 1.51 | $8.54\times10^{-7}$ | $1.23\times10^{-6}$ |
| 4081, 4125 | 1.51 | 1.50 | $1.50\times10^{-6}$ | $1.36\times10^{-6}$ ($R^2=0.98$) |
| 4081, 4135 | 1.60 | 1.55 | $1.10\times10^{-6}$ | $1.25\times10^{-6}$ |

Both constants rise $+29\text{–}44\%$ when $f_1$: $4071 \to 4081$ MHz — same spectral
weighting of the coupling tensor; consistency check for the SM narrative.
(SM only: no mode numbers anywhere — BOSS RULE; the NMM-based absolute identification
stays out until the eigenmode-vs-bias-field fit.)

## Data provenance ledger (fill as figures land)

| Figure | Source file | Script | Status |
|---|---|---|---|
| Fig 1 | tmp/fig1_extract/*.npy | [to write] | skeleton |
| Fig 2 | fig2 extracts (prior session) | [locate] | skeleton |
| Fig 3 | meas3_depth_repeat.h5 | meas3_depth_repeat.py | tonight |
| Fig 4 | meas1_longgap_memory.h5 | meas1_longgap_memory.py | tonight |
| Fig 5 | meas2_threshold_vs_phase.h5 | meas2_threshold_vs_phase.py | tonight |
| Fig 6 | delay.h5 + delay2.h5 + detune | delay-recheck2 + july map | partial |
