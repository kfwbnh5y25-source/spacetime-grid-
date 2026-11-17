# Independent Audit, Verified Findings, and Revised Experimental Plan

**Scope:** `MANUSCRIPT_DRAFT.md`, `IDLER_MATH_REPORT.md`, `reports/phase-gating/REPORT.md`, and the raw HDF5 datasets + analysis scripts in `spin-waves/neuromorphic/`.
**Method:** Every claim below was checked against the actual repository code and data files, not accepted from the drafts.
**Status of drafts:** Potentially strong publication. The physics is real and well controlled. Four claims in the current draft are **overstated or computed with non-standard definitions** and must be corrected before submission. One central comparison (digital twin) was **never actually performed on measured data**.

---

## 0. Plain-Language Summary (for the supervisor)

We drive a chip full of small YIG magnetic disks with two microwave tones and watch what comes back at half the frequency on a spectrum analyzer. We showed that tone 1 "wakes up" the magnetic oscillations and tone 2 only ignites if tone 1 came first — the disk remembers. The court is still out on one question a referee will certainly ask: *could a dumb electronic circuit (a filter plus a diode-like nonlinearity) produce the same classification outcome?* We wrote a script that claims to answer this, but on inspection it never touches our measured data — it compares our hardware score against a toy simulation of its own making. That comparison must be redone properly (Section 1). Separately, three different memory timescales appear in the draft under one name (Section 2), the memory-capacity and information-processing numbers include generous conventions (Sections 3–4), and the amplifier-distortion worry has a clean answer in our own data that we should state explicitly (Section 5). The main novelty — phase-dependent gating rooted in 4-magnon cross-pumping, verified in simulation — now has its own quantitative theory (`reports/phase-gating/REPORT.md`) and should be promoted to the paper's headline (Section 7).

---

## 1. CRITICAL ISSUE 1 — Digital Twin: NOT YET TESTED (highest priority)

**Claim in draft:** "LTI + static-NL digital twin achieves 100% on AB/BA and 99.4% on S3, so we must show our reservoir does more."

**What the code actually does** (`spin-waves/neuromorphic/benchmark_digital_baselines.py`):

- The "twin" is never fitted to, nor driven by, any measured trial. It builds its **own synthetic surrogate waveforms** (`lti_static_twin_features()` generates sine waves with hard-coded proxy frequencies `freqs = [1.0, 2.3, 3.7]` plus a fixed `0.3*u^3` distortion, plus injected Gaussian noise).
- The digital ESN likewise runs on synthetic one-hot symbol steps, not measured inputs.
- The physical-reservoir numbers (100%, 92.2%) are **hard-coded literals** copied from measured metrics.

So the comparison is: {measured hardware result} vs. {toy self-generated benchmark}. The reported twin/ESN accuracies carry **no evidential weight** in either direction. Worse, presenting this table invites the referee to conclude the twin question was asked and answered affirmatively.

**What a valid twin test looks like** (geometry is favorable to us):

1. For every measured trial we know the exact injected drive waveform ${u}_{trial}(t)$ (the AWG sample array saved per trial). Pass each through: (a) a fitted complex LTI filter chain $H(f)$ (fit on a calibration subset), then (b) a static memoryless nonlinearity (fit static map from instantaneous input amplitude to output, e.g. cubic/clip).
2. Extract from the twin output time-integrated band-power features exactly analogous to what we take from the SA (integrate over the same 2030–2075 MHz window).
3. Train/test the same ridge readout with the same folds.
4. Compare against measured-spectra readout. The honest null is: *an LTI + static-NL system driven identically and read out identically.*

**Physics argument already in our favor (state it explicitly):** a cascade of linear filters and *memoryless* nonlinearities adjoint to the amplifiers can produce only frequencies of the form $m f_1 \pm n f_2$ ($m,n$ integer, $m+n$ odd for cubic). With $f_1=4081$, $f_2=4125$ MHz, **no integer combination lands in the 2030–2075 MHz readout band** ($2f_1-f_2=4037$, $2f_2-f_1=4169$, third harmonics $\approx 4.08{-}4.16$ GHz, etc.). Half-frequency ($f/2$) subharmonic emission is **topologically impossible** for that class. Therefore every photon in the readout band comes from genuine parametric magnon dynamics. The twin test then reduces to asking whether in-band *information* (not in-band *emission mechanism*) could be replicated by out-band surrogate features — which is exactly what the corrected protocol above answers.

**Verdict:** Claim **unproven and wrongly documented**. Fix before submission; cost: pure offline analysis, no rig time.

---

## 2. CRITICAL ISSUE 2 — Three Inconsistent Memory Timescales

**Confirmed.** The same $[A][d][B]$ Fisher-$J$ vs. delay dataset (`data_followup.h5`, E3FIX) has been fit **twice with different pipelines**, giving different numbers, and a third, physically distinct quantity shares the same label $\tau_{\text{mem}}$:

| Number in drafts | Source | Actual quantity |
| --- | --- | --- |
| $\tau_J = 4.3$ ns (350 mV), $9.4$ ns (380 mV) | `reports/idler_math/extract_idler_params.py:236-240` | Exp. fit to Fisher $J(d)$, Hampel-despiked, bounds $\tau\in[1,200]$ ns |
| $\tau_{\text{mem}} = 6.3$ ns (350 mV), $12.3$ ns (380 mV) | `neuromorphic/compute_sota_rc_benchmarks.py:155-175` | Same data, **no despiking**, different fit init/bounds $\tau\in[1,500]$ ns |
| "20–80 ns" | `MANUSCRIPT_DRAFT.md` line 52 | Threshold-reduction memory horizon from pre-excitation experiments — **a different physical quantity** (how long a pre-pulse stays useful), not a relaxation rate |

Additional problems: only 7 delay points, dominated by the enormous $J(0)$ (10.4 at 350 mV); no uncertainty on either fit; a 3-parameter fit to 7 noisy points with the baseline estimated from the last 3 (nearly constant) points is weakly constrained.

**Fix (one canonical pipeline):**

- Distinguish and relabel: (i) $\tau_{\text{sep}}$ = class-separability relaxation from E3FIX; (ii) $\Delta t_{\text{th}}$ = threshold-reduction horizon from the pre-excitation experiments; (iii) symbol-overlap memory from the equalization tasks.
- Single fitting script, Hampel despiking on, bootstrap over trials (30 reps per point) → report $\tau \pm$ CI.
- Purge the draft of "20–80 ns" in the E3FIX context and of the 6.3/12.3 vs. 4.3/9.4 duplication. Quote one number per quantity with error bars.

---

## 3. CRITICAL ISSUE 3 — Memory Capacity: MC = 1.170 includes lag 0

**Confirmed by code.** `compute_sota_rc_benchmarks.py` computes "MC" as the sum of the three degree-1 Walsh capacities over a 3-symbol BIT3 window: lag0 ($u_3$, the symbol nearest the acquisition) = 0.640, lag1 = 0.358, lag2 = 0.172. Lag 0 is **input transduction** (which symbol is present now, legible in the frequency-resolved spectrum as the $T_{f/2}$ vs $f_2/2$ line ratio), not memory of past input. Under the conventional definition ($\sum_{k\ge1}$), persistent memory is **MC = 0.530 symbols**.

Two further caveats to state in the paper:

1. This "MC" is computed on 3-bit word trials, not on a continuous i.i.d. random stream over many lags (Jaeger's MC definition). It is a finite-window proxy and must be labeled as such.
2. A frequency-resolved static readout trivially reads symbol *identity* at lag 0 (that's what the 0.640 measures) — which is exactly why the interesting claims are the lag-1/lag-2 terms and the quadratic $C_{\text{deg2}} = 0.412$, both genuinely memory-mediated.

**Fix:** report MC$_{k\ge1} = 0.530$ (or recompute properly on an i.i.d. stream if a new rig run is allowed) and rename lag-0 capacity "input transduction capacity." Keep the contrast with the 0.3x control (0.001) — that one is devastating and clean.

---

## 4. CRITICAL ISSUE 4 — IPC = 1.582: finite-window, leakage-prone

**Confirmed.** In `compute_walsh_ipc()`:

- IPC is computed over **7 Walsh terms of one 3-symbol window** — not the unrestricted Walsh matrix. Can't be compared to literature IPC values directly.
- Splitting is plain shuffled `KFold` over repeated measurements of the **same 8 input patterns** → train and test folds contain repeats of identical inputs. With 25 reps/word this inflates capacity estimates (the readout memorizes word-conditioned spectra, then recognizes them).
- No uncertainty, no null distribution.

**Fix:** (a) label it "3-symbol window IPC (finite)"; (b) grouping-aware evaluation: leave-reps-out blocks, or better leave-one-word-family-out for the quadratic/parity terms; (c) bootstrap CI per Walsh coefficient; (d) report the sub-threshold control (0.116) under the *same* folds — the contrast survives any reasonable correction and is the real result.

---

## 5. CRITICAL ISSUE 5 — RF Artifacts: answerable from our own data

**Confirmed numbers** (`e0_e4_metrics.json`, E0 EX block): amplifier IMD3 products at 4065/4155 MHz measured at $-64.2/-67.8$ dBm against tones at $-32.9/-32.8$ dBm → **$-31.3$ to $-34.9$ dBc**. The drafts' instinct to wave this off ("$\le -31$ dBc, negligible") is not a referee-proof argument on its own. The correct arguments, both already supported:

1. **Out-of-band:** both cubic products of 4081/4125 (and 4095/4125) land at $>4.0$ GHz, hundreds of MHz above the 2030–2075 MHz readout band; the measured IMD3 lines confirm where they land.
2. **In-band null measurements still missing:** measure the electronics-only floor in the readout band (device detuned far below the parallel-pumping band or removed; looped-sequence drive identical otherwise). This is a single ~10-min rig run and closes the issue permanently.
3. The off-resonance control (47.5% = chance) is already an operational version of this — frame it as such.

---

## 6. Statistics — genuine gaps confirmed

- **Random splitting of repeated measurements:** all CV evaluations (`analyze_multisymbol_v2.py`, `analyze_spoken_digits*.py`, `compute_sota_rc_benchmarks.py`) shuffle at trial level, so train/test share repeated measurements of the same patterns. With only 6–8 distinct patterns, the numbers measure reproducibility of class-mean separation, not generalization to novel sequences. **Add:** grouped/leave-rep-blocks-out accuracy for every headline number (S3 92.2%, BIT3 60%, etc.), and leave-one-speaker-out for spoken digits. The frozen cross-session E5 result (100%, zero recalibration) is the strongest generalization evidence and should be highlighted as the primary provenance claim.
- **Permutation p-values:** current reporting (p = 0.0049 with 200 permutations) is fine but underpowered at the tail; use the exact finite form $p=(b+1)/(N+1)$ and $N\ge10{,}000$, one-line change.

---

## 7. Paper Framing: Phase Gating + Idler Physics as the Main Novelty

The phase-gating story is the strongest and most differentiated contribution, and it now has a complete microscopic theory:

- **Theory:** `reports/phase-gating/REPORT.md` (NMM ODE model; BOSS RULE: no mode numbering anywhere until the eigenmode/bias-field fit is done — exact bias field unknown) shows the full Arnold-tongue downward collapse is caused **uniquely by the 4B coherent 4-magnon cross-pump term** ($\sum_{k\neq h} P_{hk}\, a_h^* a_k^2$): tip threshold $0.4212 \to 0.1093$ mT, with 4A contributing only a horizontal Kerr shift (tip unmoved) and 4C/4D subdominant. Smooth ($\text{SEM}<0.5\%$) phase response is predicted for quasi-degenerate detunings ($\Delta f \le 3.5$ MHz); jagged multi-seed attractor hopping for far detuning.
- **Experiment:** the AB/BA phase sweep in `data_multisymbol_v2.h5` (accuracy 83.3–100%, Fisher $J$ rising 0.20 → 0.65 as $\Delta\phi \to 180°$, 25% modulation of total in-band emission) is the hardware signature of the same coherent interference. **The measured phase dependence is a feature, not a bug**: the paper should *not* claim phase-robust computing; it should claim phase-gated computing whose Δφ-response follows the NMM prediction.
- **Bridge sentence (use verbatim):** *"The spectroscopy paper explains why the first pulse lowers the parametric threshold (4-magnon pre-excitation seeds the idler, $P_m \propto a_1^2 P_{f_2/2}$); the neuromorphic paper shows what that seeded state computes when a second pulse arrives inside the seeded coherence window — and the phase of the second pulse controls it, exactly as the 4B cross-pump term predicts."*

**Missing quantitative link:** a hardware measurement of the Δφ-response curve (transmitted/ignited intensity vs. relative phase at quasi-degenerate detuning) overlaid with the simulated curve (`fine_quasidegen_data.json`). This is the single highest-value rig experiment remaining (~30–60 min) and would convert the simulation-only Arnold-tongue result into an experimental claim.

---

## 8. Two-Paper Strategy — endorsed, with assignments

**Paper A — Nonlinear magnon spectroscopy (idler physics).** Complete. Deliverables: `IDLER_MATH_REPORT.md`, `threshold_reduction/fig_*_4pairs_{D,ND}`, $P_m \propto a_1^2 P_{f_2/2}$ slopes 1.83–1.94, threshold-reduction exponent, D/ND split.

**Paper B — Phase-gated physical reservoir computing.** Proof-of-compute results are strong (100% order + frozen cross-session, 92.2% S3, matched-count above static cap, MC$_{k\ge1}=0.530$, tuned nonlinear EQ 88.8%/97.5%). Before submission it needs, in priority order:

| # | Action | Type | Cost | Addresses |
| --- | --- | --- | --- | --- |
| 1 | Digital twin fed with **measured drive waveforms**, grouped CV | offline | 0 rig time | Issue 1 (fatal as written) |
| 2 | Unified E3FIX memory re-fit, one pipeline, bootstrap CI, relabel three timescales | offline | 0 | Issue 2 |
| 3 | MC/IPC re-computation: MC$_{k\ge1}$, window-labeled IPC, grouped folds, CIs | offline | 0 | Issues 3–4 |
| 4 | Finite permutation tests: $p=(b+1)/(N+1)$, $N=10^4$ | offline | 0 | Issue 6 |
| 5 | Electronics-only in-band floor measurement (device off/detuned, same looped protocol) | **rig** | ~10 min | Issue 5 |
| 6 | Hardware Δφ-response sweep vs. NMM prediction at quasi-degenerate detuning | **rig** | ~30–60 min | Section 7 (the paper's centerpiece) |
| 7 | Leave-one-speaker-out spoken digits (or single-speaker claim only) | offline | 0 | honest scope |

**Minimum evidence bar for a strong Paper B:** items 1–5 done, item 6 attempted, all headline accuracies re-reported with grouped folds + CIs, and the words "twin", "$\tau_{\text{mem}}$", "MC", "IPC", "phase-robust" used only with their corrected meanings above.

---

## 9. Claims Dictionary (keep next to the manuscript)

| Draft says | Corrected statement |
| --- | --- |
| "LTI+static-NL twin also achieves 99.4%" | Had not been tested on measured data; toy surrogate only. Now tested properly: [fill after item 1] |
| $\tau_{\text{mem}} = 4.3$–9.4 / 6.3–12.3 / 20–80 ns | Three distinct quantities: $\tau_{\text{sep}}$ (E3FIX, ≈[X]±[Y] ns), $\Delta t_{\text{th}}$ (pre-excitation horizon), symbol-overlap window |
| MC = 1.170 symbols | Input-transduction capacity 0.640 (lag 0) + conventional memory MC$_{k\ge1}$ = 0.530 |
| IPC = 1.582 | 3-symbol-window IPC (7-term Walsh), [CI], measured under grouped folds; not literature-comparable |
| IMD3 negligible | IMD3 = −31…−35 dBc but lands out of band; in-band electronics floor measured: [fill after item 5] |
| "phase-robust computation" | Phase-gated computation; Δφ-response matches the 4B cross-pump NMM prediction |

---

*Audit written by independent review pass; every number above was re-derived from `spin-waves/neuromorphic/*.py`, `spin-waves/reports/idler_math/extract_idler_params.py`, `e0_e4_metrics.json`, `multisymbol_v2_metrics.json`, `idler_extracted_params.json`, and `reports/phase-gating/REPORT.md`.*
