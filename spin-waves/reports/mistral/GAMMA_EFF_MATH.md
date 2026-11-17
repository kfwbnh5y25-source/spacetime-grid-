# GAMMA_EFF_MATH.md — Effective damping, cross-pump threshold, and phase memory

Closed-form backbone for Paper A. Every derivation here is done BY US from the linearized
normal-mode equations (the two AI-generated reports are treated as hypothesis generators only;
where we quote a number from data it is tagged VERIFIED, from theory PREDICTION).

BOSS RULE (standing): NO MODE NUMBERS anywhere in paper/figures/captions. Modes are called
$k$ (seed/pre-excitation mode) and $m$ (target/idler mode).

---

## 1. Notation and assumptions

- $b_m$: complex amplitude of mode $m$ (rotating frame), $|b_m|^2$ = modal power.
- $\Gamma_m$: intrinsic relaxation rate of mode $m$ [Hz]. From the NMM report:
  $\Gamma = 8.2$–$9.8\;\text{MHz}$ for $\alpha = 0.002$, $f_M = 4.932\;\text{GHz}$.
  UNVERIFIED simulation number, used only as scale.

  **$\Gamma$ SCALE TENSION (flagged 1117, not resolved):** three incommensurate scales are in
  play. (i) NMM simulation: $8.2$–$9.8$ MHz. (ii) Measured HWHM of the half-frequency lines:
  $0.6$–$1.6$ MHz (`extract_idler_params`, `threshold_reduction/data.h5`) — would predict
  $\tau_{\rm mem} \sim 100{+}$ ns. (iii) DYNAMIC decay rate of the phase contrast:
  $1/(2\tau) \approx 185$–$233$ MHz (delay-recheck2 + "3" dataset, §4.2) — $20$–$30\times$
  ABOVE (i) and $>100\times$ above (ii). Consequence: the observed memory decay CANNOT be the
  small-signal $\Gamma_k$ of linear theory; it must include nonlinear seed ring-down and/or
  dephasing across the mode manifold. The memory law (10) is therefore used as a DEFINITION
  ($\Gamma_{\rm eff,mem} := 1/2\tau_{\rm mem}$, measured), not as an ab-initio prediction.
  Units convention: rates in Hz $=$ cycles/s throughout ($\tau = 1/(2\Gamma)$, $\Gamma$ in
  Hz); the delay fits use $e^{-t/\tau}$ with $\tau$ in ns, so (iii) is free of any $2\pi$
  period/conversion slip.
- $f_1, f_2$: tone frequencies; tone 1 drives mode $k$ at $f_1 = 2f_k$, tone 2 drives mode
  $m$ at $f_2 = 2f_m$. Detuning $\Delta = f_2 - f_1$.
- $V_m$: parametric coupling coefficient [rate/field]. DEFINED so the off-diagonal pump
  element in the growth matrix is exactly $i V_m h$ (this sets $h_{\rm th} = \Gamma/V$;
  factor-2 differences with other papers are pure convention — see §2).
- $\Psi$: RELATIVE PHASE between the external tone-2 pump and the internal 4B cross-pump at
  the moment that matters. $\Psi = \phi - \phi_W$ with $\phi_W$ the internal pump phase.
- $t_{\rm gap}$: time gap between tone-1 pulse end and tone-2 pulse start (non-overlap).
- $\phi$: AWG-programmed relative phase of tone 2 ($t_0$-referenced, hw_core_spectro).

Assumption set (all standard, all stated so referees can attack them):
- (A1) Amplitude equations from linearized LLG/NMM are valid (small $|b|$, sub- or near-threshold).
- (A2) Single dominant seed at a time for analytic results; multi-seed enters only via $W$ as a sum.
- (A3) Thermal seeding is homogeneous; it sets the absolute floor, not the phase dependence.
- (A4) $\alpha^2$ corrections ($\alpha = 0.002$) are negligible; no back-action on $\Gamma$ assumed.

## 2. Single-mode parametric threshold: $h_{\rm th}(0) = \Gamma / V$

Lab-frame amplitude equation for mode $m$ under a pump at $f_p \approx 2 f_m$:

$$\frac{d b_m}{d t} = -\left(\Gamma_m + i\,\frac{2\pi f_p}{2}\right) b_m + (\text{noise}) \qquad \text{(lab frame)}$$

Rotating frame $b_m \to b\, e^{-i\pi f_p t}$:

$$\dot b = -(\Gamma + i\delta)\, b + i V h_p\, e^{-i\phi}\, b^*, \qquad \delta = 2\pi\!\left(f_m - \frac{f_p}{2}\right). \tag{1}$$

Growth matrix and eigenvalues:

$$M = \begin{pmatrix} -(\Gamma + i\delta) & i V h\, e^{-i\phi} \\ -i V h\, e^{+i\phi} & -(\Gamma - i\delta) \end{pmatrix}, \qquad
\lambda_\pm = -\Gamma \pm \sqrt{(V h)^2 - \delta^2}.$$

Exponential growth iff

$$(V h)^2 > \Gamma^2 + \delta^2 \quad\Longrightarrow\quad h_{\rm th}(\delta) = \frac{\sqrt{\Gamma^2 + \delta^2}}{V}, \tag{2}$$

and at exact resonance ($\delta = 0$): $\boxed{h_{\rm th} = \Gamma / V}$. This is the ONLY place
$\Gamma_{\rm eff}$ enters: any mechanism that subtracts from the effective $\Gamma$ below (2)
lowers the threshold.

Convention note: several texts write $h_{\rm th} = 2\Gamma/V$; their $V$ is HALF ours (pump
element $iVh/2$). We always use (2). The NMM report's threshold expression (its §3.1) is the
same object with $V_h$ absorbing the same convention.

## 3. Two tones: internal 4B cross-pump, $h_{\rm eff}$, and the threshold rule

### 3.1 Internal pump

Mode $k$, driven by tone 1 below ITS threshold, settles at steady power $|b_k|^2$ (VERIFIED
scaling law, fig2: $P_m \propto a_1^2 P_{f_2/2}$, log-log slopes $1.83$–$1.94$). The 4B term

$$\mathcal{F}_m^{(4B)} = \sum_{k \neq m} P_{mk}\, b_k^2\, b_m^* \tag{3}$$

is then a parametric drive on $m$ at frequency $f_1$ with amplitude
$W = |P_{mk}|\,|b_k|^2$ and phase $\phi_W = 2\theta_k + \arg P_{mk}$. In the tone-2 rotating
frame the total parametric drive is

$$A(t) = V h_2\, e^{-i\phi} + W\, e^{-i(2\pi \Delta\, t + \psi_0)}. \tag{4}$$

CRITICAL: the second term is a stationary pump only if $|\Delta| \lesssim \Gamma_m$ (the beat
must not average out). This splits the physics into the two regimes of §5.

### 3.2 Threshold rule (the closed-form the paper hangs on)

Threshold condition: $|A| = \sqrt{\Gamma_m^2 + \delta_m^2}$. With $\Psi = \phi - \psi_0$ the
relative angle between the two drives at turn-on, (4) gives a quadratic in $V h_2$:

$$(V h_2)^2 + 2 W \cos\Psi\, (V h_2) + W^2 - (\Gamma^2 + \delta^2) = 0 \tag{5}$$

$$\Longrightarrow\quad h_{\rm th}(\Psi) = \frac{\sqrt{\Gamma^2 + \delta^2 - W^2 \sin^2\Psi} - W \cos\Psi}{V}. \tag{6}$$

Linearized in $W/V$ ($W \ll \Gamma$, our regime):

$$\boxed{h_{\rm th}(\Psi) = h_{\rm th}^{(0)} - \frac{W}{V}\cos\Psi, \qquad h_{\rm th}^{(0)} = \frac{\sqrt{\Gamma^2 + \delta^2}}{V}} \tag{7}$$

Consequences, each falsifiable:
1. Phase modulation of the threshold is SINUSOIDAL with period $2\pi$ (not $\pi$) — because the
   4B pump is phase-locked to tone 2 in this regime. [Fig 3 overlay tests this.]
2. Modulation amplitude $= W/V =: h_W$, the "equivalent field" of the internal pump. Measuring
   the peak-to-valley of the threshold (or gain) vs $\phi$ MEASURES $W$ in field units — no
   tensor knowledge needed.
3. Phase-averaged threshold is UNCHANGED to linear order ($\langle\cos\Psi\rangle = 0$); the
   spread is $\pm h_W$. Phase-averaged data cannot see gating; phase-resolved data must.
4. On resonance the ignition condition $h_2 > h_{\rm th}(\Psi)$ for fixed $h_2$ defines the two
   phase states ON/OFF; contrast is largest when $h_2$ sits between $h_{\rm th}(\Psi_{\max})$
   and $h_{\rm th}(\Psi_{\min})$.

### 3.3 Phase-gating depth (the money number, corrected by the boss)

The headline claim is about GATING DEPTH: relative phase should modulate the idler response by
$\sim 50\%$ or more ("the more the better") between favorable and unfavorable phase.

Working strictly BELOW threshold at resonance, the sub-threshold idler power follows the
degenerate-parametric-amplifier gain for the phase-matched quadrature, with the total drive
$|A(\Psi)|^2 = (V h_2)^2 + W^2 + 2 V h_2 W \cos\Psi$:

$$P_m(\Psi) \;\propto\; \frac{P_{\rm seed}}{\left(1 - |A(\Psi)|/\Gamma\right)^2} \qquad\Longrightarrow\qquad \frac{P_{\max}}{P_{\min}} = \left(\frac{1 - x - w}{1 - x + w}\right)^{2}, \tag{8a}$$

with $x = V h_2/\Gamma$ and $w = W/\Gamma$ (drives taken aligned $+$ favorable). Exact, no
expansion. For $w \ll 1$: $P_{\max}/P_{\min} \approx 1 + 4w/(1-x)$ — the depth GROWS as the
operating point approaches the bare threshold ($x \to 1$) and diverges once $h_2$ exceeds the
favorable-phase threshold while staying below the unfavorable one (ignition vs dormancy:
arbitrarily large contrast).

Clean statements for the paper:
- The DEPTH is set by $W/\Gamma$ (and by how close $h_2$ sits to threshold). The threshold
  picture (threshold moves by $W/V$ across phases) and the gain picture are the same knob.
- $P_{\max}/P_{\min} \ge 1.5$ (i.e. $\ge 50\%$ depth) is achievable either with
  $4w/(1-x) \gtrsim 0.4$ sub-threshold or with a monotonically rising above-threshold operating
  point.
- NOTE: the earlier "$50\%$ = threshold reduction" reading was a misinterpretation of the
  boss's prompt; the correct claim is modulation DEPTH $\ge 50\%$, measured directly from the
  fig3 overlay fit ($A_1/A_0$ evaluated at $\phi_0$ vs $\phi_0 + \pi$). Two-point contrasts at
  the recorded phases $1.31/4.55\;\text{rad}$ are LOWER BOUNDS because neither is a cosine
  extremum ($\phi_0 \approx -155^\circ$).
- Status: fig3 overlay fit DONE (two-phase design, $\Delta\phi = 185.7^\circ \approx$
  antipodal): pooled gain ratio $A_1/A_0 = 0.011$ ($\Rightarrow |h_{4B}/h_2| \approx 0.005$),
  $\hat\phi_0 = 77.9^\circ$ CI $[13^\circ, 160^\circ]$; $R^2 = 0.010$ pooled, residual rms
  10.5% vs 1.1% modulation. VERDICT: the two-phase design CANNOT demonstrate
  smooth-cos sufficiency (residuals unstructured, $|r| < 0.001$ — so no evidence for
  multi-attractor hopping either). $\phi_0 = -155^\circ$ is EXCLUDED (OFF $-$ ON $=
  -2.15\% \pm 0.58\%$, $t = -3.7$), with the caveat that two antipodal phases make the
  $A_1 < 0$ solution degenerate and the exclusion rests on the sign convention.
  RECOMMENDED REMEASURE: 16-phase sweep at $\phi_0$, $\phi_0 + \pi$, and the two
  quadratures. The delay data (§4.2) support the depth scale independently:
  $1 - P_{\min}/P_{\max} = 49$–$60\%$ at $t_{\rm gap}=0$, $\Delta = 50$–$70$ MHz
  (definition-sensitive; $(P_{\max}-P_{\min})/(P_{\max}+P_{\min}) \approx 32$–$43\%$).

### 3.4 $K_{2B}$ and $C_m$, and the $\Gamma_{\rm eff}$ definition

Writing the threshold modification as a damping correction (compare the idler report's
self-energy $\Sigma_k = \kappa\kappa^* h_1^2 / (\Gamma_m + j\delta_m)$):

$$\boxed{\Gamma_{\rm eff}(h_1, \Psi) = \Gamma_m - C(h_1, \delta) - W(h_1)\cos\Psi} \tag{8}$$

with TWO physically distinct terms:
- $C h_1^2$: PHASE-INSENSITIVE shift, real part of the self-energy; time-averaged, survives at
  any detuning. This is the $K_{2B} h_1^2$ of the idler report.
- $W(h_1)\cos\Psi$: PHASE-SENSITIVE, coherent; requires $|\Delta| \lesssim \Gamma_m$ (or the
  burst-wise regime, §4.3).
- Self-consistency headline (idler report claim, VERIFIED at the 10–20 % level): $C_m \sim
  K_{2B}$ — the phase-insensitive and sub-threshold idler-gain descriptions agree: two
  measurements of the same tensor element scaled by $|b_k|^2$.
- NOTE the idler report's $\Sigma_k$ has NO phase dependence: it cannot explain gating. Only
  the coherently phase-locked part (§3.1–3.2) can. State this in the paper; it preempts the
  referee question "isn't this just the self-energy?"

## 4. Propagation gap: $\phi_{\rm eff}$ and the memory factor

### 4.1 Phase accumulation across the gap (derive, don't assume)

After pulse 1 ends, the seed rings freely at $f_k$ (its driven phase trajectory continues).
Tone 2 turns on at $t_{\rm on}$ with programmed phase $\phi$. What matters for gating is the
phase of the internal pump $W \propto b_k^2$ AT THE MOMENT tone 2 turns on — tone 2's
parametric drive multiplies $b_m^*$, so the beat between the two drives' phases is the angle
in (5). In the lab frame at $t_{\rm on}$:

$$\Psi(t_{\rm on}) = \big[2\pi f_2 t_{\rm on} + \phi\big] - \big[2\pi (2 f_k)\, t_{\rm on} + 2\theta_0\big], \qquad 2 f_k = f_1,$$

$$\boxed{\Psi(t_{\rm on}) = \phi + 2\pi \Delta\, t_{\rm on} + \text{const}} \tag{9}$$

Writing $t_{\rm on} = t_{\rm gap}$ (+ fixed offsets): $\phi_{\rm eff} = \phi + 2\pi \Delta\,
t_{\rm gap}$. VERIFIED (the delay_delta fit, M5): delay and phase are ONE knob. The same
formula for OVERLAP pulses reduces to the instantaneous relative phase — overlap and
non-overlap experiments are the same equation with different $t_{\rm on}$.

### 4.2 Memory: $W$ decays at rate $2\Gamma_k$

During the gap tone 1 is off; the seed amplitude decays $b_k(t) = b_k(0)\, e^{-\Gamma_k t}$.
Since $W \propto |b_k|^2$:

$$\boxed{W(t_{\rm gap}) = W(0)\, e^{-2\Gamma_k t_{\rm gap}} \quad\Longrightarrow\quad \text{phase-contrast amplitude decays at } 2\Gamma_k, \quad \tau_{\rm mem} = \frac{1}{2\Gamma_k}} \tag{10}$$

PREDICTION (before measurement): with $\Gamma = 8.2$–$9.8$ MHz (NMM scale),
$\tau_{\rm mem} = \frac{1}{2\Gamma_k} = 51$–$61$ ns.

**MEASURED (recheck of raw `delay_f1.h5` + `delay_delta.h5`, independent DFT-$|c_1|$
extraction from the 16-phase sweeps, all $(t_{\rm gap},\Delta)$; reading `ro == "1"` only
— the files also contain 4,745 `1EX` sweeps with different SA settings: span 60 MHz,
center 4142 MHz, 10 ms sweep — those must never enter any analysis):**

- The phase contrast $A_1$ is a plateau $\approx 95$–$100$ nW over $t_{\rm gap} = 0$–$2$ ns,
  then collapses to $\le 10$ nW by $\sim 5$ ns. Exponential $A_1 = a\,e^{-t/\tau} + c$
  (per-$f_1$, hot band $f_1 = 4066$–$4078$): median $\tau = 1.55$ ns [IQR 1.08–1.84]
  (fit from $2$ ns) and 1.80 ns [1.28–2.09] (fit from 1.5 ns, the fig4 window, which
  includes the plateau edge). **fig4's $\tau = 2.4 \pm 0.7$ ns is CONFIRMED** as the
  correct order of magnitude for its conditions; the memory-phase contrast is unambiguously
  GONE within $\sim 5$ ns — NOT the predicted $51$–$61$ ns.
- **BUT $\tau$ is strongly detuning-dependent** (`delay_delta.h5`, $\Sigma = 4106$ MHz,
  $\Delta = 50$–$88$ MHz, $t_{\rm gap} \le 7$ ns, all measured in randomized order):

  | $\Delta$ (MHz) | contrast $t=0$ | contrast at 4.5–7 ns | $\tau$ (ns) |
  |---|---|---|---|
  | 50–58 | depth 32–36% | persists (no decay in window) | $> 7$ |
  | 62–70 | depth 39–43% | decays | $\sim 1.5$–4 |
  | 74–78 | depth 29–39% | decays | 1.1–2.5 ($\approx$ fig4 value) |
  | 80–88 | depth 4–16% | gone by 2 ns | 0.1–0.8 |

- The $\phi_{\rm eff}$ slew $d\phi_0/dt_{\rm gap}$ (where $A_1$ is measurable, $\Delta \le 58$)
  is 74–79% of $2\pi\Delta$ — the same ~20–25% shortfall as fig4's 20.9 vs 25.2 deg/ns.
  So Eq. (9) holds in FORM but the effective detuning is reduced by
  $\delta_{\rm eff} \approx (0.26\text{–}0.37)\Delta$ — a candidate signature of mode
  pulling/hybridization, NOT a data error in (9).
- Candidate interpretations of the collapse at the $\pi$-line ($\Delta = 74.5$ MHz, the
  $c_2$ identity-switch line): (i) the decay we see is the NONLINEAR ringdown of the
  over-threshold driven seed (extra nonlinear damping channels), not small-signal
  $\Gamma_k$ — most likely for the $\Delta \approx 74$ data; (ii) near the mode-pair
  crossing the hybrid eigenmode damping changes (levelling the tension at low $\Delta$,
  where the $>7$ ns lower bound is COMPATIBLE with $51$–$61$ ns). Honest phrasing for the
  paper: the memory time is detuning-dependent, $\tau_{\rm mem}(\Delta)$, spanning from
  $> 7$ ns at $\Delta \le 58$ MHz to $< 1$ ns at $\Delta \ge 80$ MHz; the simple
  $e^{-2\Gamma_k t}$ law (10) is the LOW-$\Delta$ limit, and the collapse toward the
  identity-switch line is an open measurement we propose (16-phase sweeps at
  $\Delta = 60$, 62, 64 with per-phase repeat averaging to beat the baseline-recalibration
  steps).

**MEASURED (delay-recheck2, independent re-extraction 1117 of the same files; per-$(t_{\rm gap},\Delta)$
16-phase sweeps with noise tracking):**
- $\Delta = 74$: $A_1$ decay $\tau = 2.13 \pm 0.08$ ns, but NOT a clean exponential — a
  Gaussian-onset fit form does better and a secondary plateau (8–16 µW) remains at 8–16 ns
  (a residual slow component; report honestly, do not describe the decay as single-$\tau$).
- $\Delta = 61$: NO decay through 7.4 ns — memory persists at half-contrast to the window
  edge (compatible with the low-$\Delta$ limit of (10), and now the sharpest statement of
  the low-$\Delta$ memory scale).
- Slew at $\Delta = 74$: $20.7 \pm 0.1$ deg/ns $= 77.7\%$ of $2\pi\Delta$ — the earlier
  "zero-slew anomaly on the $\pi$-line" was a fit-window artifact and is RETRACTED; Eq. (9)
  holds in form at all measured $\Delta$ with the $(0.26\text{–}0.37)\Delta$ shortfall.
- Per-sweep noise $\sigma$: 5.3 µW ($\Delta=61$), 7.0 µW ($\Delta=74$); step-corrupted
  blocks excluded: 8.0% ($\Delta=74$), 3.1% ($\Delta=61$).

**MEASURED ("3"-dataset + July map, 1117):**
- "3" = $10\times$ repeat pulse/gap trains at 50% duty, $\tau_p = 0.31$–$9.2$ ns:
  threshold $a_{\rm th}$ falls 590 $\to$ 320 mV and SATURATES for $\tau_p \ge 8$ ns
  $\Rightarrow$ population build-up completes within $\sim 8$ ns of tone-2 drive; dynamic
  ring-up $\tau_b = 2.15 \pm 0.08$ ns $\Rightarrow 1/(2\tau_b) \approx 232$ MHz, EXCLUDING
  the NMM 8–10 MHz damping scale (feeds the §1 $\Gamma$-tension note).
- July $a_1 \times f_1$ map: $A_1$ peaks at $\Delta = 70$–$80$ MHz (NOT 50–58) under the
  July settings — the strong-zone position moves with the operating point; treat
  "optimal $\Delta$" as setup-dependent.
- $c_2$ sign flips at $\Delta \approx 76 \pm 2$ and $\approx 103 \pm 2$ MHz under
  $\Sigma \approx 4113$ MHz: TWO identity-switch lines bracketing the old 74.5 MHz line
  ($\Sigma = 4106$) — the jump line MOVES with $\Sigma$ (i.e. with tone-1 frequency), it is
  not a fixed detuning scalar.

**MEASURED (1117 sanity, edge-ramp effect — new physics knob):**
- With 32 ns tukey edges the idler line PERSISTS ($A_0 = 183$ µW, raw-h5 re-extraction:
  $P(\phi) = 183 \pm 2$ µW flat, $A_1 \le 1$ µW, $A_1/A_0 \le 0.005$; live verdict agreed:
  S3a = 0.00, S3c slew 63.7 deg/ns = noise) but the phase modulation is GONE. With 1 ns
  edges (June/`delay2` geometry) the live verdict recovered: $A_0 \approx 170$ µW,
  $A_1/A_0 = 0.34$, decay ratio $A_1(8)/A_1(0.5) = 0.37$ vs June anchor 0.12, slew
  $10.0$ deg/ns vs June 20.7 (QC of the GO-run h5 PENDING — not yet copied off the setup
  PC; the archived local h5 is the pre-fix NO-GO file).
- Interpretation candidate: the turn-on edge width SMEARS the phase reference that fixes
  $\Psi$ at tone-2 turn-on (GEM 9 assumes a sharp $t_{\rm on}$); ramp $\gg$ build-up scale
  $\Rightarrow$ phase averaging washes the gating. Standing rule: 1 ns edges (June scripts'
  64-sample tukey). meas5 (edge sweep 0.5–10 ns, overnight) measures the falloff shape;
  treat its 2–10 ns "expectations" as interpolation scaffolding, NOT theory.
- S1 caveat: excitation-path peak levels ran 18–27 dB above the morning probe references
  (f1 $-32.5$ vs $-59.1$ dBm; f2 $-25.6$ vs $-43.3$ dBm) — both lines robustly seen;
  sample-path physics unaffected; check cables/attenuation/level-currency of the probe
  reference next time the EX path is used.

**MEASURED (meas5, edge sweep 0.5–10 ns, 1117 21:40, OPERATING POINT = sanity OP
$\Delta = 74$ MHz; numbers read off `meas5_edges.png`, h5 re-fit PENDING morning):**
- $A_1/A_0$ vs $t_{\rm edge}$ (ns): 0.38 (0.5), 0.35 (1), 0.19 (2), 0.045 (4), 0.073 (7),
  0.042 (10); end-control 1 ns block re-run at the end reproduces 0.35 (no drift).
  $A_1$ itself: 63 $\to$ 60 $\to$ 36 $\to$ 9 $\to$ 14.5 $\to$ 9 µW while $A_0$ stays
  $\approx$ constant — the modulation dies, the line does not.
- **VERDICT: Outcome A.** The turn-on ramp is the phase-smearing knob: gating survives
  edges $\lesssim 1$ ns, halves by $\sim 2$ ns, floor ($\approx 5$–$7\%$, non-zero
  residual!) by 4 ns. Interpolation scaffold of meas5 matched the data; the 32 ns point
  ($A_1/A_0 \le 0.005$) sits on the same floor.
- The falloff scale ($\sim 2$ ns) coincides with BOTH the memory $\tau(2.13$ ns at
  $\Delta\!=\!74)$ and the "3"-dataset ring-up $\tau_b = 2.15$ ns — consistent with the
  smearing picture: the ramp must be sharp compared with the seed dynamics for $\Psi$ to be
  a well-defined turn-on phase (GEM 9 assumes delta-function $t_{\rm on}$). The nonzero
  floor above 4 ns is worth one sentence + a measurement Jurassic candidate (residual
  overlap of tone-2 tail with ringing seed? $Q$-smearing? report as open).

Data-quality note (both files): baseline recalibrates every $\sim 90$ s while a
(delay, $f_1$) phase block takes $\sim 8$ s; $\sim 9\%$ of blocks carry a mid-block step
that fakes large $|c_1|$. Clean blocks are identifiable by cos-residual $\gg$ noise;
all numbers above come from manual inspection of clean hot-band blocks plus robust medians.

Second caveat for the theory figure: after subtracting $\phi_{\rm eff}$ (9), the residual phase
response must be STATIONARY in $t_{\rm gap}$ — the apparent $\phi_0$ must follow
$\phi + 2\pi\Delta\, t_{\rm gap}$ exactly (it did, M5, VERIFIED). Any deviation = additional
slow physics (e.g. heating drift) to report honestly.

### 4.3 Far-detuned CW tension (open problem the paper CAN claim)

Phase-gating report (AI): for $|\Delta f| \ge 50$ MHz the 4B pump averages OUT in steady
state — yet f1_4069_5 ($\Delta \approx 58$–$64$ MHz, overlap CW) shows a clean rank-1
$\cos\phi$ response with large contrast. Burst-wise hypothesis: the pulsed acquisition
sequence re-selects the attractor basin EVERY burst; each burst's outcome is a probabilistic
function of the instantaneous phase at tone-2 turn-on; the SA integrates probability-weighted
averages over the burst train. If confirmed (pulse-width tomography, report Protocol 3), the
paper contributes the steady-state-vs-burst regime map itself: the 4B pump needs
$|\Delta| \lesssim \Gamma$ for CW coherence, but burst-wise re-selection restores phase
sensitivity at ANY $\Delta$ — bounded only by the attractor-selection probability (the
far-detuned multi-attractor hopping of the report becomes the noise floor of the gating).

---

## 5. Regime map and figure mapping (one look-up, nothing more)

| Regime | Condition | Physics | Figure |
|---|---|---|---|
| Single tone | one pump near $2f_m$ | $h_{\rm th} = \Gamma/V$, Arnold tongue | Fig 1 (overview map, zoom, hinge) |
| Overlap CW, near-degenerate | $\Delta \lesssim \Gamma$ | 4B locked: Eq. (7) | Fig 3 overlay |
| Overlap CW, far-detuned | $\Delta \gg \Gamma$, pulsed | burst-wise re-selection; memory via (9) | Fig 3 + open-tension text |
| Non-overlap (gap) | $t_{\rm gap} > 0$ | $W = W_0 e^{-2\Gamma_k t_{\rm gap}}$; $\phi_{\rm eff}$ (9) | Fig 4 (delay axis) |
| Non-overlap fine-detune map | $\Delta$ 56–89 MHz, 16 phases | $c_1$ dominates strong zone; $c_2$ $\pi$-jump at $\Delta = 74.5$ MHz (identity switch, VERIFIED cross-dataset) | Fig 5 |
| Scaling / maps | $a_1$ vs $a_2$ 2D | $P_m \propto a_1^2 P_{f_2/2}$; $\Gamma_{\rm eff}$ shift $C h_1^2$ | Fig 2 |

Synthesis figure (LAST theory figure, PREDICTS the Fig 1/2 contours): overlay of (6)+(8)+(10):
$\Gamma_{\rm eff}(h_1, \Psi, t_{\rm gap}) \to$ predicted threshold surface
$h_{\rm th}(f_1, f_2, \phi, t_{\rm gap})$ with the measured points from Figs 2–5 plotted on
top. ONE equation per axis, no numerics in main text.

## 6. Open items / status ledger

- VERIFIED: $P_m \propto a_1^2 P_{f_2/2}$ (slopes 1.83–1.94); $\phi_{\rm eff} = \phi +
  2\pi\Delta\, t_{\rm gap}$ (M5); $c_2$ $\pi$-jump line $\Delta = 74.5$ MHz (fig5 + detune
  copy 2, cross-dataset); $\phi_0 = -155^\circ \pm 4^\circ$ (031126 data); $c_1$ continuous
  across the jump line; $c_1$ dominant in the strong zone.
- THEORY → NUMERICALLY SELF-TESTED (tmp/test_gamma_eqs.py v2, ALL PASS): (2), (6), (7),
  (8a) — exact eigenvalue bisection, exact algebraic steady state, plus unbiased
  RK4 cross-check (asymptotic log-slope). v1 test failures were harness bugs
  (fixed 5 us detection window = growth-rate contour, +1 sign error in the real form).
  Still THEORY for the physics mapping: (8) full $\Gamma_{\rm eff}$, (10) memory law
  (see MEASURED below — the memory law FAILS quantitatively near $\Delta\approx74$).
- MEASURED (delay recheck, raw-data independent): fig4 $\tau = 2.4 \pm 0.7$ ns CONFIRMED
  in its regime; $\tau_{\rm mem}(\Delta)$ detuning-dependent: $>7$ ns (no decay in window)
  at $\Delta \le 58$, 1.1–2.5 ns at $\Delta = 74$–$78$, $<1$ ns at $\Delta \ge 80$;
  $\phi_{\rm eff}$ slew = 74–79% of $2\pi\Delta$ (both datasets); contrast
  $1 - P_{\min}/P_{\max} = 49$–$60\%$ at $t_{\rm gap} = 0$ (the >=50% depth money number,
  reproduced with the $1 - P_{\min}/P_{\max}$ definition). CAUTION: 1EX sweeps in the same
  files must be excluded; some blocks carry baseline-step artifacts.
- RETRACTED: "zero slew on the $\pi$-line" — fit-window artifact; slew $= 77.7\%$ of
  $2\pi\Delta$ at $\Delta = 74$ (delay-recheck2, §4.2).
- IDLER_MATH_REPORT $\S$2.2 INDEPENDENTLY RE-DERIVED (1117, by us, from the 3-wave coupled
  pair; d'Aquino (42)-(43) $\to$ rotating frame $\to$ adiabatic elimination of $c_m^*$,
  valid for $|\kappa h_1| \ll \sqrt{\Gamma_m^2+\delta_m^2}$): the report's
  $\Sigma_k = \kappa_{km}\kappa_{mk}^* h_1^2/(\Gamma_m + j\delta_m)$,
  $\Gamma_{k,\rm eff} = \Gamma_k - \Re\Sigma\,h_1^2$, $\delta'_k = \delta_k -
  \Im\Sigma\,h_1^2$, and $K_{2B} = \Re(\kappa_{km}\kappa^*_{mk})\,\Gamma_m /
  [\Gamma_k(\Gamma_m^2+\delta_m^2)]$ are ALL CORRECT. Caveats: (a) its closed threshold
  keeps only the real branch — exact form is $\sqrt{\Gamma_{\rm eff}^2 + \delta'^{\,2}}/|V_k|$
  (fine at $\delta_k \approx 0$, state it); (b) $\phi_1$ cancels in
  $\kappa_{km}\kappa_{mk}^*$ $\Rightarrow$ self-energy is PHASE-INSENSITIVE, consistent with
  the GEM (8) split — the report cannot explain gating, as already noted in $\S$3.4.
- IDLER $\S$3 (p-channel) VERIFIED (1118, sympy, tmp/check_idler3.py): the upper-side
  channel is the EXACT MIRROR of the idler channel under $1\leftrightarrow2$,
  $k\leftrightarrow h$: tone 2 scatters tone-1-half magnons into
  $f_p = f_2 - f_1/2$; the resonant rotating term for $c_p^*$ is slow iff
  $\omega_p = \omega_2-\omega_1/2$, same algebraic structure as the $m$ channel — hence
  $P_p \propto a_2^2\,P_{f_1/2}$ with $C_p$ the same coupling object. The report's
  "why $P_p$ is $50\times$ weaker" logic (same $C$, smaller seed $|c_h|^2$) is
  internally consistent.
- IDLER $\S$6 self-consistency VERIFIED (1118, sympy, same script):
  $C_m(\beta_{mk}{=}1) \equiv K_{2B}$ exactly
  ($K_{2B}/C_m = \beta_{mk}^{-2}$ in general, both $\propto |\kappa|^2 k_{\rm eff}^2 /
  (\Gamma_m^2+\delta_m^2)$ with $\kappa_{km}=\kappa_{mk}$, $\delta_m$ common). The
  measured 10–20% agreement is therefore a genuine two-measurements-one-tensor
  check, not an algebraic accident.
- IDLER $\S$4 quote-check (1118, pdftotext p.10): d'Aquino's amplitude-dependent detunings
  $\epsilon'_h = \epsilon_h - 2\omega_h N_{hn}|a_n|^2$, $\epsilon'_n = \epsilon_n -
  2\omega_n N_{nh}|a_h|^2$ VERIFIED as PDF equations (51) — the idler report's
  $\epsilon'_k = \epsilon_k - 2\omega_k N_{kh}|a_h|^2$ is the same object (relabelled).
  4B cross-pump $W_{4B,k} = \sum_{h\neq k} P_{kh}a_h^2$ matches GEM §3.1 form.
- IDLER $\S$5 (NMM ODE verification) RE-RUN LIVE 1118 (node,
  reports/idler_math/verify_idler_nmm.mjs): reproduces the printed table EXACTLY —
  RWA-fast $P_{17}\approx2.1\times10^{-8}$ flat while full-NMM
  $\Delta P_{17}/(h_1^2P_{28}) = 0.00549/0.00605/0.00574/0.00547\ \text{mT}^{-2}$
  across $16\times$ in $h_1^2$; D-vs-ND total power $0.0059 > 0.0038$ (non-commutative
  direction as claimed). Internal verification only — mode indices NEVER leave this
  ledger (BOSS RULE).
- THRESHOLD-REDUCTION NUMBER (was hard to find; make it explicit): the idler report states
  it only in RELATIVE form, $h_{2,\rm th}(h_1) = h_{2,\rm th}(0)\,(1 - K_{2B}h_1^2)$
  (its $\S$2.2, last display). Absolute shift at the working point
  $a_1 = 327$ mV: $K_{2B} a_1^2 = (0.94\text{–}1.36)\times10^{-6}\,\text{mV}^{-2}\times
  (327\,\text{mV})^2 = 10\text{–}15\%$ reduction. The paper should state BOTH the relative
  form and this absolute number. Do NOT confuse with the depth claim: $P_{\max}/P_{\min}$
  ($\ge 50\%$) is the GATING DEPTH (§3.3), a different observable.
- PENDING (downgraded): mode-number language stays out (BOSS RULE) until the eigenmode vs
  bias-field fit is done (simulated eigenmodes on boss's PC; bottom of todo).
- BOSS RULE: no mode numbers until the eigenmode vs bias-field fit is done (slow, bottom of
  the list). The "identity switch" language is deliberate.
- LaTeXSTYLE: equations in every .md use LaTeX (works: pandoc/tectonic later).
- CITATION CONFIRMED (web, 2 independent hits): d'Aquino, Perna, Merbouche, de Loubens,
  "Nonlinear interaction theory for parametrically-excited spin-wave modes in confined
  micromagnetic systems", arXiv:2601.07798 (v2, cond-mat.mes-hall). Matches the local PDF.
- BONUS FIND (cite!): arXiv:2601.05775, Munoz, ..., Anane, Perna, d'Aquino, Serpico,
  Merbouche, de Loubens, "Nonlinear mode interactions under parametric excitation in a
  magnetic microdisk" — EXPERIMENTAL two-tone parallel pumping in a YIG microdisk; steady
  state depends on mode pair, pump detuning, and TEMPORAL SEQUENCE of the two tones
  ("non-commutative behaviors"). Direct prior/companion experiment for our Fig 3/4 results;
  must be cited and differentiated before the LaTeX draft.

## Verified citations (page + equation)

All numbers below verified against the local PDFs with `pdftotext -f N -l N -layout` (per-page
extraction). For BOTH papers the printed page number equals the PDF page number (d'Aquino:
50-pp arXiv preprint, page number printed in footer; Perna: 12-pp JMMM article, page 1
unnumbered, pages 2-12 printed). Quotes are verbatim from the layout extraction.

### d'Aquino, Perna, Merbouche, de Loubens 2026 — arXiv:2601.07798v2, "Nonlinear interaction theory for parametrically-excited spin-wave modes in confined micromagnetic systems" (50 pp)

NOTATION WARNING: the paper uses mode indices $h, n$ (not $k, m$) and coefficients
$N_h$ (self-NFS) / $N_{hn}$ (mutual-NFS) — the strings "$N_{hk}$", "$n_{mj}$", "$P_{mk}$" DO
NOT APPEAR anywhere in the PDF. Index.html's cited equation numbers map onto the paper as follows.

**index.html claim check:**
- "Eqs. (55)-(60)": EXISTS but is NOT the amplitude equation of motion — (55)-(56) are the
  coupled steady-state amplitude formulas $|a_h|^2_{\pm\pm,\mathrm{coup}}$, $|a_n|^2_{\pm\pm,\mathrm{coup}}$
  [p. 10]; (57)-(60) define the four steady-state regimes z, $u_h$, $u_n$, $c_{\pm\pm}$ [p. 11].
  The amplitude equations WITH pump terms are actually (42)-(43) [p. 9] and their averaged
  form (49)-(50) [p. 10].
- "Table I Rules 1-12": EXISTS — Table I, p. 14, caption: "Summary of qualitative features and
  algorithm to construct the phase diagram in the $(\omega_{\mathrm{rf}1},\omega_{\mathrm{rf}2})$
  control plane..."; text (p. 14): "The rules 1-4 refer to existence whereas the rules 5-12 refer
  to stability of interacting steady-state regimes."
- "Routh-Hurwitz Appendix C.3": EXISTS — Appendix C, subsection 3 "P-mode stability under
  two-tone excitations" (heading p. 30; body refs "appendix C 3", pp. 12, 13, 21). The words
  "Routh-Hurwitz" appear verbatim p. 33: "Among various technical conditions, the Routh-Hurwitz
  criterion states that a necessary condition for the stability of P-mode $c^{++}$ is that all
  the coefficients of the polynomial are positive." Resulting stability conditions (C30)-(C32),
  p. 34; Hopf condition $S_2 = S_3/S_1 + S_4 S_1/S_3$, p. 34.
- "Eq. (33)/(34) single-mode PWM steady state": EXISTS and is correct — see below.
- "Eq. (9) 4B term": PARTIALLY correct — Eq. (9) contains the third-order (four-magnon) term
  $\sum_{i,j,k} d_{hijk}\, a_i a_j a_k / 2$, but no explicit cross-pump form $P_{mk} b_k^2 b_m^*$;
  the quartic coupling is packaged later as the mutual NFS matrix, Eq. (40).

**1. Amplitude equation of motion for a mode amplitude with pump terms:**
- Eq. (9), p. 4 — Normal Modes Model amplitude equation (truncated at third order):
  $\dot a_h = j\omega_h\Big(b_h + \sum_i b_{hi} a_i + \sum_{i,j} c_{hij} a_i a_j +
  \sum_{i,j,k} \tfrac{d_{hijk}}{2} a_i a_j a_k\Big)$;
  "In eq.(9), the third and fourth terms resemble those arising from three- and four-magnon
  scattering interactions [43]" (p. 4).
- Eq. (31), p. 6 — single-mode amplitude equation with pump:
  $\dot a_h = j\omega_h\big[(1 + N_h|a_h|^2 - \delta h_a(t)\|\varphi_h\|^2)a_h
  - \delta h_a(t)(\varphi_h,\varphi_h^*)\, a_h^*\big] - \lambda_h a_h$, and c.c.
- Eqs. (42)-(43), p. 9 — two-modes model (amplitude equations with the two pump tones):
  "$\dot a_h = j\omega_h[1 + N_{hh}|a_h|^2 + N_{hn}|a_n|^2 - \delta
  h_a(t)\|\varphi_h\|^2]a_h + \frac{j\omega_h}{2}\delta h_a(t)(\varphi_h,\varphi_h^*)a_h^*
  - \lambda_h a_h$, and c.c." (42)
- Eqs. (49)-(50), p. 10 — averaged two-modes model in slow variables $(A_h,\Phi_h,A_n,\Phi_n)$:
  "$\dot A_h = -[\frac{\omega_h(\varphi_h,\varphi_h^*)\delta h_1}{2}\cos(2\Phi_h) - \lambda_h
  A_h]$, (49)"; $\dot\Phi_h + \frac{\epsilon_h}{2} = \omega_h(N_{hn}A_n^2 + N_{hh}A_h^2) +
  \frac{\omega_h(\varphi_h,\varphi_h^*)\delta h_1}{2}\sin(2\Phi_h)$ (both inside (50)).

**2. Parametric threshold (pump amplitude vs damping and detuning):**
- Eq. (20), p. 5 — minimum threshold at zero detuning:
  $\delta h^{\mathrm{thres},h} = \dfrac{2\lambda_h}{\omega_h|(\varphi_h,\varphi_h^*)|} =
  \dfrac{2\alpha\omega_h\|\varphi_h\|^2}{|(\varphi_h,\varphi_h^*)|}$; with $\lambda_h =
  \alpha\omega_h^2\|\varphi_h\|^2$ defined p. 5 below Eq. (14).
- Eq. (21), p. 5 — threshold in terms of ellipticity:
  $\delta h^{\mathrm{thres},h} = 2\alpha\omega_h\frac{1+e_h^2}{1-e_h^2}$.
- Eq. (24), p. 6 — detuning cone (full threshold vs detuning):
  $\delta_h > \sqrt{\frac{(\omega_{\mathrm{rf}}-2\omega_h)^2}{\omega_h^2
  |(\varphi_h,\varphi_h^*)|^2} + (\delta h^{\mathrm{thres},h})^2}$ — same
  $\sqrt{\Gamma^2+\delta^2}/V$ structure as Eq. (2) of THIS document (V $\leftrightarrow$
  $\omega_h|(\varphi_h,\varphi_h^*)|/2$; note this convention matches ours, no factor 2).
- Eqs. (25)-(26), p. 6 — critical detuning $\epsilon_{\mathrm{crit},h}(\delta_h) =
  \omega_h|(\varphi_h,\varphi_h^*)|\sqrt{\delta_h^2 - (\delta h^{\mathrm{thres},h})^2}$ and
  parametric excitation condition "$|\epsilon_h| = |\omega_{\mathrm{rf}} - 2\omega_h| <
  \epsilon_{\mathrm{crit},h}(\delta_h)$. (26)". Eq. (32), p. 7 — steady-state version
  $(\omega_{\mathrm{rf}} - 2\omega_h')^2 = \epsilon_{\mathrm{crit},h}^2$.

**3. Four-magnon nonlinear frequency-shift coefficients and the "4B cross-pump":**
- Eqs. (9)-(10), p. 4 — four-magnon tensor definition:
  $d_{hijk} = -(\varphi_h, \psi_{jk}C\varphi_i) + (\varphi_h, \mathbf{m}_0\cdot C\psi_{jk}\,
  \mathbf{m}_0\varphi_i)$, with $\psi_{ij} = \varphi_i\cdot\varphi_j$.
- Eqs. (28)-(29), p. 6 — self-NFS (nonlinear frequency shift, 4B self-term):
  "= $j\omega_h N_h|a_h|^2 a_h$, (28)" with $\tilde N_h = d_{h,h,-h,h} +
  \frac{d_{h,-h,h,h}}{2}$ (29); used in the steady-state equation via $\omega_h'(|a_h|) =
  \omega_h(1 + N_h|a_h|^2)$ (27), p. 6.
- Eqs. (40)-(41), p. 9 — mutual (cross) NFS matrix — the index.html "$P_{mk}$"
  cross-pump coefficient corresponds to THIS object:
  "$\tilde N_{hn} = d_{h,n,h,-n} + d_{h,h,n,-n} + d_{h,-n,h,n}$ (40)", text: "the real part
  of the coefficient $\tilde N_{hn}$ would act as a pure mutual nonlinear frequency shift ...
  the imaginary part would play the role of a nonlinear damping-like term"; coefficients
  $N_{hn} = \Re\{\tilde N_{hn}\}$, $\lambda_{hn} = \omega_h\Im\{\tilde N_{hn}\}$ (unnumbered,
  just below (40)), symmetry "$\lambda_{hn} = \lambda_{nh} = 0,\ N_{hn} = N_{nh}$ (41)".
  NO coefficient named $P_{mk}$ exists in the PDF (checked: "P_{mk}"/"cross-pump" 0 hits);
  the explicit quartic cross coupling is $N_{hn}|a_n|^2 a_h$ inside (42)-(43), p. 9.
- Coupled steady-state amplitudes, Eqs. (55)-(56), p. 10:
  "$|a_h|^2 = L_{hn}\big[|a_h|^2_{\pm,\mathrm{unc}} - \frac{N_{hn}}{N_{hh}}|a_n|^2_{\pm,
  \mathrm{unc}}\big]$ (55)" with $L_{hn} = L_{nh} = \frac{N_{hh}N_{nn}}{N_{hh}N_{nn} -
  N_{hn}^2}$ (defined p. 11); regime list (57)-(60), p. 11: "$z: |a_h|^2 = 0, |a_n|^2 = 0$
  (57)"; $u_h$ (58); $u_n$ (59); "$c = c_{\pm\pm}: |a_h|^2 = |a_h|^2_{\pm\pm,\mathrm{coup}},
  |a_n|^2 = |a_n|^2_{\pm\pm,\mathrm{coup}}$ (60)".
- Table I, p. 14 — rules 1-12 (rules 1-4 existence, 5-12 stability), e.g. rule 6: "the coupled
  P-modes $c^{++}$, given by (60), are stable in $T^{++}$ close to the lines $c^{++}\equiv u_h$
  and $c^{++}\equiv u_n$ if $\mathrm{sign}(N_{hh})\mathrm{sign}(N_{nn})(N_{hh}N_{nn} -
  N_{hn}^2) > 0$".
- Routh-Hurwitz analysis in Appendix C 3, p. 33: eigenvalues $\zeta_{1..4}$ of the 4x4 Jacobian
  (C23), p. 31; "$S_4 = \zeta_1\zeta_2\zeta_3\zeta_4 = \det(J)$" p. 33; stability conditions
  (C30)-(C32), p. 34: "$c^{++}$ is stable close to $c^{++}\equiv u_h$ or $c^{++}\equiv u_n$ in
  $T^{++}$ if $s_h^+ s_n^+(N_{hh}N_{nn} - N_{hn}^2) > 0$, (C30)". Caveat stated p. 34:
  "the above conditions guarantees the instability ... but the converse is not true since they
  have been derived from a necessary condition for stability."

### Perna, Bruckner, Serpico, Suess, d'Aquino 2022 — J. Magn. Magn. Mater. 546 (2022) 168683, "Computational micromagnetics based on normal modes: Bridging the gap between macrospin and full spatial discretization" (12 pp; journal, DOI 10.1016/j.jmmm.2021.168683)

**1. Normal-mode projection formalism (micromagnetic modes -> oscillator amplitudes):**
- Eq. (2), p. 2 — modal expansion of transverse magnetization:
  $\delta\mathbf{m}_\perp(\mathbf{r},t) = \sum_h a_h(t)\,\boldsymbol\varphi_h(\mathbf{r})$, (2);
  text (p. 2): "the complex amplitudes $a_h(t)$ are the degrees of freedom of the magnetization
  dynamics and their evolution is mathematically formulated by the following system of coupled
  nonlinear ordinary differential equations, here and on referred to as Normal Modes Model (NMM)".
- Eq. (3), p. 2 — the NMM amplitude equation (conservative, with forcing via $b_h$, up to
  5-amplitude nonlinear terms):
  $\dot a_h = \frac{j\omega_h}{1+\alpha_G^2}\Big(b_h + \sum_i b_{hi}a_i +
  \sum_{i,j} c_{hij}a_i a_j + \frac{1}{2}\sum_{i,j,k} d_{hijk}a_i a_j a_k +
  \frac{1}{4}\sum_{i,j,k,l} e_{hijkl}a_i a_j a_k a_l + \frac{1}{4}\sum_{i,j,k,l,m}
  f_{hijklm}a_i a_j a_k a_l a_m + \dots\Big)$, (3);
  note (p. 2): "The highest order nonlinear terms involve the product of 5 complex amplitudes."
- Projection identity, Eq. (20), p. 9 — how the projection works:
  "$(\varphi_h, \frac{\partial\delta\mathbf{m}_\perp}{\partial t})_{A^\perp_0} =
  \sum_k (\varphi_h, \varphi_k)_{A^\perp_0} \dot a_k = \dot a_h$" leading to the conservative
  single-mode ODE "$\dot a_h = j\omega_h a_h$. (21)", p. 9.

**2. Normalization / orthogonality relation:**
- Eq. (18), p. 9 (also Table 1, p. 3): $(\boldsymbol\varphi_h,\boldsymbol\varphi_k)_{A^\perp_0}
  = (\boldsymbol\varphi_h, A^\perp_0\boldsymbol\varphi_k) = \frac{1}{V_\Omega}\int_\Omega
  \boldsymbol\varphi_h^*\cdot A^\perp_0\boldsymbol\varphi_k\, dV = \delta_{hk}$, (18);
  self-adjoint equivalent form Eq. (19), p. 9: $(\varphi_h,\varphi_k)_{A^\perp_0} =
  (\!A^\perp_0\varphi_h,\varphi_k) = j\omega_h(\mathbf{m}_0\times\varphi_h,\varphi_k) = \delta_{hk}$;
  generalized eigenproblem Eq. (16), p. 9: $A^\perp_0\boldsymbol\varphi = \omega B_0
  \boldsymbol\varphi$, with $B_0\boldsymbol\varphi = -j\mathbf{m}_0\times\boldsymbol\varphi$.

**3. Mode damping (the $\Gamma_m$ / $K_p$-analog):**
- NO modal damping rate "$\Gamma_m$" or $K_p$ is defined in Perna 2022. Damping enters through
  $\alpha_G$-proportional coefficients: "$b_{hi} = \delta_{hi} + \alpha_G b^\alpha_{hi}$, ..."
  (Table 1, p. 3; Eqs. (41), p. 10, with $b^\alpha_{hi} = j(\,\omega_i\,\varphi_h,\varphi_i\,)$),
  and the projected damping torque, p. 10: "we can derive the projection of the damping torque
  in the normal modes basis expressed by the following terms" (Eqs. (37)-(38) p. 10).
- Modal decay TIME (the $\tau_m = 1/\Gamma_m$ analog), Eq. (5), p. 5 — closed-form linear
  driven/free response under the diagonal approximation "$b_{hi} \approx \delta_{hi}(1 +
  b^\alpha_{hi})$" (p. 5):
  $a_h(t) = A_h e^{-t/\tau_h}(\cos\omega_h t + j\sin\omega_h t) + \frac{\omega_{\mathrm{rf}}}
  {j\omega_h - 1/\tau_h}(\sin\omega_{\mathrm{rf}}t - \cos\omega_{\mathrm{rf}}t)$,
  $A_h = \frac{-\omega_h^2 B_h}{\omega_{\mathrm{rf}}^2 - \omega_h^2 - 2j\omega_h/\tau_h}$,
  $\boxed{\tau_h = |b^\alpha_{hh}|^{-1} = \big(\alpha_G\,\omega_h^2\|\varphi_h\|^2\big)^{-1}}$
  — i.e. modal damping rate $\lambda_h = \alpha_G\omega_h^2\|\varphi_h\|^2$, consistent with
  d'Aquino 2026's $\lambda_h$ below Eq. (14), p. 5.

**4. LLG-to-amplitude reduction (full damped + forced NMM in Methods):**
- Eq. (31), p. 10 — conservative nonlinear NMM (parabolic approximation): "$\dot a_h = j\omega_h
  a_h + j\omega_h\sum_{i,j} c^0_{hij}a_i a_j + \frac{j\omega_h}{2}\sum_{i,j,k} d^0_{hijk}a_i
  a_j a_k + \frac{j\omega_h}{4}\sum_{i,j,k,l} e^0_{hijkl}a_i a_j a_k a_l$, (31)".
- Eqs. (33)-(34), p. 10 — LLG with Gilbert damping: "$\partial_t\mathbf{m} - \alpha_G\mathbf{m}
  \times\partial_t\mathbf{m} = -\mathbf{m}\times\mathbf{h}_{\mathrm{eff}}[\mathbf{m}]$, (33)";
  "$(1+\alpha_G^2)\partial_t\mathbf{m} = -\mathbf{m}\times\mathbf{h}_{\mathrm{eff}}[\mathbf{m}]
  - \alpha_G\mathbf{m}\times(\mathbf{m}\times\mathbf{h}_{\mathrm{eff}}[\mathbf{m}])$, (34)".
- Eq. (40), p. 10 — damped nonlinear NMM (generalizes (31)): "$(1+\alpha_G^2)\dot a_h =
  j\omega_h\big(\sum_i b_{hi}a_i + \sum_{i,j} c_{hij}a_i a_j + \frac{1}{2}\sum_{i,j,k}
  d_{hijk}a_i a_j a_k + \frac{1}{4}\sum_{i,j,k,l} e_{hijkl}a_i a_j a_k a_l + \frac{1}{4}
  \sum_{i,j,k,l,m} f_{hijklm}a_i a_j a_k a_l a_m\Big)$, (40)".
- Eqs. (44)-(45), p. 11 — projection of forcing ($\mathbf{f}_c = \mathbf{h}_{\mathrm{rf}} +
  \alpha_G\beta\mathbf{p}$, $\mathbf{f}_d = \alpha_G\mathbf{h}_{\mathrm{rf}} - \beta\mathbf{p}$)
  onto the mode basis; Eq. (46), p. 11 — full forced+nonlinear+damped NMM; pairwise projection
  "$b_h \approx -(\varphi_h, \mathbf{f}_c) = B_h\cos(\omega_{\mathrm{rf}}t)$" unnumbered above
  Eq. (4), p. 4.

**Explicitly NOT verifiable / NOT present (do not cite):**
- d'Aquino 2026: any equation numbered (62) or beyond — the last numbered equation found in the
  main-text numbering used here is (60) p. 11 (appendices restart at (A1),(B1),(C1); highest
  appendix numbers seen: (A10) p. 27, (B5) p. 28, (C32) p. 34).
- d'Aquino 2026: coefficient symbols "$P_{mk}$", "$N_{hk}$", "$n_{mj}$" — not in the PDF.
- Perna 2022: parametric threshold formula — not derived there (linear + nonlinear FMR/waveguide
  applications only); use d'Aquino 2026 Eqs. (20)/(24) for thresholds.
