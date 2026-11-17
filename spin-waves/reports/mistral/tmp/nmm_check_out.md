# nmm_check_out.md — NMM solver (index.html) vs GAMMA_EFF_MATH.md: internal verification

**BOSS RULE respected:** internal verification only; the only mode identifiers below are
internal solver indices of the embedded IIFE (`mIdx` into the 29-mode NMM table). Nothing
here is paper-facing. No data files or `index.html` were modified.

---

## (a) The solver's exact definitions (verbatim from `index.html`, units as coded)

All rates the threshold functions use are in **GHz (cycles)**, fields in **mT**. `F_M = 4.932391670735886` GHz, `MU0_MS_MT = 176.0` mT, `damping = 0.0005` (line 1092).

| Quantity | Solver definition (line) | Units |
|---|---|---|
| $\Gamma_m$ | `modeGammaGHz(m) = max(1e-7, damping * Kp[(m*M+m)*2] * F_M)` (l. 2170) | GHz (also modeGamma = Γ/F_M per τ) |
| $v_m$ | `vParamGHz_mT(m) = 0.5 * hypot(Q2x[(m*M+m)*2], Q2x[..+1]) * F_M / MU0_MS_MT` (l. 2172) | GHz/mT |
| $\delta$ | `delta = 0.5*(fPump - 2*(f[m] + fShift))` (l. 2254, 2266) | GHz (solver $\delta = +\tfrac12(f_p-2f_m)$; doc $\delta = 2\pi(f_m-f_p/2)$ → sign/2π-convention, see (d)) |
| bare threshold | `arnoldThresholdMT(m,fp) = (Γ/v)·hypot(1, δ/Γ) = hypot(Γ,δ)/v` (l. 2250) | mT |
| 4A NFS shift | `modeNfsShiftGHz(m) = Σ_j N4A[m,j]·F_M·|a_j|²` (self incl.); `modeMutualNfsShiftGHz` excludes self (l. 2184) | GHz |
| 4B cross-pump | `crossPump4BGHz(m, fp, φ, fRef)`: $\tilde W = \sum_{k\neq m} P_{4B}[m,k]\,F_M\,a_k^2$, resonance-filtered per $k$ by `bw/hypot(bw, δ_k)`, `bw = max(8 MHz, 3Γ_m)` (l. 2204) | GHz |
| phase threshold | `instantaneousArnoldThresholdMT(m, fp, φ, fRef) = max(0, sqrt(req² − wPerp²) − wPar)/v`, `req = hypot(Γ, δ̃)` with $\tilde\delta = \tfrac12(f_p - 2(f_m + \Delta_{4A}))$ (l. 2260) | mT |
| 2A/2B sideband | rotating sideband quadrant $u = -i\,e^{i\theta_p} Q_{2x}/|Q_{2x}|$ multiplies $a_h^*$; `wPar = Re(W̃ ū)`, `wPerp = |Im(W̃ ū)|` | GHz |
| steady state (1 mode) | `steadyStateAmp`: $|a|^2 = (\delta + \mathrm{sgn}(n_{kk})\,r)/n_{kk}$, $r=\sqrt{(vh)^2-\Gamma^2}$ (l. 2291) | a.u. |
| 2-mode SS/stability | `analyzeTwoMode` (Eqs. 55–60 + Routh–Hurwitz, l. 2309); delayed path selected by `tone2WaitSS || tone2DelayNs > 0` | — |

The ½ in $v_m$ is exactly the $\sin\theta_p \to \tfrac{1}{2i}(e^{i\theta_p}-e^{-i\theta_p})$ RWA half: the raw RHS coupling is `hx·Q2x·a*` (rate $|Q_{2x}|F_M h_{mT}/176$), whose slow rotating component carries the ½. Hence **`v·h` is the doc's $Vh$ in Eq. (1)`s convention** — no extra factor 2.

Globals stubbed to run the functions headless (all values verbatim from the IIFE):
`damping = 0.0005` (l. 1092), `simTime = 0`, `dualTone = true` (forced), `tone2WaitSS = true` (forced; app default `false`, l. 1082), `tone2DelayNs = 3.0` (l. 1081), `carrierPhase1/2 = 0`, `ramp1/ramp2.active = false`, `driveFreq/2`, `driveAmp`, `drivePhase1/2` = verbatim defaults, `terms` all toggles forced `true` (app default l. 1288: `g1,t2A,t4A` true). `channelStateAt` copied verbatim (l. 1336); its ramp branch (`evalRampVal`, `integrateRampFreq`) is dead with `active:false` and stubbed to throw if ever reached. `aRe/aIm` = `Float64Array(M)` verbatim (l. 1114). The 4-magnon tensors in `NMM_DATA` are named `N4A, P4B, C4C, D4D` ("P4A" does not exist under that name; checked l. 981–1000).

Verbatim-copy audit: every code line of all copied functions (`b64ToBuf`, `modeGammaGHz`, `vParamGHz_mT`, `modeNfsShiftGHz`/`mutualNfs`, `crossPump4BGHz`, `arnoldThresholdMT`, `instantaneousArnoldThresholdMT`, `steadyStateAmp`, `routhHurwitzStable`, `checkCoupledStab`, `analyzeTwoMode`, `steadyState2Mode`, `closestArnoldMode`, `channelStateAt`) was diffed line-by-line against machine-extracted blocks of `index.html` (in `tmp/verbatim/`): **all code lines identical**; only stripped trailing comments and 4 intentionally unused constructs (`pp/pm/mm` arrays, `syncCarrierPhaseToSimTime`, `symPairIdx`) differ.

## [1] Mode table (from `NMM_DATA`), all 29 modes

Full table in `nmm_check_stdout.txt`. Summary at `damping = 0.0005`:

- $f_m$ ∈ [1.703297, 2.068889] GHz (internal idx 0…28; `selIdx` = npz 6…34)
- $\Gamma_m$ ∈ [1.0412, 1.5478] **MHz** — i.e. $\Gamma/\alpha \in [2.08, 3.10]$ GHz $= (0.42\text{–}0.63)f_M$
- $v_m$ ∈ [0.009063, 0.015749] GHz/mT
- bare threshold $\Gamma/v$ ∈ [0.0973, 0.1269] mT

## [2] Single-mode Arnold threshold vs doc Eq. (2): **EXACT MATCH**

Eq.(2) is $h_{\rm th}(\delta) = \sqrt{\Gamma^2+\delta^2}/V$. Solver `arnoldThresholdMT` returns mT; the comparison quantity $\Gamma\,\mathrm{hypot}(1,\delta/\Gamma)$ is its **numerator in rate units (GHz)** — it must be divided by $v$ (or the solver value multiplied by $v$) to compare. Both forms checked for **all 21 internal indices with $f_m\in[1.9,2.1]$ GHz** × δ ∈ {0, 3, 6, 10} MHz (84 rows in stdout). Representative rows:

| idx | $f_m$ GHz | δ MHz | solver [mT] | Eq.(2)$/v$ [mT] | rel. dev |
|---|---|---|---|---|---|
| 8 | 1.944505 | 0 | 0.102268809 | 0.102268809 | 0 |
| 8 | 1.944505 | 10 | 0.759804869 | 0.759804869 | 8.8e-16 |
| 16 | 1.999151 | 3 | 0.223111920 | 0.223111920 | 3.0e-14 |
| 28 | 2.068889 | 10 | 0.938885073 | 0.938885073 | 2.1e-14 |

**Max relative deviation: 3.5e-14** (float32 storage of `Kp` — machine precision). Also `parametricThresholdMT` ≡ $\Gamma/v$ exactly. No factor-2 discrepancy: the ½ in $v_m$ is the RWA convention and already matches Eq. (2).

## [3] Phase-dependent two-mode threshold vs doc Eq. (6)/(7)

Experimental pumps $f^{p}_1 = 4.070$ GHz, $f^{p}_2 = 4.142$ GHz.
`closestArnoldMode` (matches $|f_p - 2f_m|$) selects the two-photon partners:
- pump 1 → internal idx 22, $f_m = 2.035551$ GHz, $|f_p - 2f_m| = 1.10$ MHz (>the pair partner idx 21 by margins >3 MHz)
- pump 2 → internal idx 28, $f_m = 2.068889$ GHz, $|f_p - 2f_m| = 4.22$ MHz

Parameters: $v_1h_1 = 1.3\,\Gamma_1$ → $h_1 = 0.1604$ mT (normalized 9.115e-4). Γ₁ = 1.2920 MHz, Γ₂ = 1.3205 MHz; v₂ = 0.010743 GHz/mT.
`steadyStateAmp(22, 4.070, h1, 0, pwm) = 6.6192e-3` ($|a|^2 = 4.3814\times10^{-5}$) — nonzero because the self-NFS ($n_{11}$ = 11.9 GHz per |a|²) regexp-resonances the branch; **note $v_1h_1 = 1.3Γ_1$ alone would be above the bare single-mode threshold** (PWM resonance requires |δ₁| < r₁; δ₁ = −0.5515 MHz, r₁ = 1.073 MHz), so this operating point is the sub-ignition edge, consistent with the doc's "seed below ITS threshold" only in the PWM sense.

4B cross-pump onto the target with only the seed excited:
- $W = |\tilde W| = 5.2499\times10^{-5}$ GHz (resonance-filtered; unfiltered total 2.351e-4 GHz; dominant driver = the seed, index 22)
- **$W/v_2 = 4.887\times10^{-3}$ mT**; $W/\Gamma_2 = 0.0398$ (≤ 4 % → deep linear regime, as assumed in Eq. (7))
- internal pump angle χ = atan2(Im, Re) = π/2; arg(Q2x_diag(target)) = π/2

16-phase scan (φ = k·2π/16) of `instantaneousArnoldThresholdMT(28, 4.142, φ, 4.142)` vs doc:
- **Eq. (6) exact** $h_{\rm th}(\Psi) = [\sqrt{\Gamma_2^2+\tilde\delta_2^2 - W^2\sin^2\Psi} - W\cos\Psi]/v_2$ with $\tilde\delta_2$ = δ₂ + 4A shift: **max deviation 2.1e-11 mT** (machine precision) after a single fitted constant offset φ₀ = 1.5708 rad (= π/2, exactly the $-iQ/|Q|$ sideband rotation — pure phase convention; residual RMS 1.5e-11 mT). Peak-to-valley = 9.773e-3 mT = 2W/v₂ exactly.
- **Eq. (7) linearized** (in W): max dev 5.8e-5 mT, RMS 3.6e-5 mT ≈ 0.6 % of the modulation depth — the expected O(W²/req²) ≈ 6e-4 truncation.
- **Crucial:** if the doc formula uses the *bare* δ₂ = 2.111 MHz (ignoring the seed's 4A mutual shift of the target, $\Delta_{4A} = n_{21}|a_1|^2 = 0.340$ MHz → $\tilde\delta_2$ = 1.771 MHz), the curve misfits by up to **2.6e-2 mT ≈ 13 % of h_th(0) = 0.2056 mT and 2.7× the modulation depth**. The doc's Eq. (7) must be written with the NFS-shifted detuning whenever the seed carries power.
- Phase-averaged threshold 0.205621 mT vs bare 0.205650 mT: unchanged to linear order — doc consequence 3 confirmed (residual = the O(W²) term).

**Verdict: Eq. (6)/(7) is exactly the solver's `instantaneousArnoldThresholdMT` structure; $h_{4B}/h_2 = W/v_2/h_2 = 3.5\%$ at this seed drive (in field units, comparable to the doc's fig3 bound $|h_{4B}/h_2|\approx0.005$ only within an order of magnitude — depends on seed amplitude).**

## [4] Delayed-tone protocol (tone2WaitSS / tone2DelayNs > 0) and the memory problem

Regime map from `analyzeTwoMode` (delayed path), seed at steady state, Δ = $f_2^{p}-f_1^{p}$, target = `closestArnoldMode(f_2^p)`:

| Δ MHz | h₂=0.10 | 0.14 | 0.20 | 0.30 | 0.50 | 0.80 | 1.50 mT |
|---|---|---|---|---|---|---|---|
| 40 | c++ | un | un | un | un | un | un |
| 44 | uh | uh | c++ | un | un | un | un |
| 48 | uh | uh | uh | c++ | un | un | un |
| 52 | uh | uh | uh | uh | un | un | un |
| 56 | uh | uh | uh | uh | uh | un | un |
| 60 | uh | uh | uh | uh | un | un | un |
| 64 | uh | uh | uh | uh | un | un | un |
| 68 | uh | un | un | un | un | un | un |
| 72 | uh | uh | uh | un | un | un | un |
| 76 | uh | uh | uh | uh | un | un | un |
| 80 | uh | uh | uh | uh | uh | un | un |
| 84 | uh | uh | uh | uh | uh | un | un |
| 88 | uh | uh | uh | uh | uh | uh | un |
| 100 | uh | uh | uh | uh | uh | uh | Q/osc |

Mechanics extracted from the code:
1. The delayed path is a **static steady-state selection rule** — it contains **no `t_gap`, no time integration, and no $W(t_{\rm gap})$ decay term**. `tone2DelayNs` only gates which branch is *selected*; nothing in `analyzeTwoMode` evolves in time.
2. Its only Δ-dependence is through $d_2^{\rm bare} = \tfrac12(f_2^p - 2f_{\rm tgt}) = \tfrac12(\Delta - \Delta_c)$ with $\Delta_c = 2f_{\rm tgt} - f_1^p$ (67.78 MHz for the pump-2 partner) plus the discrete `closestArnoldMode` switch at **Δ = 53.49 MHz** (fp₂ midpoint between the two highest two-photon lines; seed fixed at pump 1). The measured collapse near Δ = 74.5 MHz does **not** coincide with any switch line of the solver (Δ_c + shifts → engagement window sits near Δ ≈ 68–71 MHz at h₂ ≈ 0.14 mT).
3. Engagement condition along the delayed path: $|d_2 - n_{21}|a_1|^2| \lessgtr r_2(h_2)$, $r_2 = \sqrt{(v_2h_2)^2-\Gamma_2^2}$. So for fixed $h_2$ there is a threshold-like critical Δ*, and $h_{2,\min}(\Delta^*)\propto|\Delta-\Delta_c|/2$ — a **threshold curve, not a lifetime**. Regimes found: `uh` (seed-only SS, tone 2 disengaged), `un` (tone-2 self-oscillation), `c++` (S_det > 0 coupled), `Q/osc` (coupled branch unstable). `c+-`/`c-+` never selected in this configuration (S_det > 0 for all pairs probed).
4. **Answer to the task question:** *no regime boundary of the delayed path can reproduce a memory lifetime that decreases monotonically with Δ.* The solver has (i) no decaying variable, (ii) Γ independent of Δ, (iii) boundaries that move linearly in Δ with slope ½ rather than producing an exponential with Δ-dependent τ. What it *can* offer is a qualitative reinterpretation consistent with the doc's §4.3 burst-wise picture: crossing $h_{2,\min}(\Delta)$ at fixed $h_2$ switches the post-tone-1 attractor between seed-dominated (`uh`) and tone-2-dominated (`un`) states — i.e. phase-contrast collapse by attractor re-selection, with the collapse line movable by amplitude; the measured Δ-dependent τ_mem (PS >7 ns at Δ≤58, 1–2.5 ns at 74–78, <1 ns at ≥80 MHz) must then come from dynamics absent from this static rule (nonlinear ring-down of the over-threshold driven seed, or hybrid-mode damping near the pair crossing), as the doc already conjectures. With solver Γ values, Eq. (10)'s τ_mem = 1/(2Γ_k) = 387 ns (α = 5e-4) or 97 ns (α = 2e-3, idx 22) — far above measured 1–7 ns at any Δ, reinforcing that Eq. (10) is not the operative decay even at low Δ.

## (d) Factor-2 / convention findings

1. **Units of the Eq.(2) comparison:** $\Gamma\,\mathrm{hypot}(1,\delta/\Gamma)$ is a *rate* (GHz); the mT threshold is that **divided by V**. `arnoldThresholdMT` returns mT. Comparison then matches to 3.5e-14 (float32 floor). No factor-2 error.
2. **½ in `vParamGHz_mT` is the RWA half of the sinusoidal pump, not a convention factor of 2 opposite to the doc's** — solver and doc Eq. (2) use the same convention ($Vh$ = rotating pump element). Doc's caveat "their V is half ours" applies to textbooks, not to this solver.
3. **GHz (cycles) vs rad/s:** the solver works entirely in GHz-cycles ($\delta$, Γ, $v h$, W). The doc's Eq. (1) writes $\delta$ in rad/s ($2\pi(f_m-f_p/2)$). As long as all rates in one equation share the unit, thresholds are unaffected (all-in-cycles ≡ all-in-rad/s up to a global 2π); but the doc's §1 statement "Γ [Hz]" should say *cycles* (Hz, not rad/s) to be drop-in compatible with `arnoldThresholdMT`.
4. **Γ scale discrepancy (≈ factor 2), needs doc correction:** doc quotes $\Gamma = 8.2$–$9.8$ MHz for α = 0.002 (⇒ Γ/α = 4.1–4.9 GHz ≈ 0.83–1.0 f_M). Solver: Γ/α = 0.42–0.63 f_M ⇒ at α = 0.002 the solver gives 4.16–6.19 MHz — a factor ≈ 1.6–2.0 smaller. Either the doc's NMM-report number used a different linewidth convention (e.g. full resonance width vs rate) or a different α; resolve before printing Figs 1–4 theory overlays. With the solver's own Γ, h_th(0) ≈ 0.10–0.13 mT per the mode table.
5. **4A shift is not optional in Eq. (7):** with the seed excited the target's detuning is shifted by $n_{21}|a_1|^2$ (0.34 MHz here; self-clock $n_{22}$ = 6.40 GHz per |a|²). Using the bare detuning breaks the 16-phase fit by 13 % of h_th(0) — larger than the phase modulation itself (2.4 %). Doc Eq. (7) should read $h_{\rm th}(\Psi) = [\sqrt{\Gamma_2^2 + \tilde\delta_2^2 - W^2\sin^2\Psi} - W\cos\Psi]/v_2$ with $\tilde\delta_2 = \delta_2 + \Delta^{(4A)}_2$.
6. **2B/P2x term exists in the RHS but not in any threshold formula:** the diagonal $P_{2x}$ is 1.65–1.67× larger than $Q_{2x}$ (and parallel to it, both purely imaginary at these indices). It multiplies $a$ (not $a^*$) — a pump-phase-locked frequency shift, i.e. an extra detuning shift of order $\tfrac12 h_2\,|P_{2x}|F_M/\mu_0M_s \approx$ MHz-scale at h₂ ~ 0.1–0.2 mT — that neither `arnoldThresholdMT` nor `instantaneousArnoldThresholdMT` includes (only the 4A shift). If experiments show h_th(δ) curves horizontally shifted vs Eq. (2) beyond the 4A shift, this 2B term is the candidate.
7. **`closestArnoldMode` ↔ target assignment is discrete:** under a pump sweep the engaged target jumps at fp₂ = midpoints of $2f_m$ lines (for fp₁ = 4.070 GHz: Δ = 33.10, 38.21, 53.49 MHz within the experimental band). Any data interpretation invoking "the mode nearest each pump" should quote these switch lines, not the 74.5 MHz line.

## (e) Paths

- Data: `reports/mistral/tmp/nmm_data.json` (extracted from `index.html` line 963 via regex `const NMM_DATA=(\{.*?\});`)
- Script: `reports/mistral/tmp/nmm_check.mjs` (run: `/opt/homebrew/bin/node nmm_check.mjs`)
- Full stdout: `reports/mistral/tmp/nmm_check_stdout.txt` (248 lines; all mode tables, all 84 Eq.(2) rows, 16-phase threshold table, regime maps)
- Verbatim audit blocks: `reports/mistral/tmp/verbatim/*.txt` (machine-extracted from `index.html`, lines 965–1000, 1336–1350, 2169–2413)
- Reference doc: `reports/mistral/GAMMA_EFF_MATH.md`
