# Spin-Waves Scientific Audit & Full-NMM Findings

**Repository / Branch:** `opencode-work` (baseline `b5f9208` -> Subtask 1 verified)
**Scope:** [spin-waves/index.html](file:///Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/index.html) (29-mode d'Aquino et al. 2026 NMM, `eigen_27mT_35modes.npz`, modes `6..34`, $f_h \in [1.7033, 2.0689]\,\text{GHz}$).
**Reference Documentation:** See [spin-waves/FUNCTION_ORIGINS.md](file:///Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/FUNCTION_ORIGINS.md) for the line-by-line ODE RHS derivation with exact page/section/equation citations to d'Aquino *et al.* (2026), Perna *et al.* (2022), and Brown (1963).

## 1. Confirmed Physics, Tensors & Solver Audit

### Units, Constants & Basis Normalization
- $B_{\text{DC}} = 27.0\,\text{mT}$, $\mu_0 M_s = 176.0\,\text{mT}$, $\omega_M = 3.09911309\times 10^{10}\,\text{rad/s}$, $f_M = 4.93239167\,\text{GHz}$, `FREQ_TO_RAD = 1 / F_M`, $\tau = \omega_M t$.
- $\Theta_{\text{per K}} = k_B / (\mu_0 M_s^2 V) = 1.58478\times 10^{-9}\,\text{K}^{-1}$, $\Theta(300\,\text{K}) = 4.7543\times 10^{-7}$.
- Direct verification of embedded `NMM_DATA` tensors against `eigen_27mT_35modes.npz` (`npz = 6..34`) & theory:
  - Eigenfrequencies: `max |fGhz/fM - wBar| = 2.98e-8` and `max |npz["omega"][sel]/omegaM - wBar| = 1.46e-8` ($\bar\omega_h \in [0.345329, 0.419449]$).
  - Real-gauge condition (d'Aquino Eq. (12), p. 4): `max |Im(Km_mm)| = 8.23e-17`, `Km_mm / (wBar^2 * qDiag) = +1.0000001`, ensuring $K^+_{mm} - K^-_{mm} \ge 0.15916 > 0$ and positive-definite Brown FDT diffusion variances.
  - Linear `g1` completeness & torque sign: Summing $2\operatorname{Re}[b_{yL,h}\langle\boldsymbol{\varphi}_h\rangle] = [5.2\times 10^{-8}, -4.7\times 10^{-9}, -0.9053]$ and $2\operatorname{Re}[b_{yA,h}\langle\boldsymbol{\varphi}_h\rangle] = [6.1\times 10^{-6}, +0.8936, -1.7\times 10^{-9}]$. Note: `Perna (2022) Table 1` (p. 3) omits a minus sign in $b_h$ vs `Perna Eq. (43)` (p. 11, $-\mathbf{m}\times(\mathbf{m}\times\mathbf{f}^d)$); `spin-waves/index.html` uses the exact LLG sign in `bxA/byA`.
  - Inhomogeneous equilibrium $\mathbf{m}_0(\mathbf{r})$ in `g2` (`P2x/Q2x/P2y/Q2y`): Edge mode `m0` (`edgeFrac = 0.990`) gives `|P2x[0,0]|/(wBar*norm2) = 0.9608` and `|Q2x[0,0]|/(wBar*qDiag) = 0.9624` (vs `0.9970`/`0.9981` for bulk `m13`), proving exact projection per d'Aquino Eq. (10) / Perna Eq. (44) rather than the $\mathbf{m}_0\approx\mathbf{e}_x$ approximation of d'Aquino Eq. (13) (`3.8%` edge threshold difference, `<0.2%` bulk).
  - Symmetric NFS `4A` (`N4A`): `max |N4A[h,k]/wBar[h] - N4A[k,h]/wBar[k]| = 1.53e-5` (float32), confirming $N_{hk}=N_{kh}$ (d'Aquino Eq. (41), p. 9). `C4C[h,h,1] === 0` prevents double-counting secular $|a_h|^2 a_h$.
- Legacy `tone2DelayNs` stores **microseconds** (`usToTau(tone2DelayNs)`).

## 2. Full-NMM vs RWA Fast Nonreciprocity Measurements (Confirmed)

Tested `nonrecip` pair `(m13 = 1.978 GHz, m16 = 1.999 GHz)` at `alpha = 0.0010`, `delay = 3.0 us`, detunings $+0.05 r_1$ (`m13`) and $-0.20 r_2$ (`m16`):

### RWA Fast (`setNmmMode(false)`, 0 K, seed 42, nominal `2.4x / 2.0x h_th`)
- `T1->T2`: `a13 = 0.0465 (2us), 0.0465 (3us), 0.0418 (5us), 0.0355 (8us), 0.0360 (12us)`, `a16 = 0.00000`.
- `T2->T1`: `a16 = 0.0580 (2us), 0.0581 (3us), 0.0568 (5us), 0.0659 (8us), 0.0642 (12us)`, `a13 <= 0.00005`.
- **Conclusion (Confirmed):** RWA Fast produces persistent 2-mode nonreciprocity to `12 us`.

### Full NMM (`setNmmMode(true)`, all 29 modes, `G1+G2+G3+4ABCD`, `dt = 0.65 tau`)
1. **Nominal high drive `(2.4x, 2.0x h_th)` (mode-competition regime):**
   - **0 K (seed 42):**
     - `T1->T2`: `m13` dominates at `2 us` (`0.0635`) and `3 us` (`0.0623`), then competing **`m14`** steals dominance (`0.0659` at `5 us`, `0.0663` at `8 us`, `0.0634` at `12 us`), suppressing `m13` to `0.0146–0.0238` (`m16 <= 0.0093`).
     - `T2->T1`: dense neighbor **`m15`** (`df = 0.8 MHz`) dominates at `2–3 us` (`0.0297–0.0298` vs `m16 = 0.0140–0.0141`), then **`m10`** takes over at `5–12 us` (`0.0451–0.0490`), suppressing `m16` to `0.0017–0.0061`.
   - **Single-Term Ablation at 6 us (0 K):** `4B`, `2B`, `2C` individually preserve `(m13, m16)`; `4C` elevates `m15` (`0.0350`); **`G3`** and **`4D`** each independently trigger **`m10`** (`0.053–0.055`) takeover in `T2->T1`; `FULL_ALL` combines to produce `m14` (Fwd, `0.0628`) and `m10` (Swap, `0.0572`).
2. **Moderate / safe drive `(1.35x, 1.35x h_th)` and `(1.60x, 1.50x h_th)` (Full-NMM nonreciprocity SURVIVES to 12 us!):**
   - **0 K (`seed=42`), `(1.35x, 1.35x)`:**
     - `T1->T2`: `a13 = 0.03749 (2us), 0.03745 (3us), 0.03901 (5us), 0.03586 (8us), 0.03326 (12us)`, `a16 = 0.00000`, `dom = m13` across all `2–12 us`.
     - `T2->T1`: `a16 = 0.05603 (2us), 0.05518 (3us), 0.04697 (5us), 0.04126 (8us), 0.04767 (12us)`, `a13 <= 0.00486`, `dom = m16` across all `2–12 us`.
   - **0 K (`seed=42`), `(1.60x, 1.50x)`:**
     - `T1->T2`: `a13 = 0.05243 (2us), 0.05188 (3us), 0.04746 (5us), 0.04298 (8us), 0.04408 (12us)`, `a16 <= 0.00003`, `dom = m13` across all `2–12 us`.
     - `T2->T1`: `a16 = 0.05517 (2us), 0.05524 (3us), 0.05507 (5us), 0.04465 (8us), 0.04457 (12us)`, `a13 <= 0.00373`, `dom = m16` across all `2–12 us`.
   - **300 K multi-seed (`(1.35x, 1.35x)`, `t = 12 us`):** Pairwise `(m13, m16)` nonreciprocity holds cleanly across all 4 independent thermal seeds:
     - Seed `42`: `FWD dom=m13 (0.03286, a16=0.00302)` | `SWP dom=m16 (0.05414, a13=0.00831)`
     - Seed `20260101`: `FWD dom=m13 (0.02980, a16=0.00356)` | `SWP dom=m16 (0.05385, a13=0.00902)`
     - Seed `777`: `FWD dom=m13 (0.03237, a16=0.00068)` | `SWP dom=m16 (0.05430, a13=0.00769)`
     - Seed `9991`: `FWD dom=m13 (0.03073, a16=0.00346)` | `SWP dom=m16 (0.05402, a13=0.00723)`
- **Scientific Conclusion:** Full-NMM `(m13, m16)` arrival-order nonreciprocity is real and thermally robust at `1.35x–1.60x h_th` out to `12 us` (both `0 K` and `300 K` across 4 seeds), whereas raising the preset to `2.4x h_th` crosses a secondary multi-mode instability threshold (`G3 + 4D + 4C`) where `m14/m10` steal dominance after `~3 us`.

### Timestep Convergence (`Full NMM`, 0 K, `dt` halving)
- **`0.6 us` (`2.4x`):** `a13(dt=0.65) = 0.074733`, `a13(dt=0.325) = 0.074634` (`0.13%` rel diff).
- **`3.0 us` (`1.35x`, clean nonreciprocal window):**
  - `dt = 0.65`: `a13 = 0.037448`, `a14 = 0.006691`, `a11 = 0.009845`
  - `dt = 0.325`: `a13 = 0.037303`, `a14 = 0.006670`, `a11 = 0.009503` (`0.39%` vs `0.65`)
  - `dt = 0.1625`: `a13 = 0.037294`, `a14 = 0.006671`, `a11 = 0.009478` (`0.024%` vs `0.325` — `16.2x` error reduction on halving, matching $O(h^4)$ RK4 scaling).
- **`3.0 us` (`2.4x`, onset of `m14` multi-mode competition):**
  - `dt = 0.65`: `a13 = 0.06226` | `dt = 0.325`: `a13 = 0.07077` (`12.02%` rel diff due to phase/growth-rate sensitivity right as `m14` destabilizes). Note: testing `2.4x` with default `damping=0.0005` left active before rescaling `h_th` produces `NaN` at `dt=0.65`; always set `damping` prior to computing `parametricThresholdMT`.

## 3. Remediation Status
1. **[DONE — Subtask 1] State Continuity, Continuous Carrier Phase & Control Separation:** Implemented exact piecewise-linear carrier phase integral (`evalRampVal`, `integrateRampFreq`, `channelStateAt`, `advanceCarrierPhases`), configurable smooth parameter ramps (`paramRampNs`), `Live Smooth` vs `Staged` + `Apply Now`, separation of `Restart Sequence` (`t=0` rewind preserving spins) vs `Reset Simulation`, live state-preserving swap (`#btnLiveSwapTones`), and removal of the blind `3.0 us` `Wait T1 SS` force-release. Verified by `Test AB` (`29/29 PASS` in `tests/physics.test.mjs`).
2. **[NEXT — Subtask 2] Multi-Window Convergence & Genuine Continuous CW Hysteresis Sweep.**
3. **[PENDING — Subtask 3] Configurable Welch PSD Spectrum, Physical-Time Arnold Averaging & Extrema-Preserving Fast-Mode Time Plots.**
4. **[PENDING — Subtask 4] Synchronized Clamped Numeric Inputs, Dashboard Responsive QA (Desktop & Phone) & Final Verification.**
