# Function Origins — `spin-waves/index.html`

This document classifies every function in [spin-waves/index.html](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html) by its mathematical and engineering origin after the full **d'Aquino *et al.* (2026) Normal Modes Model (NMM)** overhaul (`arXiv:2601.07798v2` / PRB `10.1103/k77w-mslr`, `eigen_27mT_35modes.npz` pre-extracted tensors for the $3\,\mu\text{m}\times 50\,\text{nm}$ YIG disk at $B_{\text{DC}} = 27\,\text{mT}$, retaining strictly the $M = 29$ modes with $f_h > 1.6\,\text{GHz}$, `npz = 6..34`).

Every function falls into one of four categories:

1. **d'Aquino (2026) NMM** — Exact generalized eigenmodes (`Eq. (2)`), symplectic amplitudes $a_h = \bar{\omega}_h(\boldsymbol{\varphi}_h, \mathcal{B}_0 \delta\mathbf{m}_\perp)$, NMM equation of motion `Eq. (9)` ($G_1 + G_2 + G_3 + G_4$), RWA single-mode steady states (`Eqs. (33)–(34)`), 2-mode steady-state coupled analysis & Routh–Hurwitz stability (`Eqs. (55)–(57)`, Table I Scenarios 1 & 2), and Arnold instability tongues under a **strictly uniform external field** ($\nabla\mathbf{h}_{\text{ext}} = 0$).
2. **Brown (1963) FDT Thermal Bath** — Classical fluctuation–dissipation thermal field noise ($\Theta(T) = \frac{k_B T}{\mu_0 M_s^2 V}$, quadrature Langevin variances derived directly from $K^+_{hh}$ and $K^-_{hh}$).
3. **Experimental AWG Protocol & DSP** — 64 GS/s AWG `tukey(128, 1)[:64]` 1.0 ns cosine turn-on ramp, continuous dual-tone carrier phases $(\phi_1, \phi_2)$ / $\Delta\phi$, $3.0\,\mu\text{s}$ delayed Tone-2 sequencing (`Wait T1 SS` auto-release), and interpolated uniform-`SPEC_DT` 0–6.5 GHz Hann-windowed DFT.
4. **UI / Canvas Rendering / Guide / Telemetry Logging** — Slow-motion speed control (`0.25×`–`150×`), exact dynamic tensor profile $\delta\mathbf{m} = \mathbf{m} - \mathbf{m}_0$ rendering, interactive $6\times 6$ tensor inspector, 8-step button-spotlight coach guide with live press announcer, and structured telemetry + Supabase/IndexedDB evaluation logging.

> [!NOTE]
> **All custom / non-NMM code removed:** Spatial field gradients (`pumpGradProfile`, `mutScale`, `mutualPumpThreshold`), phenomenological matrices ($\mu_{ij}$, `KAPPA_MU`, `pair4`/`tri4`/`scat4`/`xmod4`, `nikCustom`), artificial NFS/threshold hacks (`chi4m`, `chi4nr`), `twoModeIsolate`, `projectOutMode`, and analytical Bessel proxies have been completely removed from `spin-waves/index.html`.

---

## 1. Tensor Unpacking & Thermal FDT Noise

| Function | Lines | Origin | Mathematical / Engineering Basis |
| :--- | :--- | :--- | :--- |
| `b64ToBuf` | [L363–L369](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L363-L369) | Data Unpacking | Decodes embedded Base64 precomputed 29-mode NMM tensors (`fGhzArr`, `Kp`, `Km`, `Q2x`, `Q2y`, `C3A`–`C3D`, `N4A`, `P4B`, `C4C`, `D4D`, `phiRe/Im`, `m0`, `actPix`) and builds `m0xFull`, `m0yFull`, `m0zFull` for `eigen_27mT_35modes.npz` (`npz = 6..34`, $f \in [1.7033, 2.0689]\,\text{GHz}$). |
| `thermalTheta` | [L499](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L499) | **Brown FDT (1963)** | Dimensionless thermal energy ratio $\Theta(T) = \frac{k_B T}{\mu_0 M_s^2 V}$ for the $3\,\mu\text{m}\times 50\,\text{nm}$ YIG disk ($\Theta(300\,\text{K}) = 4.7543\times 10^{-7}$). |
| `thermalRmsAmp` | [L500](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L500) | **Brown FDT (1963)** | Equipartition RMS modal amplitude $\sqrt{\langle |a_h|^2\rangle_{\text{eq}}} = \sqrt{\Theta(T)}$ ($6.895\times 10^{-4}$ at $300\,\text{K}$). |
| `reseedRng` | [L501](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L501) | PRNG Control | Seeds the deterministic 32-bit `Mulberry32` generator for reproducible shot-to-shot stochastic realizations. |
| `nextRand01` | [L502–L507](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L502-L507) | PRNG (`Mulberry32`) | Fast uniform $[0, 1)$ deviate generator used by Box–Muller Gaussian thermal sampling. |
| `ensureThermalSigmas` | [L510–L520](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L510-L520) | **Brown FDT + d'Aquino NMM** | Caches diagonal Langevin noise standard deviations $\sigma_{R,m} = \sqrt{\alpha\Theta(K^+_{mm} - K^-_{mm})\,h}$ and $\sigma_{I,m} = \sqrt{\alpha\Theta(K^+_{mm} + K^-_{mm})\,h}$ matching linear relaxation $\dot{a}_m = -i\bar\omega_m a_m - \alpha(K^+_{mm} a_m + K^-_{mm} a_m^*)$ so $\langle|a_m|^2\rangle_{\text{eq}} \equiv \Theta$. |

---

## 2. Mode Toggle, Experimental AWG Pulse Timing & Sampling

| Function | Lines | Origin | Mathematical / Engineering Basis |
| :--- | :--- | :--- | :--- |
| `isRwaFast` | [L523–L525](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L523-L525) | Solver Mode Selector | Returns `true` when only the secular RWA core (`g1`, `t2A`, `t4A`) is active, bypassing off-diagonal/non-RWA/3-wave loops for fast integration. |
| `setNmmMode` | [L526–L532](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L526-L532) | Solver Mode Selector | Toggles between **RWA Fast** (`g1+t2A+t4A`) and **Full NMM** (`g1+t2A+t2B+t2C+g3+t4A+t4B+t4C+t4D`). |
| `simTimeNs`, `simTimeUs` | [L533–L534](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L533-L534) | Unit Conversion | Converts dimensionless Larmor time $\tau = \gamma\mu_0 M_s t$ to nanoseconds and microseconds via `FREQ_TO_RAD`. |
| `nsToTau`, `usToTau` | [L535–L536](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L535-L536) | Unit Conversion | Converts nanoseconds and microseconds to dimensionless time $\tau$. |
| `tukeyRamp` | [L538–L544](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L538-L544) | **AWG Protocol (`tukey(128,1)[:64]`)** | Exact $1.0\,\text{ns}$ ($64$ samples at $64\,\text{GS/s}$) cosine taper $\tfrac{1}{2}\!\left(1 - \cos\frac{\pi\Delta\tau}{\tau_{\text{ramp}}}\right)$ matching `scipy.signal.windows.tukey(128, 1)[:64]`. |
| `sampleProbeSignal` | [L551–L563](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L551-L563) | **d'Aquino NMM** Probe Observable | Evaluates uniform-averaged $\langle \delta m_y\rangle = \frac{1}{N_{\text{act}}}\sum_{m,p} (a_{R,m}\varphi^y_{R,m,p} - a_{I,m}\varphi^y_{I,m,p})$ for the inductive probe readout and DFT. |
| `recordSpecSamples` | [L565–L583](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L565-L583) | Standard DSP | Linearly interpolates `sampleProbeSignal()` onto an exact uniform `SPEC_DT` grid ($2.21\,\tau$, Nyquist $13.38\,\text{GHz}$) inside every substep, eliminating step-size aliasing. |
| `resetModeHistory` | [L594–L600](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L594-L600) | Diagnostics Reset | Clears rolling mode history, `visualPeakDyn`, Tone-2 release flags, and DFT ring buffer. |
| `recordModeSample` | [L601–L608](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L601-L608) | Diagnostics Logging | Stores $|a_m|^2 = a_{R,m}^2 + a_{I,m}^2$ into the circular buffer. |
| `cycleAvgModeAmp`, `modeRmsAmp` | [L609–L617](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L609-L617) | Diagnostics Smoothing | Period-averaged / recent-window RMS modal amplitude $\sqrt{\langle |a_m|^2\rangle}$. |
| `reconstructSpatialSpins` | [L620–L635](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L620-L635) | **d'Aquino (2026) NMM** | Reconstructs normalized spatial spin field $\mathbf{m}(\mathbf{r}) = \frac{\mathbf{m}_0(\mathbf{r}) + \operatorname{Re}\sum_h a_h \boldsymbol{\varphi}_h(\mathbf{r})}{|\cdots|}$ on the disk ($|\mathbf{m}| \equiv 1$). |

---

## 3. Analytical Thresholds, $N_{hk}$ Shifts & Coupled Steady States

| Function | Lines | Origin | Mathematical / Engineering Basis |
| :--- | :--- | :--- | :--- |
| `modeFreqGHz` | [L638](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L638) | **d'Aquino Eq. (2)** | Exact generalized eigenfrequency $f_h$ from `eigen_27mT_35modes.npz` ($1.7033$–$2.0689\,\text{GHz}$). |
| `modeGammaGHz`, `modeGamma` | [L639–L640](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L639-L640) | **d'Aquino Eq. (9)** | Exact linear damping rate $\Gamma_h = \alpha K^+_{hh}$. |
| `vParamGHz_mT` | [L641–L643](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L641-L643) | **d'Aquino Eq. (17)/(32)** | Uniform parallel-pumping efficiency $v_h = \tfrac{1}{2}|Q^{(x)}_{2,hh}|$. |
| `parametricThresholdMT`, `parametricThreshold` | [L644–L649](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L644-L649) | **d'Aquino Eq. (32)** | Minimum parallel-pumping threshold $h_{\text{th},h} = \Gamma_h / v_h$. |
| `mutualNfsCoeffGHz`, `mutualNfsCoeff`, `selfNfsCoeff` | [L650–L652](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L650-L652) | **d'Aquino Eq. (29)/(52)** | Exact symmetric 4-wave secular Kerr/NFS tensor entries $N_{hk} = N_{kh}$ (`N4A`). |
| `modeNfsShiftGHz`, `modeMutualNfsShiftGHz` | [L653–L669](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L653-L669) | **d'Aquino Eq. (52)** | Total ($\sum_j N_{hj}|a_j|^2$) and mutual ($\sum_{j\neq h} N_{hj}|a_j|^2$) nonlinear frequency shifts. |
| `arnoldThresholdMT`, `arnoldThreshold` | [L670–L678](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L670-L678) | **d'Aquino Eq. (31)** | Detuned Arnold instability tongue $h_{\text{crit}}(\omega_p) = h_{\text{th},h}\sqrt{1 + (\Delta/\Gamma_h)^2}$. |
| `steadyStateAmp` | [L681–L694](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L681-L694) | **d'Aquino Eqs. (33)–(34)** | Analytical single-mode PWM saturated amplitude $|a_h|^2 = \frac{\Delta + \operatorname{sgn}(N_{hh})\sqrt{(v_h h_p)^2 - \Gamma_h^2}}{N_{hh}}$. |
| `routhHurwitzStable` | [L696–L699](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L696-L699) | **d'Aquino Sec. IV / Eq. (56)** | Quartic Routh–Hurwitz criterion ($c_1..c_4 > 0$, $c_1 c_2 > c_3$, $(c_1 c_2 - c_3)c_3 > c_1^2 c_4$) for linearized 2-mode phase-amplitude Jacobian stability. |
| `checkCoupledStab` | [L700–L713](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L700-L713) | **d'Aquino Eqs. (55)–(56)** | Evaluates linearized 2-mode coexistence Jacobian characteristic polynomial coefficients and tests `routhHurwitzStable`. |
| `analyzeTwoMode`, `steadyState2Mode` | [L714–L798](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L714-L798) | **d'Aquino Sec. IV & Table I** | Computes uncoupled vs $c_{++}/c_{-+}/c_{+-}/c_{--}$ coupled fixed points, $S_{\text{det}} = \operatorname{sgn}(N_{11} N_{22})(N_{11}N_{22} - N_{12}^2)$, Scenario 1 coexistence vs Scenario 2 ($S_{\text{det}} < 0$) solo bistability, and solo mutual-detuning lockout. |
| `closestArnoldMode` | [L799–L806](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L799-L806) | Diagnostic Helper | Returns mode $m$ minimizing $|f_p - 2f_m|$. |

---

## 4. Exact 29-Mode NMM Integration & AWG Sequencing

| Function | Lines | Origin | Mathematical / Engineering Basis |
| :--- | :--- | :--- | :--- |
| `symPairIdx` | [L808–L810](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L808-L810) | Tensor Indexing | Packed symmetric pair offset for sparse 3-wave/4-wave tensors (`C3*`, `C4C`, `D4D`). |
| `tone1EnvAt`, `tone2EnvAt` | [L812–L826](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L812-L826) | **AWG Protocol (`tukey(128,1)[:64]`)** | Applies 1.0 ns Tukey turn-on from $\tau = 0$ (Tone 1) and after `tone2DelayNs` or `tone2ReleaseT` under `Wait T1 SS` (Tone 2). |
| `fieldComponentsAt` | [L828–L838](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L828-L838) | **Uniform Field + Carrier Phase** | Uniform AWG RF field $h_{\text{ext}}(\tau) = e_1(\tau)h_{p1}\sin(\omega_1\tau + \phi_1) + e_2(\tau)h_{p2}\sin(\omega_2\tau + \phi_2)$ along $\hat{\mathbf{x}}$ (parallel) or $\hat{\mathbf{y}}$ (transverse). |
| `evalNmmRhs` | [L840–L1009](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L840-L1009) | **d'Aquino (2026) Eq. (9)** | Exact 29-mode NMM RHS $\dot{a}_h = (1+\alpha^2)^{-1}[G_1 + G_2(h_x,h_y) + G_3(a,a^*) + G_4(a,a^*)]$ covering linear $\bar{\omega}_h$, $(K^+, K^-)$, 1M direct $(b_x, b_y)$, 2M parallel/transverse $(P_2, Q_2)$, 3-wave $(C_{3A\text{--}D})$, and 4-wave $(N_{4A}, P_{4B}, C_{4C}, D_{4D})$. |
| `stepLlg` | [L1011–L1070](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1011-L1070) | **Heun (RK2) SDE + Wait T1 SS** | Advances coupled complex amplitudes $(a_{R,m}, a_{I,m})$ via Heun's method + Stratonovich/Itô quadrature Brown FDT thermal kicks $\mathcal{N}(0,\sigma_{R,m}^2), \mathcal{N}(0,\sigma_{I,m}^2)$, samples DFT at exact `SPEC_DT`, and triggers physical `tone2WaitSS` release once Tone 1 saturates. |
| `seedSpins`, `exciteEigenmode` | [L1072–L1110](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1072-L1110) | Initialization | Seeds zero (`0`), deterministic, or Gaussian thermal initial amplitudes $\mathcal{N}(0, \tfrac{1}{2}\Theta)$ cleanly via $a_m$ and updates spatial spins; clears preset highlights when selecting an eigenmode. |
| `pluckAt` | [L1112–L1151](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1112-L1151) | **d'Aquino Eq. (S12) + UI** | Projects a spatial Gaussian pointer kick $(\delta m_y, \delta m_z)$ directly onto all 29 modes via exact symplectic inner product $a_h = \bar{\omega}_h(\boldsymbol{\varphi}_h, \mathcal{B}_0\delta\mathbf{m}_\perp)$ and announces excited dominant modes. |

---

## 5. Spectrum, Tensor Inspector, Arnold Diagram & Dynamic Tensor Rendering

| Function | Lines | Origin | Mathematical / Engineering Basis |
| :--- | :--- | :--- | :--- |
| `resize` | [L1153–L1173](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1153-L1173) | Canvas Plumbing | Sizes full-viewport DPR-scaled field canvas and mobile drawer defaults. |
| `buildModeBars` | [L1175–L1191](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1175-L1191) | UI Telemetry | Populates clickable DOM bars for all 29 eigenmodes (`m0..m28`) with frequency, bulk/edge tags, and plain-English press announcements. |
| `updateSpectrum` | [L1193–L1257](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1193-L1257) | Standard DSP + Rendering | 128-point Hann-windowed DFT of uniform-`SPEC_DT` probe series over $0$–$6.5\,\text{GHz}$ with $1\text{G..}6\text{G}$ grid ticks, $f_p/2$, $f_p$, $3f_p/2$ harmonic markers, and `lastSpecPeakGHz` tracking. |
| `updateEigenmodeReadout` | [L1259–L1289](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1259-L1289) | UI Telemetry | Updates HUD dominant tensor mode badge (`TENSOR φ_m` or `TENSOR Σaφ`), modal bar widths, $h/h_{\text{th}}$, $\Delta f_{\text{tot}}$, and $N_{mm}$. |
| `buildNikMatrix`, `syncNikConfigUi`, `updateNikMatrix` | [L1300–L1418](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1300-L1418) | Read-Only Tensor Inspector | Renders interactive $6\times 6$ bulk mode inspector (`m10, m13, m16, m20, m23, m27`) for exact $N_{hk}$, $|P_{hk}|$, $K^+_{hk}$, and Table I $S_{\text{det}}$ scenarios with one-click `(i, k)` pumping and press explanations. |
| `drawArnoldTongues` | [L1426–L1551](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1426-L1551) | **d'Aquino Diagrams** | Plots Arnold instability tongues (`Eq. (31)`), nonlinear fold lines (`Eq. (33)` + mutual NFS shift), steady-state $|a|_{\text{ss}}(f_p)$ curves, and Tone 1 / Tone 2 operating markers. |
| `drawHInputView` | [L1555–L1604](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1555-L1604) | UI Waveform View | Draws Tone 1 + delayed Tone 2 $1.0\,\text{ns}$ Tukey turn-on profiles, $\Delta\phi$ phase tag, `SWAP`/`WAIT` status, and live $\tau$ playhead. |
| `selectArnoldAt` | [L1606–L1631](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1606-L1631) | Interactive UI | Maps pointer click/drag on Arnold canvas to `(driveFreq, driveAmp)` or `(driveFreq2, driveAmp2)`. |
| `drawScene` | [L1633–L1708](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1633-L1708) | **Exact Tensor Profile Rendering** | Renders strictly dynamic tensor perturbations $(\delta m_y, \delta m_z) = (m_y - m_{0y}, m_z - m_{0z})$ from `NMM_DATA.phi50` clipped cleanly to the circular YIG disk, smooth auto-scaled contrast (`visualPeakDyn`, zero clipping), and cell-bounded dynamic spin vectors. |

---

## 6. 8-Step Coach Guide, Slow Playback & Telemetry Database Logging

| Function | Lines | Origin | Mathematical / Engineering Basis |
| :--- | :--- | :--- | :--- |
| `isCompactViewport`, `syncCtrlHint`, `makeDraggablePanel`, `syncUi` | [L1715–L1816](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1715-L1816) | Responsive UI & State Sync | Mobile drawer detection, pointer-drag panel positioning, and synchronization of controls, speed pills (`0.25×`–`150×`), and presets. |
| `clearCoachSpotlights`, `applyCoachSpotlights` | [L1916–L1936](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1916-L1936) | Interactive Coach Guide | Applies pulsing `.tutSpotlight` (`PRESS`/`SLIDE`) glow to each step's target controls and auto-scrolls newly revealed sliders into view. |
| `swLog`, `captureLogState`, `updateLogUiBadges`, `recordLogEvent`, `recordLogError`, `recordTelemetrySample`, `buildEvaluationReport` | [L1955–L2110](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1955-L2110) | **Evaluation Telemetry** | Captures full session timeline: coalesced UI/slider actions, $1.5\,\text{s}$ physics trajectory samples, FPS stats, NaN/runtime errors, and live analytical `d'Aquino Eq. (33) / Eqs. (55)–(57)` vs simulated mode comparisons. |
| `openLogIdb`, `saveLogToIdb`, `uploadSwLog`, `handleSaveLogDb`, `exportLogJsonFile` | [L2112–L2186](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L2112-L2186), [L2565–L2597](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L2565-L2597) | **Database & Local Persistence** | Persists structured JSON evaluation reports (`spin-waves-evaluation-v1`) to Supabase (`spacetime_feedback_runs`, `app_version: 'spin-waves-nmm-29m-v2'`), IndexedDB (`spacetime-spinwaves-logs`), and downloadable JSON files (with `pagehide`/`visibilitychange` `keepalive` auto-upload). |
| `announcePress`, `renderSwTut`, `openSwTut`, `closeSwTut`, `notifySwTutorial`, `runSwTutDemo` | [L2188–L2278](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L2188-L2278) | Interactive 8-Step Guide | Explains every pressed button/slider in plain English (`#pressToast` + `#swTutPress`), advances the 8-step guide (`0.25×` Kittel, `m13` rings, threshold, $\Delta\phi$, `COEX c++`, `3µs NONRECIP`, `SWAP` + $300\,\text{K}$, Full NMM), and logs + auto-uploads on completion. |
| `applyPreset`, `swapTonesAction`, `resetSimTime` | [L2295–L2449](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L2295-L2449) | Presets & Protocol Actions | Configures `param`, `sub`, `direct`, `coex` (`m13,m23`), `opp` (`m13,m27`), and `nonrecip` (`m13,m16`, 3.0 µs `Wait T1 SS`); swaps Tone 1 ↔ Tone 2 including phases $(\phi_1, \phi_2)$. |
| `fastForwardStepsPerFrame`, `fastForwardSteadyState`, `toggleFastForward` | [L2508–L2543](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L2508-L2543) | Acceleration Helpers | Advances `stepLlg` synchronously or per frame to steady state, preserving spontaneous thermal nucleation from `a = 0` when $T > 0$. |
| `SPEED_FRAME_STEPS`, `frame` | [L2635–L2708](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L2635-L2708) | Slow-Motion Loop + Telemetry | Main `requestAnimationFrame` loop defaulting to slow visual playback (`1× SLOW = 0.06 τ/frame`, `0.25× = 0.025 τ/frame`, up to `150×`), live Tone-2 wait/delay release logging, FPS tracking, periodic IndexedDB saves, DFT spectrum, and rendering. |
