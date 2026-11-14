# Function Origins — `spin-waves/index.html`

This document classifies every function in [spin-waves/index.html](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html) by its mathematical and engineering origin after the full **d'Aquino *et al.* (2026) Normal Modes Model (NMM)** overhaul (`arXiv:2601.07798v2` / PRB `10.1103/k77w-mslr`, `eigen_27mT_35modes.npz` pre-extracted tensors for the $3\,\mu\text{m}\times 50\,\text{nm}$ YIG disk at $B_{\text{DC}} = 27\,\text{mT}$, retaining strictly the $M = 29$ modes with $f_h > 1.6\,\text{GHz}$, `npz = 6..34`).

Every function falls into one of four categories:

1. **d'Aquino (2026) NMM** — Exact generalized eigenmodes (`Eq. (2)`), symplectic amplitudes $a_h = \bar{\omega}_h(\boldsymbol{\varphi}_h, \mathcal{B}_0 \delta\mathbf{m}_\perp)$, NMM equation of motion `Eq. (9)` ($G_1 + G_2 + G_3 + G_4$), RWA single-mode steady states (`Eqs. (33)–(34)`), 2-mode steady-state coupled analysis & Routh–Hurwitz stability (`Eqs. (55)–(57)`, Table I Scenarios 1 & 2), physical 4-magnon cross-pumping (`4B` $P_{hk} a_h^* a_k^2$), and Arnold instability tongues under a **strictly uniform external field** ($\nabla\mathbf{h}_{\text{ext}} = 0$).
2. **Brown (1963) FDT Thermal Bath** — Classical fluctuation–dissipation thermal field noise ($\Theta(T) = \frac{k_B T}{\mu_0 M_s^2 V}$, quadrature Langevin variances derived directly from $K^+_{hh}$ and $K^-_{hh}$).
3. **Experimental AWG Protocol & DSP** — 64 GS/s AWG `tukey(128, 1)[:64]` 1.0 ns cosine turn-on ramp, continuous dual-tone carrier phases $(\phi_1, \phi_2)$ / $\Delta\phi$, $3.0\,\mu\text{s}$ delayed Tone-2 sequencing (`Wait T1 SS` auto-release), and interpolated uniform-`SPEC_DT` 4096-point radix-2 FFT (`3.42 MHz/bin`) with `1.5–2.3 GHz` mode-resolved view + `0–6.5 GHz` harmonic view.
4. **UI / Canvas Rendering / Verification / Telemetry Logging** — Slow-motion speed control (`0.25×`–`150×`), exact dynamic tensor profile $\delta\mathbf{m} = \mathbf{m} - \mathbf{m}_0$ rendering, smooth non-instantaneous time-plot Y-axis slewing, 29-mode unclipped Arnold diagram (`ALL 29` / `BULK`) with interactive ZOOM/PAN/AUTO viewport, interactive $6\times 6$ tensor inspector, 8-step button-spotlight coach guide, dismissible launch `"What's New"` modal (`spacetime-shell-v79`), interactive 8-stage Benchmark / Calibration / Verification Protocol + Live `4A` vs `4B` Probe + `BENCHMARK HANDOUT` viewer, and Supabase/IndexedDB evaluation logging with HTTP/RLS diagnostics.

> [!NOTE]
> **All custom / non-NMM code removed:** Spatial field gradients (`pumpGradProfile`, `mutScale`, `mutualPumpThreshold`), phenomenological matrices ($\mu_{ij}$, `KAPPA_MU`, `pair4`/`tri4`/`scat4`/`xmod4`, `nikCustom`), artificial NFS/threshold hacks (`chi4m`, `chi4nr`), `twoModeIsolate`, `projectOutMode`, and analytical Bessel proxies have been completely removed from `spin-waves/index.html`.

---

## 1. Tensor Unpacking & Thermal FDT Noise

| Function | Lines | Origin | Mathematical / Engineering Basis |
| :--- | :--- | :--- | :--- |
| `b64ToBuf` | [L470–L476](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L470-L476) | Data Unpacking | Decodes embedded Base64 precomputed 29-mode NMM tensors (`fGhzArr`, `Kp`, `Km`, `Q2x`, `Q2y`, `C3A`–`C3D`, `N4A`, `P4B`, `C4C`, `D4D`, `phiRe/Im`, `m0`, `actPix`) and builds `m0xFull`, `m0yFull`, `m0zFull` for `eigen_27mT_35modes.npz` (`npz = 6..34`, $f \in [1.7033, 2.0689]\,\text{GHz}$). |
| `thermalTheta` | [L608](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L608) | **Brown FDT (1963)** | Dimensionless thermal energy ratio $\Theta(T) = \frac{k_B T}{\mu_0 M_s^2 V}$ for the $3\,\mu\text{m}\times 50\,\text{nm}$ YIG disk ($\Theta(300\,\text{K}) = 4.7543\times 10^{-7}$). |
| `thermalRmsAmp` | [L609](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L609) | **Brown FDT (1963)** | Equipartition RMS modal amplitude $\sqrt{\langle |a_h|^2\rangle_{\text{eq}}} = \sqrt{\Theta(T)}$ ($6.895\times 10^{-4}$ at $300\,\text{K}$). |
| `reseedRng` | [L610](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L610) | PRNG Control | Seeds the deterministic 32-bit `Mulberry32` generator for reproducible shot-to-shot stochastic realizations. |
| `nextRand01` | [L611–L616](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L611-L616) | PRNG (`Mulberry32`) | Fast uniform $[0, 1)$ deviate generator used by Box–Muller Gaussian thermal sampling. |
| `ensureThermalSigmas` | [L619–L629](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L619-L629) | **Brown FDT + d'Aquino NMM** | Caches diagonal Langevin noise standard deviations $\sigma_{R,m} = \sqrt{\alpha\Theta(K^+_{mm} - K^-_{mm})\,h}$ and $\sigma_{I,m} = \sqrt{\alpha\Theta(K^+_{mm} + K^-_{mm})\,h}$ matching linear relaxation $\dot{a}_m = -i\bar\omega_m a_m - \alpha(K^+_{mm} a_m + K^-_{mm} a_m^*)$ so $\langle|a_m|^2\rangle_{\text{eq}} \equiv \Theta$. |

---

## 2. Mode Toggle, Experimental AWG Pulse Timing & High-Res DSP

| Function | Lines | Origin | Mathematical / Engineering Basis |
| :--- | :--- | :--- | :--- |
| `isRwaFast` | [L632–L634](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L632-L634) | Solver Mode Selector | Returns `true` when only the secular RWA core (`g1`, `t2A`, `t4A`) is active, bypassing off-diagonal/non-RWA/3-wave loops for fast integration. |
| `setNmmMode` | [L635–L639](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L635-L639) | Solver Mode Selector | Toggles between **RWA Fast** (`g1+t2A+t4A`) and **Full NMM** (`g1+t2A+t2B+t2C+g3+t4A+t4B+t4C+t4D`). |
| `simTimeNs`, `simTimeUs` | [L642–L643](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L642-L643) | Unit Conversion | Converts dimensionless Larmor time $\tau = \gamma\mu_0 M_s t$ to nanoseconds and microseconds via `FREQ_TO_RAD`. |
| `nsToTau`, `usToTau` | [L644–L645](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L644-L645) | Unit Conversion | Converts nanoseconds and microseconds to dimensionless time $\tau$. |
| `tukeyRamp` | [L647–L651](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L647-L651) | **AWG Protocol (`tukey(128,1)[:64]`)** | Exact $1.0\,\text{ns}$ ($64$ samples at $64\,\text{GS/s}$) cosine taper $\tfrac{1}{2}\!\left(1 - \cos\frac{\pi\Delta\tau}{\tau_{\text{ramp}}}\right)$ matching `scipy.signal.windows.tukey(128, 1)[:64]`. |
| `runSpecFft4096` | [L694–L722](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L694-L722) | Standard DSP (`Cooley–Tukey FFT`) | Iterative in-place 4096-point radix-2 FFT yielding `3.42 MHz` frequency bin spacing across `0–7.0 GHz`. |
| `sampleProbeSignal` | [L724–L732](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L724-L732) | **d'Aquino NMM** Probe Observable | Multi-pixel asymmetric rim + bulk inductive pickup (`probeCoeffRe/Im`) ensuring every mode $m_0\text{..}m_{28}$ has strong SNR. |
| `recordSpecSamples` | [L734–L756](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L734-L756) | Standard DSP | Linearly interpolates `sampleProbeSignal()` onto an exact uniform `SPEC_DT` grid ($2.21\,\tau$, $f_s \approx 14.01\,\text{GHz}$) inside every substep and tracks `probeCount` / `specSeq`. |
| `SPEED_FRAME_STEPS`, `quantizePlotTop`, `recordPlotSample`, `updatePlotWindow` | [L758–L835](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L758-L835) | Moving Time-Window History | Configures visual speed tiers (`0.25×`–`150×`), quantized Y-axis target stepping (`quantizePlotTop`), decimated ring buffer (`PLOT_HIST_CAP = 512`), and adaptive window span tracking. |
| `resetModeHistory`, `recordModeSample`, `cycleAvgModeAmp`, `modeRmsAmp` | [L837–L865](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L837-L865) | Diagnostics & Reset | Resets rolling RMS/time-window/4096-pt FFT buffers cleanly and computes modal amplitudes $|a_m|$ and window RMS values. |
| `reconstructSpatialSpins` | [L867–L882](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L867-L882) | **d'Aquino (2026) NMM** | Reconstructs normalized spatial spin field $\mathbf{m}(\mathbf{r}) = \frac{\mathbf{m}_0(\mathbf{r}) + \operatorname{Re}\sum_h a_h \boldsymbol{\varphi}_h(\mathbf{r})}{|\cdots|}$ on the disk ($|\mathbf{m}| \equiv 1$). |

---

## 3. Analytical Thresholds, $N_{hk}$ Shifts, `4B` Cross-Pumping & Coupled Steady States

| Function | Lines | Origin | Mathematical / Engineering Basis |
| :--- | :--- | :--- | :--- |
| `modeFreqGHz`, `modeGammaGHz`, `vParamGHz_mT`, `parametricThresholdMT` | [L885–L896](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L885-L896) | **d'Aquino Eqs. (2, 9, 17, 32)** | Exact generalized eigenfrequencies ($1.7033$–$2.0689\,\text{GHz}$), linear damping $\Gamma_h = \alpha K^+_{hh}$, parallel-pumping efficiency $v_h$, and distinct 29-mode thresholds $h_{\text{th},h} = \Gamma_h / v_h$. |
| `mutualNfsCoeffGHz`, `modeNfsShiftGHz`, `modeMutualNfsShiftGHz` | [L897–L918](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L897-L918) | **d'Aquino Eqs. (29, 52)** | Exact symmetric 4-wave secular Kerr/NFS tensor $N_{hk}$ (`N4A`), gated by `terms.t4A` so unchecking `4A N_hk` pins Arnold tongues horizontally. |
| `crossPump4BGHz`, `instantaneousArnoldThresholdMT`, `arnoldZoomAt`, `arnoldPanBy`, `arnoldResetAuto` | [spin-waves/index.html](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html) | **d'Aquino Eq. (9) `4B` Term + UI Viewport** | Evaluates coherent phase-resolved 4-magnon cross-pumping $W_{4B,h}(\phi_p) = \sum_{k\neq h} P_{hk} a_k^2$ against the $2A$ pumping sideband (`cross4b` prior-excitation downward shift + `slPhase1`) and powers the cursor-anchored Zoom/H/V/AUTO Arnold diagram controls. |
| `stepLlg` | [spin-waves/index.html](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html) | **Adaptive RK4 + Brown FDT Noise** | Integrates $(a_{R,m}, a_{I,m})$ via classical RK4 with nonlinear-rate-bounded adaptive substepping (`maxSubH ∈ [0.08, 0.40]`), Brown FDT thermal kicks, exact `SPEC_DT` sampling, and smooth `FAST→SS` pacing (`140×0.35 = 49 τ/frame`). |
| `seedSpins`, `exciteEigenmode`, `pluckAt` | [L1361–L1440](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1361-L1440) | **d'Aquino Eq. (S12) + Init** | Zero/deterministic/Gaussian thermal seeding and exact symplectic projection of spatial pointer plucks. |

---

## 5. Moving Time Plot, 1.5–2.3 GHz Spectrum, Tensor Inspector, Unclipped Arnold & Disk Rendering

| Function | Lines | Origin | Mathematical / Engineering Basis |
| :--- | :--- | :--- | :--- |
| `drawModeTimePlot` | [L1475–L1640](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1475-L1640) | **Moving Time-Window Telemetry** | Plots all nonzero modal amplitudes $|a_h(t)|$ over `[t - W, t]`, uses hysteresis + wall-`dt`-bounded upward/downward exponential Y-axis slewing (`plotYTarget`, `plotYScale`), and debounces `#modeLegend` chips (`300 ms`). |
| `updateSpectrum` | [L1650–L1811](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1650-L1811) | **High-Res 4096-Pt FFT + UI** | Zero-mean Hann-windowed 4096-pt radix-2 FFT (`3.42 MHz/bin` + quadratic peak interpolation), toggling between **Mode Band (`1.5–2.3 GHz`)** with all 29 eigenmode frequency markers/needles and **Wide Harmonics (`0–6.5 GHz`)**. |
| `updateEigenmodeReadout`, `buildNikMatrix`, `syncNikConfigUi`, `updateNikMatrix` | [L1813–L1925](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1813-L1925) | Telemetry & Tensor Inspector | Updates HUD badge and renders interactive $6\times 6$ matrix inspector (`N_hk`, `|P_hk|`, `K+_hk`, `S_det`). |
| `computeArnoldHMaxMT`, `getArnoldRange`, `drawArnoldTongues`, `drawHInputView`, `selectArnoldAt` | [L1976–L2293](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L1976-L2293) | **d'Aquino Unclipped 29-Mode Arnold Map** | Renders all 29 modes (`ALL 29` `3.38–4.18 GHz` or `BULK` `3.88–4.16 GHz`) with dynamic mT ceiling (fixing `1.40 mT` clipping), dashed static reference vs instantaneous `4A`-shifted/`4B`-reduced curves, and `h(t)` Tukey playhead. |
| `drawScene` | [L2295–L2380](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L2295-L2380) | **Exact Tensor Profile Rendering** | Renders strictly dynamic tensor perturbations $(\delta m_y, \delta m_z) = (m_y - m_{0y}, m_z - m_{0z})$ clipped to the circular YIG disk. |

---

## 6. 7-Stage Benchmark / Calibration / Verification Protocol, Slow Guide & DB Diagnostics

| Function | Lines | Origin | Mathematical / Engineering Basis |
| :--- | :--- | :--- | :--- |
| `relPhaseDeg`, `syncUi` | [L2438–L2490](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L2438-L2490) | State Sync | Wraps $\Delta\phi = \phi_2 - \phi_1$ cleanly into `[-180°, 180°]` (preserved across Tone 1 ↔ Tone 2 swap) and syncs UI + band toggles. |
| `getDbConfig`, `saveDbConfig`, `resetDbConfig`, `swLog`, `formatDbDiagText`, `buildEvaluationReport`, `listIdbSessions`, `deleteIdbSession`, `uploadSwLog`, `testDbConnection` | [L2626–L3010](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L2626-L3010) | **Evaluation Telemetry & DB Vault** | Configurable Supabase REST endpoint/key/JWT (`localStorage`), structured JSON reports (`spin-waves-evaluation-v1`) attaching `calibration`, `verification`, `benchmark`, and `lastDbAttempt`, explicit HTTP `401 / 42501` RLS diagnostics + copyable SQL fix, `GET + POST` connection test, and IndexedDB session vault management. |
| `applyPreset`, `swapTonesAction`, `fastForwardSteadyState`, `handleSaveLogDb` | [L3119–L3425](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L3119-L3425) | Presets & Protocol Actions | Sets presets (`nonrecip` calibrated for both `0 K` and `300 K` deterministic/thermal reliability), handles swaps, fast-forwarding, and DB/IDB saves. |
| `VERIFY_STAGES`, `runVerificationStage`, `runFullVerificationProtocol`, `runThroughputBenchmark`, `inspectLiveCalibration`, `probe4Avs4BCrossPump` | [L3495–L4115](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L3495-L4115) | **Interactive Verification & Calibration Engine** | 7-stage automated physics/calibration/throughput suite (`Cal 1`–`Bench 7`), live simulation-vs-d'Aquino inspector, and `4A` shift vs `4B` 4-magnon cross-pumping ablation probe inside `#verifyModal`. |
| `frame` | [L4281–L4355](file:///home/antigravity/kfwbnh5y25-source/spacetime-grid-/spin-waves/index.html#L4281-L4355) | Main Loop | Real-time slow-motion RK4 loop (`1× SLOW`), time-window trace sampling, 4096-pt FFT spectrum update, FPS telemetry, and 10 s IndexedDB autosave. |

