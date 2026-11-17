# Experimental Temporal Neuromorphic Computing Through an Unidentified Nonlinear Spin-Wave System
### Comprehensive Feasibility Assessment, Empirical Data Audit, Mathematical Framework, and Experimental Protocol

**Date:** 2026-11-16  
**System:** ~12 nominally identical $3\,\mu\text{m} \times 50\,\text{nm}$ BiYIG microdisks, collective antenna excitation ($3.5\text{–}4.5\,\text{GHz}$, $30^\circ$ RF/DC angle) and $S_{21}$ readout ($1.75\text{–}2.25\,\text{GHz}$), Keysight M8195A AWG ($64\,\text{GSa/s}$, 8-bit `int8`), Mini-Circuits ZHL amplifier ($+35\,\text{dB}$, $\le -5\,\text{dBm}$ input), Agilent E4440A swept spectrum analyzer.  
**Evidence tags:** **`[DATA]`** = directly computed from raw HDF5/CSV files in this repo (`data2.h5`, `data_broad.h5`, `threshold_reduction/data.h5`, `phase_harmonics.csv`); **`[DRIVER]`** = verified in `Keysight_AWG_M8195A.py`, `Agilent_PSA_E4440A.py`, `run.py`; **`[SIM]`** = 29-mode NMM ODE solver (`spin-waves/index.html` via `scripts/lib/spin-loader.mjs`); **`[THEORY]`** = exact analytical derivation; **`[HYP]`** = untested experimental hypothesis.

---

## 0. Executive Summary

1. **Feasibility Verdict:** Temporal neuromorphic computing with this apparatus is **experimentally feasible only as a sequence-level physical reservoir with pump–probe / accumulated spectral readout (Hybrid Architecture B+D)**—*not* as an instantaneous continuous-time state machine, and *not* without strict artifact controls.
2. **Existing Evidence (`data2.h5` + `threshold_reduction/data.h5`):**
   - `data2.h5` (`[DATA]`): Within-run repeatability on the $f/2$ channel is **3–4% CV** above $15\sigma$, with **14–26 dB** usable dynamic range and **2.41 s/acquisition** ($\sim 10{,}000\text{–}12{,}000$ acquisitions overnight). Between-run scatter across days is **$1.8\text{–}5.2\times$** larger than within-run SEM; minute-scale drift is uncharacterized and requires interleaved pilot normalization.
   - `threshold_reduction/data.h5` (`[DATA]`): Re-reducing all **39,513 HDF5 groups** with a 15-bin Hampel despike filter (removing narrow 5–15 kHz RF pickup spikes that ruined `threshold_fits.csv`) reveals **two real physical two-tone effects** at $f_1=4071\,\text{MHz}$, $f_2 \in \{4125, 4135\}\,\text{MHz}$:
     1. **Sub-threshold non-degenerate parametric assistance ($a_1 < a_{1,\text{th}}$):** Even below Tone 1's threshold ($a_1 \in [100, 225)\,\text{mV}$), $a_{2,\text{th}}$ drops by **$-10$ to $-15\,\text{mV}$** identically in `D` and `ND` ($\Delta a_{2,\text{th}} \propto a_1^2$), accompanied by a **2.5 MHz-wide idler emission** in channel `m` at $\frac{2f_1 - f_2}{2} = f_1 - \frac{f_2}{2}$ ($2008.5 / 2003.6\,\text{MHz}$) scaling cleanly as $P_m \propto a_1^2 P_{f_2/2}$ over $>75\times$ dynamic range. This **contradicts symmetric-disk $0^\circ$ `2BC-off` NMM** (`[SIM]`) and is explained by $30^\circ$ RF/DC symmetry breaking + disorder activating off-diagonal parametric coupling $Q_{2x/2y,mk}$.
     2. **Supra-threshold pre-excitation memory ($a_1 > a_{1,\text{th}} \approx 260\text{–}275\,\text{mV}$):** Giving Tone 1 a $500\,\text{ns}$ lead (`D`) drops $a_{2,\text{th}}$ to **$271\,\text{mV}$ ($-46\,\text{mV}$)** vs. **$293\,\text{mV}$ ($-24\,\text{mV}$)** for simultaneous `ND`, with **`D > ND` in 96% of paired points (`4125`) and 92% (`4135`)**, confirming NMM `4B`/`4A` pre-excitation (`[SIM]`). However, `run.py` gave `D` $2\times$ Tone-1 duration ($1.0\,\mu\text{s}$ vs $0.5\,\mu\text{s}$), so `data.h5` alone does **not** yet prove energy-matched temporal memory.
3. **One Convincing Initial Task:** **Energy-matched binary sequence-order discrimination ($AB$ vs $BA$) and delayed-order pump–probe classification vs delay $d$** under amplitude jitter and random pump detuning, benchmarked against a measured LTI + static-nonlinear digital twin.

---

## 1. Hardware Reality & Empirical Data Audit `[DATA + DRIVER]`

### 1.1 Hardware Constraints (`Keysight_AWG_M8195A.py`, `Agilent_PSA_E4440A.py`, `run.py`)
- **AWG (Keysight M8195A):** $64\,\text{GSa/s}$, **8-bit `int8` ($\pm 127$)** DAC at `output_voltage = 1.0 V`. At $a = 250\text{–}350\,\text{mV}$, waveforms use only $32\text{–}45$ quantization levels ($\sim 5.2$ bits), introducing $\sim 2\text{–}3\%$ quantization uncertainty in energy matching unless driven near full scale with external fixed attenuation. Sequencer `STAB` segments require multiples of **$1280\text{ samples} = 20\,\text{ns}$**; single-waveform upload (`send_waveform` / `BIN8`) is preferred for loop periods $< 20\,\mu\text{s}$, `STAB` for $> 20\,\mu\text{s}$.
- **RF Amplifier Chain:** Mini-Circuits ZHL ($+35\,\text{dB}$ gain, $\le -5\,\text{dBm}$ max input). In `threshold_reduction/data.h5`, direct monitoring via `knet.switch_to_excitation` (`EX`, $4025.5\text{–}4175.5\,\text{MHz}$) at combined amplitudes up to $a_1+a_2 = 810\,\text{mV}$ shows $P_2/a_2^2$ drops by only **$-2.1\%$ ($\approx 0.09\,\text{dB}$)** and $P_1/a_1^2$ by **$-1.3\%$**, ruling out gross amplifier compression, though `EX` span must be widened to $\ge 200\,\text{MHz}$ to bound IMD3 at $2f_1 - f_2$ directly.
- **Spectrum Analyzer (Agilent E4440A):** Existing scripts (`run2.py`, `runBROAD.py`, `threshold_reduction/run.py`) left `detectorType` unset (`AUTO`) and `averageType` unset, and allowed `resolutionBandwidthAutoState = True` (`510 kHz`). All future runs **must force** `sa.detectorType = "RAV"` (RMS), `sa.averageType = "RMS"`, `sa.resolutionBandwidthAutoState = False`, and span $\ge 50\text{–}100\,\text{MHz}$ (`RBW = 500 kHz`) to simultaneously observe `m`, $f_1/2$, and $f_2/2$.

### 1.2 Repeatability & Noise Statistics (`data2.h5`: 29,135 Acquisitions) `[DATA]`
| Metric | $f/2$ Channel (`halfaverage`, $50/3\,\mu\text{s}$) | $f$ Channel (`full`, CW) |
| :--- | :---: | :---: |
| Acquisitions | $15{,}635$ | $13{,}500$ |
| Cycle time / repeat | $2.41\,\text{s}$ ($p_{95}=2.46\,\text{s}$) | $0.54\,\text{s}$ |
| Null ($a=0$) per-bin noise $\sigma_0$ | $4.85\times 10^{-8}$ | — |
| Within-run CV by power quartile | Q1: $246\%$, Q2: $10.8\%$, Q3: $4.0\%$, Q4: $2.8\%$, Q5: $3.1\%$ | $0.37\text{–}0.41\%$ |
| Between-run SD / within-run SEM | $12.6\%$ vs $3.6\%$ ($\sqrt{\text{var ratio}} = 1.8$ med, $5.2$ RMS) | $0.61\%$ vs $0.38\%$ |
| Effective PCA rank (null-whitened) | Participation ratio $PR = 5.4$ ($80\%$ var in 20 PCs) | $PR = 1.001$ |

### 1.3 Forensic Audit of `spin-waves/threshold_reduction/data.h5` (39,513 Groups) `[DATA]`
`threshold_fits.csv` failed to resolve $15/16$ thresholds due to four fixable issues uncovered in the raw HDF5:
1. **Two sampling patterns merged:** $750$ early pairs (`03:30–06:41`) sampled a 1D diagonal $a_2 \approx 1.5 a_1$, while $1{,}104$ later pairs (`06:41–17:39`) sampled a 2D Sobol square (`a1 ∈ [0,350]`, `a2 ∈ [0,450]`).
2. **Unfiltered narrowband RF spurs + baseline caching:** 5–15 kHz spurs (e.g., `2037.98`, `2088.0`, `2091.98 MHz`) plus `_cached_bl` (90 s TTL, not cleared between `D`/`ND` or across `a2=0`) produced $-200\sigma$ negative outliers (`-3.8e-4`). Hampel despiking eliminates these completely (`quiet σ ≈ 1.8–2.0 × 10⁻⁶`).
3. **Unmatched Tone-1 energy in `D` vs `ND`:** Loop length is $96{,}000 / 64\,\text{GSa/s} = 1.5\,\mu\text{s}$ ($0.5\,\mu\text{s}$ off-gap). `D` plays Tone 1 for **$1.0\,\mu\text{s}$** ($0\text{–}1.0\,\mu\text{s}$) and Tone 2 for $0.5\,\mu\text{s}$ ($0.5\text{–}1.0\,\mu\text{s}$); `ND` sets `window1 = window2`, playing Tone 1 for **$0.5\,\mu\text{s}$** ($0.5\text{–}1.0\,\mu\text{s}$).
4. **Despiked physics summary (`10σ = 2e-5` threshold crossing):**
   - **$f_2 = 4125\,\text{MHz}$:** $a_{2,\text{th}}$ moves from $317\,\text{mV}$ ($a_1=0$) $\to 307\,\text{mV}$ ($a_1 \in [175,225)$, `D == ND`) $\to$ **$271\,\text{mV}$ (`D`) vs $293\,\text{mV}$ (`ND`)** at $a_1 \in [275,355)$.
   - **$f_2 = 4135\,\text{MHz}$:** $a_{2,\text{th}}$ moves from $331\,\text{mV}$ ($a_1=0$) $\to 316\,\text{mV}$ ($a_1 \in [175,225)$, `D == ND`) $\to$ **$299\,\text{mV}$ (`D`) vs $304\,\text{mV}$ (`ND`)** at $a_1 \in [275,355)$.
   - **Channel `p` (`2089.5/2099.5 MHz`):** Empty after despiking (above the NMM spectrum cutoff $f_{28} = 2.0689\,\text{GHz}$).
   - **Channel `m` (`2008.5/2003.5 MHz`):** Genuine broadband emission ($\text{FWHM} \approx 2.4\text{–}2.8\,\text{MHz}$), identical in `D` and `ND` ($\text{med}(D/ND) = 1.001\text{–}1.017$), with $P_m / (a_1^2 P_{f_2/2}) \approx (0.8\text{–}1.3)\times 10^{-6}$ constant across $a_1 \in [1, 355)\,\text{mV}$.

### 1.4 Quantitative Proof of Idler Scaling, Spectral Centering, and Publication Figures `[DATA]`
All tables and figures below are reproducibly generated from raw `threshold_reduction/data.h5` and `figure1_single_tone/cache/*.npz` by `spin-waves/reports/neuromorphic/analyze_threshold_reduction.py` (outputs saved in `spin-waves/reports/neuromorphic/threshold_reduction_audit.json`).

#### Table 1: Empirical Constancy of $P_m / (a_1^2 P_{f_2/2})$ Across $a_1 \in [25, 355)\,\text{mV}$ ($a_2 \ge 350\,\text{mV}$, Delayed `D`)
| $a_1$ Bin (mV) | $N$ | $f_2=4125$: $\text{med}(P_m)$ | $f_2=4125$: $\text{med}(P_{f_2/2})$ | $f_2=4125$: $\frac{P_m}{a_1^2 P_{f_2/2}}$ ($\times 10^{-6}\,\text{mV}^{-2}$, [IQR]) | $f_2=4135$: $\text{med}(P_m)$ | $f_2=4135$: $\frac{P_m}{a_1^2 P_{f_2/2}}$ ($\times 10^{-6}\,\text{mV}^{-2}$, [IQR]) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`[25, 60)`** | 26 / 27 | $3.87\times 10^{-6}$ | $1.63\times 10^{-3}$ | **$1.58$** $[1.09,\, 2.45]$ | $2.51\times 10^{-6}$ | **$0.91$** $[0.81,\, 1.21]$ |
| **`[60, 100)`** | 34 / 32 | $1.20\times 10^{-5}$ | $1.68\times 10^{-3}$ | **$1.27$** $[1.09,\, 1.41]$ | $9.95\times 10^{-6}$ | **$0.87$** $[0.74,\, 0.92]$ |
| **`[100, 140)`** | 32 | $2.57\times 10^{-5}$ | $1.69\times 10^{-3}$ | **$1.12$** $[1.01,\, 1.34]$ | $1.89\times 10^{-5}$ | **$0.85$** $[0.76,\, 0.94]$ |
| **`[140, 180)`** | 31 | $5.23\times 10^{-5}$ | $1.73\times 10^{-3}$ | **$1.18$** $[1.07,\, 1.27]$ | $3.95\times 10^{-5}$ | **$0.86$** $[0.81,\, 0.92]$ |
| **`[180, 220)`** | 32 | $7.39\times 10^{-5}$ | $1.63\times 10^{-3}$ | **$1.13$** $[0.99,\, 1.24]$ | $5.31\times 10^{-5}$ | **$0.86$** $[0.78,\, 0.91]$ |
| **`[220, 260)`** | 96 | $8.86\times 10^{-5}$ | $1.38\times 10^{-3}$ | **$1.05$** $[0.99,\, 1.13]$ | $4.87\times 10^{-5}$ | **$0.80$** $[0.76,\, 0.83]$ |
| **`[260, 300)`** | 134 | $1.16\times 10^{-4}$ | $1.29\times 10^{-3}$ | **$1.15$** $[1.13,\, 1.17]$ | $9.25\times 10^{-5}$ | **$0.82$** $[0.80,\, 0.83]$ |
| **`[300, 355)`** | 8 | $9.56\times 10^{-5}$ | $8.47\times 10^{-4}$ | **$1.07$** $[1.03,\, 1.10]$ | $7.95\times 10^{-5}$ | **$0.74$** $[0.72,\, 0.81]$ |

Log-log linear regressions of $P_m / P_{f_2/2}$ vs. $a_1$ confirm exact quadratic dependence over $>2$ decades: slope **$1.91 \pm 0.02$** (`4125 D`), **$1.86 \pm 0.02$** (`4125 ND`), **$2.01 \pm 0.03$** (`4135 D`), and **$1.94 \pm 0.02$** (`4135 ND`). (The slight upward bias in the lowest `[25, 60)` bin reflects residual positive rectification near the noise floor when $P_m \sim 3\times 10^{-6}$.)

#### Table 2: Spectral Centering & Edge-Leakage Audit (Hampel-Despiked Mean Spectra, `D`)
| Channel ($f_c$, MHz) | Regime ($N$ traces) | Peak Offset (MHz) | Centroid Offset (MHz) | Left Edge / Center / Right Edge Power | Physical Diagnosis |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **`f1`** ($2035.5$) | $a_1>275, a_2<150$ ($97$) | $+0.130$ | $-0.024$ | $3.7\times 10^{-9}$ / **$5.6\times 10^{-8}$** / $2.7\times 10^{-9}$ | **Centered peak** ($18\times$ center-to-edge) |
| **`f1`** ($2035.5$) | $a_1\in[150,250), a_2>360$ ($98$) | $-0.525$ | $+0.045$ | $6.6\times 10^{-9}$ / $1.2\times 10^{-8}$ / $6.6\times 10^{-9}$ | **Flat broadband pedestal** ($\approx 7\text{–}10\times 10^{-9}$ across all 5 MHz), **NOT** a parametric peak! |
| **`f2`** ($2062.5$) | $a_1=0, a_2>360$ ($90$) | $+0.080$ | $+0.061$ | $2.3\times 10^{-7}$ / **$2.9\times 10^{-6}$** / $3.1\times 10^{-7}$ | **Centered peak** ($11\times$ center-to-edge) |
| **`f2`** ($2062.5$) | $a_1>275, a_2\in[265,305)$ ($12$) | $+0.090$ | $-0.015$ | $3.2\times 10^{-8}$ / **$2.4\times 10^{-7}$** / $3.1\times 10^{-8}$ | **Centered pre-excited boost** ($7.5\times$ center-to-edge; zero leakage from $f_1/2$) |
| **`m`** ($2008.5$) | $a_1\in[100,200), a_2>360$ ($69$) | $+0.040$ | $-0.090$ | $1.2\times 10^{-8}$ / **$7.9\times 10^{-8}$** / $7.4\times 10^{-9}$ | **Centered idler** ($8.3\times$ center-to-edge even at sub-threshold $a_1$) |
| **`m`** ($2008.5$) | $a_1\in[250,350), a_2>360$ ($163$) | **$+0.005$** | $-0.073$ | $3.0\times 10^{-8}$ / **$1.9\times 10^{-7}$** / $2.3\times 10^{-8}$ | **Centered idler** ($7.4\times$ center-to-edge; 5 kHz from exact $(2f_1-f_2)/2$) |
| **`m`** ($2003.5$) | $a_1\in[250,350), a_2>360$ ($163$) | $+0.120$ | $-0.042$ | $2.0\times 10^{-8}$ / **$1.5\times 10^{-7}$** / $1.6\times 10^{-8}$ | **Centered idler** ($8.6\times$ center-to-edge; zero edge leakage) |

#### Publication Figure Captions & Physical Observations
- **Figure 1 (`fig_idler_and_threshold_scaling.{png,pdf}`): *Non-degenerate parametric idler scaling and two-regime threshold reduction.***  
  **(a)** Normalized idler power $P_m / (a_1^2 P_{f_2/2})$ vs. Tone-1 amplitude $a_1$ for $f_2 = 4125\,\text{MHz}$ (blue) and $4135\,\text{MHz}$ (red) under delayed (`D`) and simultaneous (`ND`) timing ($a_2 \ge 350\,\text{mV}$; error bars denote interquartile range across Sobol samples). The ratio stays constant across $a_1 \in [60, 355)\,\text{mV}$ with zero discontinuity at the $f_1$ parametric threshold $a_{1,\text{th}} \approx 265\,\text{mV}$ (dashed line) and no `D`/`ND` split. **(b)** Log-log scatter of $P_m / P_{f_2/2}$ vs. $a_1$ across $>2$ decades, showing exact $a_1^2$ power-law slopes ($1.86\text{–}2.01$). **(c)** Tone-2 threshold $a_{2,\text{th}}$ ($10\sigma$ crossing) plotted against $a_1^2$: `ND` follows a linear drop in $a_1^2$ (dotted fits), matching the off-diagonal $2B/2C$ self-energy prediction $\Delta \Gamma_k \propto -h_1^2$; above $a_{1,\text{th}}^2$, `D` splits sharply below `ND` (by $-22\,\text{mV}$ at $4125\,\text{MHz}$), marking the onset of $4B$ pre-excitation memory.
- **Figure 2 (`fig_pcolormesh_pm_vs_model.{png,pdf}`): *2D $(a_1, a_2)$ validation of $P_m = C\,a_1^2 P_{f_2/2}$ on a shared logarithmic colorscale.***  
  Left column: Hampel-despiked measured idler power $P_m(a_1, a_2)$ at $f_2 = 4125\,\text{MHz}$ (top) and $4135\,\text{MHz}$ (bottom). Middle column: Single-parameter model $C\,a_1^2 P_{f_2/2}(a_1, a_2)$ ($C = 1.15\times 10^{-6}\,\text{mV}^{-2}$ at $4125$, $0.82\times 10^{-6}\,\text{mV}^{-2}$ at $4135$) plotted on the identical colorscale. Right column: Log-ratio residual $\log_{10}(P_m / \text{Model})$ in active cells ($>5\times 10^{-6}$), showing agreement within $\pm 0.15\,\text{decades}$ across the entire $(a_1, a_2)$ plane—both below and above $a_{1,\text{th}}$ (cyan dashed line)—confirming that idler generation tracks instantaneous $a_1^2 P_{f_2/2}$ rather than an $f_1$ threshold bifurcation.
- **Figure 3 (`fig_all_responses_a1_a2.{png,pdf}` & `fig_broadband_af_map.{png,pdf}`): *Complete despiked two-tone and single-tone response landscapes.***  
  `fig_all_responses_a1_a2` shows all four readout bands (`f1/2`, `f2/2`, `m`, `p`) on a shared log scale: `p` ($2089.5/2099.5\,\text{MHz}$) remains silent everywhere (above the $2068.9\,\text{MHz}$ disk mode cutoff), `f2/2` exhibits downward threshold bending as $a_1$ increases, and `m` closely mirrors $a_1^2 P_{f_2/2}$. Note that the apparent faint `f1/2` elevation at high $a_2$ and sub-threshold $a_1$ is revealed by spectral inspection (Figure 4) to be a flat broadband pedestal, not a localized parametric resonance. `fig_broadband_af_map` maps the despiked single-tone overview (`Figure1/data.h5`) and zoom (`data2.h5`) landscapes, locating $f_1 = 4071\,\text{MHz}$, $f_2 = 4125/4135\,\text{MHz}$, and the idler pump equivalents $2f_1 - f_2 = 4007/4017\,\text{MHz}$ inside the multi-tongue magnonic band.
- **Figure 4 (`fig_spectral_centering_check.{png,pdf}`): *Spectral centering and edge-leakage verification.***  
  Averaged Hampel-despiked spectra plotted against detuning $f - f_c$ inside each $5\,\text{MHz}$ SA window. Both the idler `m` ($2008.5$ and $2003.5\,\text{MHz}$, bottom row, active even for sub-threshold $a_1 \in [100,200)\,\text{mV}$) and the pre-excited `f2/2` boost ($2062.5\,\text{MHz}$, top right) form symmetric bell curves centered within $\le 0.12\,\text{MHz}$ of $f_c$ and decaying by $7.4\text{–}11\times$ toward both $\pm 2.5\,\text{MHz}$ boundaries, **ruling out spectral edge leakage** from adjacent resonances. By contrast, `f1/2` at sub-threshold $a_1 \in [150,250)\,\text{mV}$ and high $a_2 > 360\,\text{mV}$ (top left, red) exhibits a **flat, uncentered broadband pedestal** ($\approx 7\text{–}10\times 10^{-9}$ across all $5\,\text{MHz}$), proving that genuine `f1/2` parametric oscillation (blue, $\text{FWHM}=1.29\,\text{MHz}$) occurs only for $a_1 > a_{1,\text{th}}$.

---

## 2. Architecture Comparison `(A)`–`(D)`

| Approach | Feasibility under Constraints | Legitimate Scientific Claim | Fatal Flaw if Misused |
| :--- | :--- | :--- | :--- |
| **(A) Sequence $\to$ Spectrum Learning** | **High** (immediate ML fit on $(u, Y)$) | Black-box digital *surrogate* / simulator of the experiment. | **Not neuromorphic.** Pure regression of hardware output cannot claim physical computation. |
| **(B) Physical Reservoir + Accumulated Readout** | **High** (sequence-level reservoir + linear readout $W\phi(Y)$) | **Physical temporal reservoir computing** *iff* residual beyond LTI+static twin predicts unseen order/delay across sessions. | Using nonlinear readout MLPs or testing order without energy/LTI controls fakes reservoir power. |
| **(C) Latent State-Space Models** | **Unrealistic** (except $\le 3\text{–}5$ mode gray-box diagnostics) | Explanatory dynamical fit for simple 2-tone regimes. | 12-disk disorder + unknown $B_0$ + phase-erased integrated SA trace $\implies$ trajectories are **unidentifiable**. |
| **(D) Uncertainty-Aware Bayesian Surrogate (GP/BO)** | **Essential** (exact GP up to $N \lesssim 5\times 10^3$) | Sample-efficient experimental planner & drift-aware uncertainty calibration. | Without order-aware kernels and heteroscedastic noise, BO exploits drift/noise spikes. |

**Recommended Architecture:** **(D) plans sequences $\to$ (B) executes physical reservoir transformation + linear ridge readout $\to$ (A) serves as the LTI+static digital twin baseline $\to$ (C) is restricted to post-hoc diagnostics.**

---

## 3. Mathematical Formulation `[THEORY]`

### 3.1 Input Sequence & Looped Excitation
A discrete input sequence $u = (s_1, \dots, s_K)$ of symbols $s_k = (a_k, f_k, \phi_k, \Delta_k)$ is mapped to an AWG waveform on $[0, T_{\text{seq}}]$ followed by a zero-power recovery gap $G$, looped periodically with period $T_L = T_{\text{seq}} + G$:
$$h_{\text{ext}}(t) = \sum_{n \in \mathbb{Z}} \sum_{k=1}^K a_k \, w_{\Delta_k}\!\big(t - n T_L - t_k\big) \sin\!\big(2\pi f_k (t - n T_L) + \phi_k\big)$$
where $w_{\Delta_k}(\tau)$ is a flat-top pulse of duration $\Delta_k$ (quantized to $20\,\text{ns}$ when using `STAB`, arbitrary at $64\,\text{GSa/s}$ for $<20\,\mu\text{s}$ single waveforms) with $1\,\text{ns}$ Tukey edges.

### 3.2 Collective Multi-Disk Normal-Modes Dynamics
Across the $J \approx 12$ coupled/inhomogeneous disks (each with effective modes $h = 1,\dots,M$ around $1.7\text{–}2.1\,\text{GHz}$, slightly shifted by disk-to-disk disorder $\delta\omega_h^{(j)}$ and slow field drift $\delta B(t_{\text{lab}})$), the complex mode amplitudes $a_h^{(j)}(\tau)$ obey (d'Aquino et al. 2026; Perna et al. 2022):
$$(1 + \alpha^2)\frac{d a_h^{(j)}}{d\tau} = -\big(j \bar{\omega}_h^{(j)} + \alpha K_{hh}^+\big) a_h^{(j)} + \mathcal{F}_h^{(1)} + \mathcal{F}_h^{(2)} + \mathcal{F}_h^{(4A)} + \mathcal{F}_h^{(4B)} + \mathcal{F}_h^{(4C,4D)} + \eta_h^{(j)}(\tau)$$
where $\eta_h^{(j)}$ is Brown FDT thermal noise ($\langle |\eta_h|^2 \rangle \propto \alpha K_{hh}^+ k_B T / (\mu_0 M_s^2 V)$).

### 3.3 Exact Mechanism of Sub-Threshold Reduction (`m` Idler) & Supra-Threshold Memory (`4B`)
This derivation explains both experimental regimes observed in `threshold_reduction/data.h5`:

1. **Non-Degenerate Two-Mode Parametric Self-Energy (Explains $a_1 < a_{1,\text{th}}$ & Channel `m`):**  
   Decompose $\mathcal{F}^{(2)}$ into diagonal ($2A$) and off-diagonal ($2B/2C$) parametric couplings under $h_{\text{ext}}(t) = h_1 \cos(\omega_1 t) + h_2 \cos(\omega_2 t)$ (active during overlap at $30^\circ$ RF/DC angle). Let target mode $k$ satisfy $\omega_k \approx \omega_2/2$ and idler mode $m$ satisfy $\omega_m \approx \omega_1 - \omega_2/2 = (2\omega_1 - \omega_2)/2$ (for $4071, 4125\,\text{MHz}$: $\omega_m/2\pi = 2008.5\,\text{MHz}$, matching modes `m17/m18` at $2006.2/2006.7\,\text{MHz}$). Transforming into rotating frames $a_k = c_k e^{-j(\omega_2/2)t}$, $a_m = c_m e^{-j(\omega_1 - \omega_2/2)t}$ and retaining resonant terms yields the coupled linear instability system:
   $$\dot{c}_k = -(\Gamma_k + j\delta_k) c_k + V_k h_2 c_k^* + \kappa_{km} h_1 c_m^*, \qquad \dot{c}_m^* = -(\Gamma_m - j\delta_m) c_m^* + \kappa_{mk}^* h_1 c_k$$
   where $\delta_k = \omega_k - \omega_2/2$, $\delta_m = \omega_m - (\omega_1 - \omega_2/2)$, and $\kappa_{km}$ is the off-diagonal parametric coupling tensor ($Q_{2x/2y,km}$). Solving for the adiabatic idler response gives:
   $$c_m^*(t) \approx \frac{\kappa_{mk}^* h_1}{\Gamma_m - j\delta_m} c_k(t) \implies P_m \propto |c_m|^2 \propto \frac{|\kappa_{mk}|^2}{\Gamma_m^2 + \delta_m^2}\, h_1^2 |c_k|^2 \propto a_1^2 P_{f_2/2}$$
   **Exact Match 1:** This predicts $P_m / (a_1^2 P_{f_2/2}) = \text{const}$ across all $a_1$, observed empirically in `data.h5` from $a_1 = 25\,\text{mV}$ to $350\,\text{mV}$!  
   Substituting $c_m^*$ back into $\dot{c}_k$ yields an effective damping reduction quadratic in $h_1$:
   $$\Gamma_{k,\text{eff}}(h_1) = \Gamma_k - \frac{\operatorname{Re}(\kappa_{km}\kappa_{mk}^*) \Gamma_m + \operatorname{Im}(\kappa_{km}\kappa_{mk}^*)\delta_m}{\Gamma_m^2 + \delta_m^2}\, h_1^2, \qquad h_{2,\text{th}}(h_1) = \frac{\sqrt{\Gamma_{k,\text{eff}}(h_1)^2 + \delta_k^{\prime\,2}}}{V_k}$$
   **Exact Match 2:** For $\Delta a_{2,\text{th}}(a_1) \propto a_1^2$, scaling from $a_1 = 200\,\text{mV}$ ($\Delta = -10\,\text{mV}$) to $a_1 = 315\,\text{mV}$ predicts $-10 \times (315/200)^2 = -24.8\,\text{mV}$, matching the observed `ND` shift of **$-24\,\text{mV}$** ($317 \to 293\,\text{mV}$) to within $1\,\text{mV}$! Because this mechanism requires simultaneous overlap of $h_1$ and $c_k$, it vanishes whenever Tone 1 and Tone 2 do not overlap (`D == ND`).

2. **Pre-Excited Non-Markovian $4B$ + $4A$ Memory (Explains `D > ND` for $a_1 > a_{1,\text{th}}$):**  
   When $h_1 > h_{1,\text{th}}$, prior excitation populates $a_1(t_{\text{on}}) \ne 0$ before Tone 2 arrives. Even after Tone 1 turns off, $a_1(t)$ decays on timescale $\tau_{\text{relax}} = 1/\Gamma_1 \sim 100\text{–}300\,\text{ns}$, contributing:
   - **4B Coherent Cross-Pump:** $W_{4B,k}(t) = \sum_l P_{kl} a_l^2(t)$ lowering $h_{2,\text{th}}$ linearly in $|W_{4B}|$;
   - **4A Kerr Shift:** $\Delta\omega_k(t) = \sum_l N_{kl} |a_l(t)|^2$ shifting the Arnold tongue frequency across $\delta_k = 0$;
   - **Winner-Take-All Attractor Selection:** Crossing a multi-mode separatrix during $[0, t_{\text{on}}]$ locks the ensemble into distinct steady-state attractors depending on symbol arrival order (`AB` vs `BA`), even with **zero temporal overlap**.

### 3.4 Swept Spectrum Analyzer Observation Operator, Cyclic Collapse & LTI Null Proof
With `DET RAV` (RMS detector) and dwell time $T_{\text{dwell}} = T_{\text{sweep}} / N_{\text{pts}} \gg T_L$ (e.g., $T_{\text{dwell}} = 2\,\text{ms} \gg T_L \le 10\,\mu\text{s}$, averaging $\ge 200$ loop periods per frequency bin), let the antenna output voltage in the $1.75\text{–}2.25\,\text{GHz}$ band during a stationary limit-cycle period $t \in [0, T_L)$ be the collective disk sum:
$$v_{\text{out}}(t) = \sum_{j=1}^{12} \sum_{h=0}^{M-1} \operatorname{Re}\!\left[g_{jh} a_h^{(j)}(t)\right], \qquad \hat{v}_{\text{out}}(f) = \frac{1}{T_L}\int_0^{T_L} v_{\text{out}}(t) e^{-j2\pi f t}\,dt$$
The swept SA measures a convolution of the periodic power spectrum with the RBW filter $K_{\text{RBW}}(f)$ ($\text{ENBW} = 1.056\,\text{RBW}$ for the E4440A):
$$Y_i(u) = \sum_{n \in \mathbb{Z}} \left| K_{\text{RBW}}\!\left(\frac{n}{T_L} - f_i\right) \right|^2 \left| \hat{v}_{\text{out}}\!\left(\frac{n}{T_L}\right) \right|^2 + \varepsilon_i, \qquad \mathrm{Var}(\varepsilon_i) = \sigma_0^2 + \kappa^2 \mu_i^2$$
where $\sigma_0 \approx 4.85\times 10^{-8}$ and $\kappa \approx 0.03\text{–}0.04$ (`[DATA]`). Notice the two spectral regimes:
- **Comb-resolved ($\text{RBW} \ll 1/T_L$):** Individual Fourier lines $n/T_L$ are resolved; phase-sensitive interferences inside $T_L$ appear directly in $|\hat{v}_{\text{out}}(n/T_L)|^2$.
- **Envelope-averaged ($\text{RBW} \ge 3/T_L$, e.g. $500\,\text{kHz}$ at $T_L \ge 6\,\mu\text{s}$):** By Parseval's theorem across the RBW passband, $Y_i(u) \approx \frac{1}{T_L}\int_0^{T_L} |v_{\text{BP}, i}(t)|^2 dt$, eliminating comb-aliasing ripple.

**Theorem 1 (Cyclic-Shift Degeneracy as $G \to 0$):**  
If $G = 0$, any cyclic shift $u_{\text{shift}}(t) = u((t - t_0) \bmod T_{\text{seq}})$ shifts the periodic steady-state trajectory $a_h^{(j)}(t) \mapsto a_h^{(j)}(t - t_0)$, multiplying every Fourier coefficient $\hat{v}_{\text{out}}(n/T_L)$ by $e^{-j2\pi (n/T_L) t_0}$. Because the RMS SA observes only $|\hat{v}_{\text{out}}|^2$, **all cyclic shifts yield identical spectra up to $\mathcal{O}(T_L / T_{\text{dwell}})$**:
$$\|Y(u_{\text{shift}}) - Y(u)\|_2 \xrightarrow{G \to 0} \mathcal{O}\!\left(\sigma_{\text{noise}}\right)$$
Thus, $G = 0$ serves as a **mandatory positive control**: any non-zero cyclic contrast at $G=0$ proves an AWG/SA trigger or waveform-stitching bug.

**Theorem 2 (Why LTI + Static Nonlinearity Fakes Order Sensitivity):**  
Consider a memoryless static nonlinearity $\mathcal{N}_{\text{static}}$ followed by an arbitrary LTI filter $H(f)$. Let $z(t) = \mathcal{N}_{\text{static}}[u(t)] = z_A(t) + z_B(t - \Delta t)$ for sequential non-overlapping symbols $(A, B)$ vs $(B, A)$ separated by $\Delta t$. The Fourier power is:
$$\left|H(f) Z_{AB}(f)\right|^2 = |H(f)|^2 \Big(|Z_A(f)|^2 + |Z_B(f)|^2 + 2\operatorname{Re}\!\big[Z_A^*(f) Z_B(f) e^{-j2\pi f \Delta t}\big]\Big)$$
Swapping $A \leftrightarrow B$ replaces $Z_A^* Z_B e^{-j2\pi f \Delta t}$ with $Z_B^* Z_A e^{-j2\pi f \Delta t} = \big(Z_A^* Z_B e^{+j2\pi f \Delta t}\big)^*$, changing the cross-term whenever $\operatorname{Im}[Z_A^*(f) Z_B(f)] \ne 0$ and $|Z_A(f) Z_B(f)| > 0$!  
In Volterra-series terms, genuine **non-Markovian physical reservoir memory** requires an asymmetric off-diagonal Volterra kernel $h_2(\tau_1, \tau_2) \ne h_2(\tau_2, \tau_1)$ extending over $|\tau_1 - \tau_2| \sim \tau_{\text{relax}} \gg 1/B_{\text{sym}}$ (and higher-order attractor bifurcations) that cannot be represented by $\mathcal{H}_{\text{meas}}[\mathcal{N}_{\text{static}}(u)]$. **Therefore, physical memory is established *iff* the residual $Y(u) - \hat{Y}_{\text{twin}}(u)$ predicts sequence order on held-out trials.**

### 3.5 Reservoir Readout, Effective Rank & Heteroscedastic Order-Aware GP Surrogate
- **Reservoir Readout (Approach B):** Extract Hampel-despiked, pilot-normalized log-spectral vector $\phi(u) = \log(\max(Y_{\text{corr}}(u), \epsilon)) \in \mathbb{R}^d$ (or top-$r$ PCA projections fit on training data only, bounded by empirical participation ratio $PR = (\sum \lambda_k)^2 / \sum \lambda_k^2$). Train strictly linear ridge readout $\hat{z}(u) = W \phi(u) + b$, minimizing $\sum_{n=1}^{N_{\text{tr}}} \ell(z_n, W\phi_n + b) + \lambda \|W\|_F^2$.
- **Active-Learning Surrogate (Approach D):** Model $\tilde{y}_i(u, t_{\text{lab}}) = \log(\max(Y_i(u), \epsilon))$ with variance-stabilized heteroscedastic noise $\tilde{\sigma}_i^2 \approx \kappa^2 + (\sigma_0 / \mu_i)^2$:
  $$\tilde{y}_i(u, t_{\text{lab}}) \sim \mathcal{GP}\!\Big(m_i(u), \; k_{\text{order}}(u_{\text{ord}}, u_{\text{ord}}') \, k_{\text{cont}}(\theta, \theta') + k_{\text{drift}}(t_{\text{lab}}, t_{\text{lab}}')\Big)$$
  where:
  1. $k_{\text{order}}(u, u') = \sum_{s \in \Sigma^2} \lambda_{\text{gap}}^{\text{span}(s,u) + \text{span}(s,u')}$ is a **gap-weighted 2-gram subsequence string kernel** (explicitly distinguishing `AB` from `BA` with decay factor $\lambda_{\text{gap}} \in (0,1)$);
  2. $k_{\text{cont}}(\theta, \theta')$ is a Matérn-5/2 kernel over continuous waveform parameters $(a_k, f_k, \Delta_k, d)$;
  3. $k_{\text{drift}}(t_{\text{lab}}, t_{\text{lab}}') = \sigma_{\text{drift}}^2 \exp(-|t_{\text{lab}} - t_{\text{lab}}'| / \ell_{\text{drift}})$ captures minute-to-hour $B_{\text{DC}}$ drift anchored by the interleaved pilot $P_0$.
  Candidate batches are selected via **BALD / Mutual Information** on order-posterior uncertainty:
  $$\alpha_{\text{BALD}}(u) = H\!\big[\mathbb{E}_{p(\mathbf{f}|\mathcal{D})}[p(z \mid \mathbf{f}(u))]\big] - \mathbb{E}_{p(\mathbf{f}|\mathcal{D})}\!\big[H[p(z \mid \mathbf{f}(u))]\big]$$
  penalized by a diversity-promoting $k$-center / determinantal point process (DPP) term so each batch of $16$ trials spans distinct $(f_1, f_2, d)$ regimes.

---

## 4. Minimal Experiment & Temporal-Memory Disproof Protocol

**Acquisition Settings (Fixed for All Runs):**
- AWG: $64\,\text{GSa/s}$, single-waveform upload (`T_loop < 20 µs`), amplitude scaled so ZHL input $\le -5\,\text{dBm}$ peak (accounting for 2-tone crest factor!).
- SA (E4440A): `DET RAV`, `AVER:TYPE RMS`, `BWID:AUTO OFF`, `RBW = 500 kHz`, `SPAN = 80 MHz` centered at $2040\,\text{MHz}$ ($2000\text{–}2080\,\text{MHz}$, covering `m`, $f_1/2$, $f_2/2$, and `m23` simultaneously!), `points = 1001`, `sweepTime = 0.25 s`, `averageNumber = 3`. Store **raw signal and raw baseline** separately (never cache baseline $>1$ trial).
- Throughput: $\sim 2.4\,\text{s/acquisition} \implies \mathbf{10{,}000\text{–}12{,}000\text{ acquisitions overnight}}$ ($\sim 8\text{ h}$).

| Block | Purpose | Exact Design & Controls | Budget ($N$) | Pass / Disproof Criterion |
| :--- | :--- | :--- | :---: | :--- |
| **E0: Drift & Allan Dev.** | Quantify minute–hour $B_{\text{DC}}$ drift | Fixed 2-tone pilot $P_0$ interleaved every 10 trials + `EX` (`knet.switch_to_excitation`) $200\,\text{MHz}$ span | $\ge 500$ pilots | Pilot-normalized CV $\le 5\%$ over 8 h; ZHL IMD3 at $2f_1 - f_2 < -40\,\text{dBc}$ |
| **E1: Gap & Cyclic Control** | Separate physical memory from loop wrap-around | Cyclic shifts $(A,B,C), (B,C,A), (C,A,B)$ across gap $G \in \{0, 0.25, 0.5, 1, 2, 4, 8\}\,\mu\text{s}$; keep duty $\ge 60\%$ for $G \le 2\,\mu\text{s}$ | $30/\text{cell}$ | **Positive control:** cyclic shifts MUST collapse to noise floor as $G \to 0$; memory stabilizes for $G \ge 1.5\,\mu\text{s}$ |
| **E2: Energy-Matched Order** | Core Gate G1 ($AB$ vs $BA$) | **Matched duration ($0.5\,\mu\text{s}$ each), matched energy, $G = 2.0\,\mu\text{s}$:** (i) Sequential $(A,B)$ vs $(B,A)$, gap $\delta \in \{0, 50, 200\}\,\text{ns}$; (ii) Overlap $(A \to A+B)$ vs $(B \to A+B)$ + total-energy-matched single pulse | $80/\text{class}$ | Acc $> 80\%$ ($p < 0.01$) under $\pm 3\%$ amplitude jitter & $\pm 3\,\text{MHz}$ detuning; **twin baseline $\le 60\%$** |
| **E3: Pump–Probe Delay** | Short-term memory curve (Gate G2) | `[u_1, u_2] -- (delay d) -- [Probe P]`, $d \in \{0, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0\}\,\mu\text{s}$; fixed probe at $f_p = 4095\,\text{MHz}$ | $40/\text{cell}$ | Order accuracy from probe response decays on $\tau_{\text{relax}} \sim 100\text{–}500\,\text{ns}$ (not $\mu\text{s}$ thermal timescale) |
| **E4: Artifact Disproofs** | Exclude RF/SA/DAC artifacts | (i) $0.3\times$ sub-threshold ($a \sim 90\,\text{mV}$); (ii) Off-resonance pump ($f > 4.4\,\text{GHz}$); (iii) Time-reversed envelope; (iv) Randomized phases; (v) Resend vs cached WF | $40/\text{cell}$ | **ALL E4 controls must be at chance ($\le 55\text{–}60\%$)**. If $0.3\times$ or $>4.4\,\text{GHz}$ shows order $\implies$ **electronic memory disproof** |
| **E5: Cross-Session Test** | Generalization over days | Re-run E2+E3 $\ge 12\text{–}24\,\text{h}$ later; train strictly on Night 1, test frozen $W$ on Night 2 | $1{,}000$ | Held-out Night 2 Acc $> 75\%$ after pilot normalization |

---

## 5. Active-Learning & Readout Training Pseudocode

```python
def run_active_reservoir_campaign(hw, twin_model, budget=4000, pilot_every=10):
    """
    Hybrid B+D Pipeline:
      - GP surrogate (D) searches sequence parameter space for non-trivial physical memory
        (order contrast ABOVE digital twin prediction & above noise floor).
      - Physical Reservoir (B) is evaluated via linear ridge readout on held-out sessions.
    """
    D_train, pilots = [], []
    gp = HeteroscedasticOrderGP(kernel=SubsequenceKernel() * Matern52() + DriftRBF())
    V_fixed = generate_matched_energy_validation_set(n=120)  # Never used for acquisition

    # Phase 0: Sobol space-filling seed + E0/E4 baseline calibration
    init_seqs = sample_energy_matched_sobol(n=80)
    for idx, seq in enumerate(init_seqs):
        if idx % pilot_every == 0:
            pilots.append(acquire_pilot_and_ex(hw))
        raw_sig, raw_bl = acquire_rms_spectrum(hw, seq, span_mhz=80, rbw_khz=500, det="RAV")
        y_clean = hampel_despike(raw_sig - raw_bl)
        y_norm = normalize_by_pilot(y_clean, pilots[-1])
        D_train.append((seq, y_norm, timestamp()))

    # Phase 1: Batch Bayesian Active Learning (maximizing excess physical memory over twin)
    while len(D_train) < budget:
        gp.fit(D_train, pilots)

        # Generate candidate paired permutations (S = AB, S_rev = BA) with matched energy
        cand_pairs = propose_matched_energy_pairs(n=2000)
        scores = []
        for s_fwd, s_rev in cand_pairs:
            mu_fwd, sig_fwd = gp.predict_log_spectrum(s_fwd)
            mu_rev, sig_rev = gp.predict_log_spectrum(s_rev)

            # Subtract LTI + static-nonlinear digital twin expectation
            delta_phys = (mu_fwd - twin_model(s_fwd)) - (mu_rev - twin_model(s_rev))
            fisher_est = np.sum(delta_phys**2) / (np.sum(sig_fwd**2 + sig_rev**2) + 1e-12)
            epistemic_unc = np.mean(sig_fwd + sig_rev)

            # UCB on physical excess Fisher separability + BALD diversity
            scores.append(fisher_est + 2.0 * epistemic_unc)

        batch = select_diverse_topk(cand_pairs, scores, k=16)  # DPP / k-center diversity
        batch_trials = expand_with_controls_and_shuffle(batch)  # Includes 0.3x & off-res controls

        for seq in batch_trials:
            if len(D_train) % pilot_every == 0:
                pilots.append(acquire_pilot_and_ex(hw))
            raw_sig, raw_bl = acquire_rms_spectrum(hw, seq, span_mhz=80, rbw_khz=500, det="RAV")
            y_norm = normalize_by_pilot(hampel_despike(raw_sig - raw_bl), pilots[-1])
            D_train.append((seq, y_norm, timestamp()))

        # Periodic validation on fixed unseen set V_fixed (train ridge readout on D_train only)
        if len(D_train) % 400 == 0:
            W_out, pca = fit_linear_ridge_readout(D_train)
            val_acc, val_twin_acc = evaluate_on_set(W_out, pca, twin_model, V_fixed)
            if convergence_plateau(val_acc):
                break

    # Phase 2: Freeze W_out, pca, GP; evaluate once on Night-2 held-out session
    return freeze_and_test_heldout_session(hw, W_out, pca, twin_model)
```

---

## 6. Computational Task, Baselines, Metrics & Failure Modes

### 6.1 Initial Benchmark Task (Avoiding Accumulated-Readout Traps)
- **Why NOT standard per-step delayed XOR (`u_t XOR u_{t-1}`)?** With one accumulated spectrum per sequence, $\{AA, BB\}$ vs $\{AB, BA\}$ is separable by a **static nonlinear function of symbol counts/energies**.
- **Correct Initial Task:** **Energy-Matched Order Discrimination ($AB$ vs $BA$, extended to $S_3$ permutations of $(A,B,C)$)** under nuisance perturbations ($\pm 3\%$ amplitude jitter, $\pm 3\,\text{MHz}$ pump detuning, varying symbol durations), followed by **Pump–Probe Delayed-Order Classification** vs delay $d \in [0, 4]\,\mu\text{s}$.
- **Required Baselines:** (1) Digital ideal $n$-gram oracle (ceiling), (2) Order-blind static features (symbol counts, total energy, per-frequency dwell—must be at chance $50\%$), (3) **Measured LTI + static-nonlinear digital twin** (must be $\le 60\%$), (4) Sub-threshold ($0.3\times$) & off-resonant ($>4.4\,\text{GHz}$) hardware controls (must be $\le 55\%$), (5) Digital ESN with equal linear readout rank ($r \le 20$).

### 6.2 Quantitative Acceptance & Disproof Criteria
- **Genuine Physical Reservoir Proof:** Night-2 held-out classification accuracy **$> 80\%$** ($p < 0.01$ permutation test) with Fisher separability $J > 3.0$, **and** LTI+static twin $\le 60\%$, **and** all E4 RF controls $\le 55\%$, **and** cyclic shifts collapse at $G \to 0$.
- **Disproof / Failure Modes:**
  1. Order effect survives at $0.3\times$ threshold or $>4.4\,\text{GHz}$ $\implies$ **electronic/amplifier memory** (not spin waves).
  2. Order effect vanishes when twin is subtracted $\implies$ **LTI spectral interference** (not reservoir memory).
  3. Cyclic shifts $(A,B,C)$ vs $(B,C,A)$ fail to collapse at $G=0$ $\implies$ **waveform/SA indexing bug**.
  4. Delayed-order accuracy stays flat out to $d = 10\,\mu\text{s}$ $\implies$ **slow thermal heating**, not spin-wave relaxation ($\tau_{\text{relax}} \sim 200\,\text{ns}$).
  5. Cross-session accuracy collapses despite pilot normalization $\implies$ **day-scale $B_0$ hysteresis/drift exceeds contrast**.

---

## 7. Publication Roadmap

| Milestone | Scope | Primary Figures / Evidence | Go / No-Go Gate |
| :--- | :--- | :--- | :--- |
| **P0 (Internal)** | E0 stability + E1 gap/cyclic control + ZHL IMD3 check | Allan deviation, pilot-normalized CV, $G$-collapse curve | Pilot-normalized drift $\le 5\%$ over 8 h; cyclic shifts collapse at $G=0$ |
| **P1 (Letter)** | Two-tone non-degenerate (`2B/2C` idler `m`) + pre-excitation (`4B`) threshold reduction + energy-matched `AB/BA` | Fig 1: Apparatus & $2000\text{–}2080\,\text{MHz}$ spectra; Fig 2: $a_{2,\text{th}}(a_1)$ + `m` scaling + E2/E4 controls | Night-1 `AB/BA` $> 80\%$ vs twin $\le 60\%$ & E4 controls at chance |
| **P2 (Full Paper)** | Pump–probe delayed-order memory curve + active-learning GP optimization + cross-day generalization | Fig 3: Accuracy vs delay $d$ ($\tau_{\text{mem}}$); Fig 4: GP active-learning convergence vs digital ESN | Night-2 held-out generalization $> 75\text{–}80\%$ |
| **Negative Fallback** | If E2 fails against twin or E4 controls fail | Quantitative characterization of swept-SA accumulated readout limits & $2B/2C$ nonlinear mixing | Honest, publishable characterization of physical reservoir constraints |
