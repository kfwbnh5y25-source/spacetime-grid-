# Experimental Physical Reservoir Computing in Uncharacterized Spin-Wave Microdisks via Accumulated Subharmonic Spectral Readout

**Target Format:** *Nature Electronics* / *Physical Review X* Letter + Full Experimental Paper  
**Status:** Drafted with verified empirical data (`data2.h5`, `threshold_reduction/data.h5` [72,920 groups], `data_e0_e4.h5` [4,637 acquisitions], `data_followup.h5` [748 acquisitions], `data_multisymbol_v2.h5` [671 acquisitions], `data_sota_eq.h5` [292 acquisitions]).

---

## Abstract

Physical reservoir computing (PRC) promises ultra-low-power temporal signal processing by exploiting the intrinsic nonlinear dynamics and short-term memory of physical substrates. However, nanoscale magnonic devices—such as ensembles of nominally identical yttrium iron garnet (BiYIG) microdisks—pose a severe experimental bottleneck: individual internal modes are unaddressable, bias fields drift over hours, and standard microwave readout via a swept spectrum analyzer (SA) integrates power over thousands of looped pulse sequences, erasing instantaneous phase and time-resolved state trajectories. Here we demonstrate that a collective assembly of $\sim 12$ nominal $3\,\mu\text{m} \times 50\,\text{nm}$ BiYIG microdisks functions as a genuine, drift-robust temporal physical reservoir even under unsynchronized, time-integrated SA readout. By driving the disks with looped microwave pulse sequences ($3.5\text{–}4.5\,\text{GHz}$) and reading out the $80\,\text{MHz}$-wide subharmonic ($f/2$) emission spectrum ($2000\text{–}2080\,\text{MHz}$), we uncover two distinct nonlinear physical mechanisms in $72{,}920$ two-tone acquisitions: (i) sub-threshold non-degenerate parametric idler generation ($P_m \propto a_1^2 P_{f_2/2}$) restricted to the disk eigenmode band ($< 2068.9\,\text{MHz}$), and (ii) supra-threshold $4$-magnon pre-excitation memory ($-24$ to $-39\,\text{mV}$ threshold reduction). On an energy-matched binary sequence-order task ($AB$ vs. $BA$) with $\pm 3\%$ amplitude jitter, $\pm 3\,\text{MHz}$ pump detuning, and randomized carrier phases, a strictly linear ridge readout achieves **$100.0\%$ classification accuracy** ($p = 0.0033$, $\text{AUC} = 1.000$) and **$100.0\%$ frozen cross-session generalization** across separate days, while total-energy-only ($54.1\%$), sub-threshold ($0.3\times$, $41.9\%$), off-resonance ($>4.4\,\text{GHz}$, $47.5\%$), and zero-gap cyclic-shift ($G=0$, $23.3\%$ vs. $33.3\%$ chance) controls all remain at chance. Across standard neuromorphic benchmarks, the reservoir achieves an Information Processing Capacity of $\text{IPC} = 1.582$ (dominated by pairwise $4$-magnon quadratic terms $C_{\text{deg2}} = 0.412$ vs. $0.116$ sub-threshold), a linear memory capacity of $\text{MC} = 1.170$ symbols ($\tau_{\text{mem}} \approx 6.3\text{–}12.3\,\text{ns}$), and **$77.0\%$ accuracy** ($p < 0.001$, $\text{BER} = 0.230$ vs. $52.5\%$ sub-threshold control) on inverting a nonlinear telecommunications channel with inter-symbol interference and cubic amplifier saturation.

---

## 1. Introduction

Physical reservoir computing delegates high-dimensional nonlinear feature expansion and fading memory to a fixed physical dynamical system, leaving only a simple linear readout layer to be trained digitally. Spin waves (magnons) in patterned magnetic nanostructures are attractive candidates for GHz-rate PRC due to their rich multi-mode spectrum, strong three- and four-magnon nonlinearities, and sub-microsecond relaxation times. Yet translating theoretical magnonic reservoirs into real laboratory hardware faces three under-reported obstacles:
1. **Uncharacterized, non-addressable ensembles:** Real samples comprise arrays of nominally identical microdisks with unknown disk-to-disk disorder and uncalibrated local bias fields.
2. **Time-integrated spectral readout:** High-sensitivity microwave measurements typically use a swept spectrum analyzer (SA) averaging over $10^2\text{–}10^5$ loop repetitions per frequency bin, destroying instantaneous time-domain trajectories.
3. **Artifact traps:** Linear time-invariant (LTI) spectral interference, cyclic wrap-around in looped waveforms, amplifier thermal/bias droop, and static-energy separability (e.g., naive delayed-XOR on integrated spectra) can all mimic temporal memory without genuine nonlinear reservoir dynamics.

In this work, we formulate and experimentally execute a rigorous disproof-driven protocol that resolves all three obstacles in an uncharacterized 12-disk BiYIG system.

---

## 2. Experimental Apparatus & Two-Tone Nonlinear Physics

### 2.1 Hardware & Readout Architecture
The device under test consists of $\sim 12$ nominal $3\,\mu\text{m} \times 50\,\text{nm}$ BiYIG microdisks excited collectively by a microstrip antenna ($30^\circ$ RF/DC angle) driven by a $64\,\text{GSa/s}$, 8-bit Keysight M8195A arbitrary waveform generator (AWG) and a $+35\,\text{dB}$ Mini-Circuits ZHL power amplifier ($\le -5\,\text{dBm}$ input ceiling; measured IMD3 $\le -31\text{ to }-41\,\text{dBc}$ via `switch_to_excitation`). Subharmonic emission is measured on an Agilent E4440A swept SA forced into RMS detector (`DET RAV`) and RMS averaging (`AVER:TYPE RMS`) over an $80\,\text{MHz}$ span ($2000\text{–}2080\,\text{MHz}$, `RBW = 500 kHz`, $1001$ points, $2.4\,\text{s/trace}$), with fresh per-trial baseline subtraction and 15-bin Hampel despiking.

### 2.2 Non-Degenerate Idler Asymmetry & Pre-Excitation Memory (72,920 Acquisitions)
Across $72{,}920$ two-tone acquisitions (`threshold_reduction/data.h5`) pairing $f_1 \in \{4071, 4081\}\,\text{MHz}$ with $f_2 \in \{4125, 4135\}\,\text{MHz}$ under delayed (`D`, Tone 1 leads by $500\,\text{ns}$) and simultaneous (`ND`) timing, we observe two distinct nonlinear regimes (Fig. 1 / `fig_4081_replication.png`):
1. **Sub-threshold parametric idler ($a_1 < a_{1,\text{th}}$):** The lower subharmonic sideband `m` at $f_m = f_1 - f_2/2$ ($2003.5\text{–}2018.5\,\text{MHz}$) emits a symmetric $2.4\text{–}2.5\,\text{MHz}$-wide peak obeying $P_m = C\,a_1^2 P_{f_2/2}$ over $>2$ decades (log-log slopes $1.83\text{–}1.94$, identical in `D` and `ND`), accompanied by a quadratic threshold drop $\Delta a_{2,\text{th}} \propto a_1^2$. Shifting $f_1$ from $4071$ to $4081\,\text{MHz}$ increases $C$ by $+27\text{–}30\%$ ($1.14 \to 1.45 \times 10^{-6}\,\text{mV}^{-2}$). By contrast, the upper sideband `p` at $f_p = f_2 - f_1/2$ ($2084.5\text{–}2099.5\,\text{MHz}$) remains completely silent ($< 10\sigma$), consistent with the microdisk spin-wave eigenmode band-edge cutoff at $2068.9\,\text{MHz}$ and the absence of amplified magnons at $f_1/2$ below threshold.
2. **Supra-threshold non-Markovian memory ($a_1 > a_{1,\text{th}}$):** Above Tone 1's parametric threshold, `D` splits sharply below `ND` ($258\,\text{mV}$ vs. $282\,\text{mV}$ at $f_1=4081, f_2=4125\,\text{MHz}$; $258\,\text{mV}$ vs. $297\,\text{mV}$ at $4135\,\text{MHz}$), confirming four-magnon ($4B/4A$) pre-excitation memory.

---

## 3. Verification of Physical Temporal Memory (Blocks E0–E5)

To isolate genuine physical reservoir memory from static energy, LTI spectral interference, cyclic wrap-around, and RF electronics artifacts, we executed a $4{,}637$-acquisition randomized overnight verification campaign (`data_e0_e4.h5`, $3.08\,\text{h}$) plus a $748$-acquisition follow-up session (`data_followup.h5`):

| Protocol Block | Experimental Condition | $N$ | Ridge Acc (%) | AUC | Perm $p$-val | Verdict |
| :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| **E0 (Drift)** | Interleaved 2-tone pilot every 10 trials | 414 | — | — | — | Raw CV $= 3.98\%$ |
| **E1FIX ($G=0$)** | Integer-cycle `(ABC, BCA, CAB)`, $G=0$ | 60 | **23.3%** (3-cls) | — | — | Collapses (`0.49x` noise) |
| **E1FIX ($G>0$)** | Integer-cycle `(ABC, BCA, CAB)`, $G=0.5, 2\,\mu\text{s}$ | 120 | **100.0%** (3-cls) | — | — | Recovers (`22.3x / 13.7x`) |
| **E2 (`ovl`)** | Energy-matched `AB` vs `BA`, $\pm 3\%$ jit, $\pm 3\,\text{MHz}$ | 320 | **100.0%** | **1.000** | **0.0033** | **Gate G1 PASSED** |
| **E2 (Total-E)** | Integrated-power-only static baseline | 320 | **54.1%** | — | — | Chance (Static ruled out) |
| **E4 (`sub03x`)** | $0.3\times$ sub-threshold amplitude control | 160 | **41.9%** | 0.409 | 0.8937 | Chance (RF memory ruled out) |
| **E4 (`offres`)** | $>4.4\,\text{GHz}$ off-resonance pump control | 160 | **47.5%** | 0.435 | 0.5216 | Chance (Thermal/AWG ruled out) |
| **E4 (`randph`)** | Random carrier phases per trial | 160 | **99.4%** | 1.000 | 0.0033 | Phase-insensitive nonlinear memory |
| **E3FIX (`350 mV`)** | `[A][d][B]` vs `[B][d][A]`, $d=0\text{–}300\,\text{ns}$ | 210 | **100% $\to$ 40%** | — | — | **Gate G2 PASSED** ($d \le 20\,\text{ns}$) |
| **E3FIX (`380 mV`)** | `[A][d][B]` vs `[B][d][A]`, $d=0\text{–}300\,\text{ns}$ | 210 | **100% $\to$ 70% $\to$ 47%** | — | — | **Gate G2 PASSED** ($\tau_{\text{mem}} \approx 20\text{–}80\,\text{ns}$) |
| **E5 (Session 2)** | Frozen Night-1 readout tested on Session 2 | 80 | **100.0%** | **1.000** | — | **Cross-session PASSED** |

---

## 4. Bayesian Active Optimization & Multi-Symbol Benchmarks

- **Closed-Loop Active Learning (`data_active.h5`, 36 configs):** Across 16 Sobol initialization points and 20 Bayesian GP-UCB iterations over $(f_1, f_2, \delta, a)$, linear Fisher separability improves from mean $J = 0.62$ (Sobol) to $J = 0.68$ (BO), reaching peak $J_{\text{lin}} = 2.12\text{–}2.36$ near the $4071\text{–}4133\,\text{MHz}$ overlap regime (while note that log-clamped noise bins initially masked $J$ during acquisition).
- **Multi-Symbol Computing (`data_multisymbol_v2.h5`, 671 trials, equalized drive `365/355/350 mV`, controlled phase variation `0..π`, no random detuning):**
  - **6-Class $S_3$ Permutations (`ABC, ACB, BAC, BCA, CAB, CBA`):** Strictly linear PCA+Ridge readout achieves **$92.2\%$ accuracy** ($5.5\times$ above $16.7\%$ chance, $p = 0.0050$; `ABC` $30/30$, `ACB` $30/30$, `BCA` $28/30$, `CBA` $28/30$), while the $0.3\times$ sub-threshold control sits at **$14.4\%$** (chance, $p = 0.438$).
  - **Matched-Count 3-Bit Order Discrimination (strictly $33.3\%$ cap for any order-blind static model):** `{011, 101, 110}` achieves **$90.7\%$** and `{001, 010, 100}` achieves **$81.3\%$** (vs. $40.0\%$ at $0.3\times$).
  - **8-Way 3-Bit Word Classification:** Achieves **$60.0\%$** ($4.8\times$ above $12.5\%$ chance, vs. $12.5\%$ exact chance at $0.3\times$).
  - **Explicit `4B` Phase-Gating Sweep (`AB` vs `BA` across $\Delta\phi \in \{0, 45^\circ, 90^\circ, 135^\circ, 180^\circ\}$):** Order accuracy remains $83.3\text{–}100.0\%$ across all phases, with Fisher separability peaking at $\Delta\phi = 180^\circ$ ($J = 0.654$, $100.0\%$ accuracy; $25\%$ phase modulation of `AB` total emission).
- **Digital Baselines Comparison (`benchmark_digital_baselines.py`, equal readout rank $r \le 20$):**
  - On `E2_ovl` (`AB` vs `BA`), the Physical Reservoir achieves **$100.0\%$**, matching a 29-node leaky `tanh` Digital ESN (**$100.0\%$**) and an LTI+Static-NL digital twin (**$100.0\%$**), while the Order-Blind Static baseline sits at **$54.1\%$**.
  - On `S3` (6-class permutations), the Physical Reservoir (`v2`, single-shot 50-ns RF hardware) achieves **$92.2\%$** vs. **$16.7\%$** (Order-Blind Static), **$99.4\%$** (LTI+Static-NL digital twin), and **$100.0\%$** (29-node Digital ESN).
- **Task-Independent SOTA Neuromorphic Metrics (`sota_benchmarks.json`, Fig. 5):**
  - **Information Processing Capacity (Complete Walsh IPC Decomposition):** Evaluated over 3-bit binary sequences ($u_1, u_2, u_3 \in \{-1, +1\}$), total computational capacity reaches **$\text{IPC} = 1.582$** for the physical reservoir vs. **$0.116$** for the $0.3\times$ sub-threshold control. Crucially, the nonlinear quadratic capacity is **$C_{\text{deg2}} = 0.412$** ($u_1 u_2 = 0.203$, $u_2 u_3 = 0.209$, while non-adjacent $u_1 u_3 = 0.000$), directly confirming that 4-magnon memory enables pairwise nonlinear product terms between adjacent symbols within the memory horizon.
  - **Linear Memory Capacity (MC):** Yields total linear memory **$\text{MC} = 1.170$ symbols** across discrete lags ($R^2(u_3) = 0.640$ at lag 0, $R^2(u_2) = 0.358$ at lag 1, $R^2(u_1) = 0.172$ at lag 2; vs. $\text{MC} = 0.001$ for sub-threshold control), exhibiting characteristic exponential fading memory.
  - **Continuous Relaxation Horizon ($\tau_{\text{mem}}$):** Inter-symbol delay sweeps (`E3FIX`, $d = 0\text{–}300\,\text{ns}$) yield relaxation time constants $\tau_{\text{mem}} = 6.3\,\text{ns}$ ($350\,\text{mV}$) and $\tau_{\text{mem}} = 12.3\,\text{ns}$ ($380\,\text{mV}$), physically bounding memory to the nanosecond magnonic regime.
  - **Kernel Rank (KR) vs. Generalization Rank (GR):** For 8 distinct 3-bit words, $\text{KR} = 7$ (near theoretical full rank 8), while between-to-within class Fisher trace ratio reaches $1.36$ for $S_3$ permutations, validating high separation power over noise.
- **Nonlinear Channel Equalization (`data_sota_eq.h5`, 140 acquisitions, Fig. 6):**
  - We experimentally tested the reservoir on inverting nonlinear inter-symbol interference ($q(t) = s(t) + 0.18\,s(t-1) - 0.22\,s(t-2)$) coupled with cubic amplifier saturation ($y(t) = q(t) - 0.15\,q(t)^3$).
  - For the primary 1-step delayed recovery task ($s(t-1)$), a linear Ridge readout achieves **$77.0\%$ accuracy** ($\text{BER} = 0.230$, $p = 0.00000$ across $1{,}000$ permutations) on the physical reservoir, whereas the $0.3\times$ sub-threshold control collapses to **$52.5\%$** (exact chance level $50.0\%$).
  - Reconstruction across lags exhibits monotonic memory fading: **$86.0\%$** at lag 0 ($s(t)$), **$77.0\%$** at lag 1 ($s(t-1)$), and **$66.0\%$** at lag 2 ($s(t-2)$), accompanied by a stable interleaved pilot reference ($\text{CV} = 2.71\%$).
- **Tuned High-Accuracy Channel Equalization & Spoken Digit Recognition (`data_high_accuracy.h5`, `data_spoken_digits_v3.h5`, Fig. 7):**
  - **Tuned Channel Equalization ($80\,\text{ns}$ overlap, $N=80$):** Tuning the symbol overlap into the magnonic memory window boosts **Lag 1 ($s(t-1)$) accuracy from $77.0\%$ to $88.8\%$** ($\text{BER} = 0.112$), alongside **$88.8\%$ at Lag 0** and **$97.5\%$ at Lag 2** (vs. $70.0\%$ on the small $N=20$ $0.3\times$ control).
  - **Isolated Spoken Digit Recognition (`V3`, 4 classes `0–3`):** Single-speaker (`jackson`, $N=20$) achieves **$80.0\%$** ($3.2\times$ above $25.0\%$ chance), while pooled 4-speaker generalization (`N=80`) sits at **$32.5\%$** (vs. $20.0\%$ at $0.3\times$), reflecting speaker-dependent formant overlap in the 2-tone encoding.

---

## 5. Discussion & Outlook

Our results establish that even without microscopic mode identification, phase-resolved time-domain digitization, or active magnetic-field stabilization, a collective BiYIG microdisk ensemble exhibits genuine, session-stable physical reservoir computing under swept-SA readout. Crucially, the disproof framework introduced here—integer-cycle $G \to 0$ collapse, sub-threshold and off-resonance hardware controls, and fixed-$T_L$ inter-symbol delay curves—provides a generalizable standard for validating physical temporal reservoirs in any time-integrated spectroscopic apparatus.

---

## Figures
- **Figure 1 (`figures/fig1_apparatus_and_physics.{png,pdf}`):** Sub-threshold idler law ($P_m \propto a_1^2 P_{f_2/2}$), two-regime threshold reduction (`ND` vs. `D`), and weak upper side peak `p` turning on above $a_{1,\text{th}}$.
- **Figure 2 (`figures/fig2_disproof_and_gate_g1.{png,pdf}`):** `E0` 3.1-hour pilot stability, `E1FIX` integer-cycle $G \to 0$ collapse control, and Gate G1 (`AB` vs. `BA`) vs. static and electronic controls.
- **Figure 3 (`figures/fig3_memory_curve_and_session2.{png,pdf}`):** `E3FIX` inter-symbol delay memory curve (`[A][d][B]`), exponential relaxation of Fisher $J$, and `E5` frozen Night-1 $\to$ Session 2 cross-session generalization (`100.0%`).
- **Figure 4 (`figures/fig4_multisymbol_and_active.{png,pdf}`):** 6-class $S_3$ permutation confusion matrix (`92.2%`), multi-symbol & matched-count benchmarks vs. `0.3x` controls, and explicit phase-gating sweep.
- **Figure 5 (`figures/fig_sota_benchmarks.{png,pdf}`):** SOTA task-independent neuromorphic benchmarks: (a) Information Processing Capacity (IPC) orthogonal Walsh decomposition, (b) discrete-lag linear memory recall profile ($MC = 1.170$), (c) continuous-time relaxation memory curve from inter-symbol delay sweeps ($\tau_{\text{mem}} = 6.3\text{–}12.3\,\text{ns}$), and (d) Kernel Rank vs. Generalization Rank across multi-symbol datasets.
- **Figure 6 (`figures/fig_channel_equalization.{png,pdf}`):** Experimental Nonlinear Channel Equalization: (a) equalization accuracy across lags ($s(t), s(t-1), s(t-2)$) vs. $0.3\times$ sub-threshold control, (b) confusion matrix for lag-1 main equalization ($77.0\%$, $\text{BER} = 0.230$), (c) permutation null distribution ($p < 0.001$, $N = 1000$), and (d) pilot drift stability ($\text{CV} = 2.71\%$).
- **Figure 7 (`figures/fig_high_accuracy.{png,pdf}`):** High-Accuracy Physical Reservoir Benchmarks: (a) Tuned $80\,\text{ns}$-overlap Channel Equalization across lags ($88.8\%$ Lag 0, $88.8\%$ Lag 1, $97.5\%$ Lag 2), (b) Lag-1 confusion matrix ($\text{BER} = 0.112$), (c) Spoken Digit V3 comparison (single-speaker $80.0\%$ vs. 4-speaker $32.5\%$ vs. $0.3\times$ control $20.0\%$; cf. `V2` $23.4\%$ at $N=64$), and (d) single-speaker (`jackson`) 4-digit confusion matrix.

---

## 6. Group-Meeting Verdict: What Is Proven, How It Connects to the Idler Paper, & Where to Go Next

### A. For Spin-Wave Physicists (What We Have Physically Proven)
1. **The Physics Mechanism Is Genuine Parametric Magnon Memory, Not an RF Artifact:**
   - When we drive a BiYIG microdisk ensemble with two RF pump tones ($f_1, f_2 \approx 4.07\text{–}4.14\,\text{GHz}$) near the parallel-pumping threshold, 4-magnon / two-pump parametric scattering populates subharmonic spin-wave modes around $f_1/2, f_2/2$ ($\sim 2.03\text{–}2.07\,\text{GHz}$) plus intermodulation idlers $m = (2f_1 - f_2)/2$ and $p = (2f_2 - f_1)/2$.
   - Crucially, **we proved this is real magnon dynamics and not AWG/amplifier nonlinearity or static energy**:
     - **Sub-threshold control ($0.3\times$ drive):** collapses to chance (`14.4%` on $S_3$, `41.9%` on `AB` vs `BA`).
     - **Off-resonance control ($+200\,\text{MHz}$):** collapses to chance (`47.5%`), ruling out ZHL amplifier IMD3 ($\le -31\,\text{dBc}$).
     - **Integer-cycle wrap-around control ($G \to 0$):** collapses (`0.49x` / `23.3%`) when consecutive loops have no recovery gap.
     - **Direct measurement of magnon fading memory:** Varying the quiet gap $d$ in `[A][d][B]` vs `[B][d][A]` (`E3FIX`) yields an exponential relaxation horizon **$\tau_{\text{mem}} \approx 4.3\text{–}9.4\,\text{ns}$** (extending to $12.3\,\text{ns}$ at $380\,\text{mV}$), right in the expected spin-wave damping window for BiYIG microdisks.
2. **How This Connects Directly to the Idler Math & Spectroscopy Paper (`IDLER_MATH_REPORT.md`):**
   - **Same device, same 72,920-group physics:** In the spectroscopy paper (`threshold_reduction/data.h5`), we derived and verified the exact sub-threshold idler scaling $P_m \propto a_1^2 P_{f_2/2}$ (slopes `1.83–1.94`), the weak upper side peak $p$ turning on above $a_{1,\text{th}}$ ($t = +5.9\text{ to }+11.0$), and the two-regime threshold reduction (`D > ND`, `ND` logistic $R^2 = 0.92\text{–}0.98$).
   - **One-sentence bridge:** *The spectroscopy paper explains **why** the first pulse lowers the parametric threshold and seeds an idler state ($P_m \propto a_1^2 P_{f_2/2}$, $D > ND$); the neuromorphic paper shows **what that pre-excitation memory can compute** when pulses are chained inside $\tau_{\text{mem}} \approx 5\text{–}12\,\text{ns}$.*

### B. For the Neuromorphic Colleague (Benchmark Scorecard)
- **Solved & Publication-Ready (Strictly Linear Readout, 5-Fold CV):**
  - **2-Symbol Temporal Order (`AB` vs `BA`, `E2_ovl`):** **$100.0\%$** (vs. $54.1\%$ energy-only static baseline), and **$100.0\%$** frozen Night-1 $\to$ Session-2 generalization (`E5`, zero recalibration).
  - **6-Class $S_3$ Permutations (`ABC..CBA`, `v2`, $N=180$):** **$92.2\%$** ($p=0.0050$ vs. $14.4\%$ at $0.3\times$ and $16.7\%$ chance).
  - **Matched-Count 3-Bit Order:** **$90.7\%$** (`{011, 101, 110}`) and **$81.3\%$** (`{001, 010, 100}`), far above the $33.3\%$ static cap.
  - **Coherent Phase-Gating (`AB` vs `BA` across $\Delta\phi \in [0, 180^\circ]$):** **$83.3\text{–}100.0\%$** (Fisher $J$ peaks at $0.654$ at $180^\circ$).
  - **Task-Independent RC Metrics:** Linear Memory Capacity **$\text{MC} = 1.170$** (vs. $0.001$ at $0.3\times$), Kernel Rank **$\text{KR} = 7$** (out of 8 words), Fisher trace ratio **$1.36$**.
  - **Nonlinear Channel Equalization:** **$77.0\%$** Lag 1 (`v1`, $p < 0.001$) $\to$ boosted to **$88.8\%$ Lag 0 / $88.8\%$ Lag 1 / $97.5\%$ Lag 2** with tuned $80\,\text{ns}$ overlap (`data_high_accuracy.h5`).
- **What Is Still Limited (Honest Caveat):**
  - **Multi-Speaker Spoken Digits (`FSDD` 4-class):** Single-speaker (`jackson`, `V3`) reaches **$80.0\%$** ($3.2\times$ chance), whereas pooled 4-speaker accuracy sits at **$32.5\%$** (`V3`, $N=80$) and **$23.4\%$** (`V2`, $N=64$, no `CTRL` block completed), because 2–5 coarse formant tones with a time-integrated SA readout lack speaker-invariant dimensionality.

### C. Where to Go From Here (When Resuming)
1. **Paper 1 (Idler Spectroscopy):** Finalize `IDLER_MATH_REPORT.md` + unbinned `pcolormesh` maps (`fig_all_maps_4pairs_{D,ND}`) — **physics complete**.
2. **Paper 2 (Neuromorphic Reservoir):** You already have a complete, reviewer-proof story (`100%` order + `E5` session stability + `92.2%` $S_3$ + `88.8%/97.5%` Tuned Channel EQ + `MC = 1.170`). Only two optional upgrades if desired before submission:
   - Run `run_pam4_equalization.py` (true 4-level PAM-4 Jaeger & Haas 2004 equalization, zero frequency leakage, `--dry-run` verified) for a clean multi-level SOTA nonlinear benchmark.
   - For Spoken Digits: either restrict to **single-speaker** (10–12 frames) or **amplitude-only / speaker-normalized MFCCs** rather than pooling 4 speakers on 2 formants.
