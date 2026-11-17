# Pulse-length $\times$ amplitude map (dataset `3/3.h5`) and $a_1\times f_1$ two-tone phase sweep (`delay_f1_JULY_varyf1a1.h5`)

All numbers computed from the raw h5 files; scripts in `reports/mistral/tmp/an_*.py`; maps as `.npy` in the same folder.

---

## 1. Dataset A — `phase gating/3/3.h5`: single tone $f=4072\,$MHz, pulse length $l\times$ amplitude $a$

### 1.1 True waveform geometry (from `RUN.py`, the acquisition script of the $l\ge 20$ keys)

The `script_content` attribute stored in the file belongs to a **later, different session** (a single-pulse short-probe run, `avg=1`, single tukey pulse centered in an 800-sample record, $l$ in ns units $0$–$19.875$). It is **not** the script that produced the main grid. The main 4698-key grid comes from `phase gating/3/RUN.py`:

- $l_{\rm steps}=l$ AWG samples at 64 GSa/s; $l = 20,30,\dots,590\ \Rightarrow\ \tau_p = l/64$ ns $= 0.3125\ldots 9.2188$ ns, step $0.15625$ ns.
- One period $= 2\,l_{\rm steps}$ samples: **pulse** ($l_{\rm steps}$, tukey with taper $\alpha=2/l_{\rm steps}$ $\Rightarrow$ 1-sample tapers, effectively rectangular) followed by a **gap of exactly $l_{\rm steps}$** $\Rightarrow$ duty cycle 50 %, gap $=$ pulse length.
- `wf = np.concatenate((period,)*10)`: 10 identical pulse-gap periods (burst $=20l/64$ ns $=0.3125\ldots 184$ ns). The AWG loops this waveform continuously during the 0.5 s SA sweep, i.e. the mode sees a **continuous 50 % duty square gate** with period $2\tau_p$; the "10 pulses" set the burst length and the comb line spacing $1/T_{\rm burst}$, but the train repeats back-to-back.
- SA: center $2036$ MHz $=f/2$, span 75 MHz, 1001 pts (0.075 MHz/bin), RBW 300 kHz, avg 10, baseline static mask no-op (span 1998.5–2073.5 lies below the 2.10–2.17 GHz mask). Baseline re-measured every $\ge 90$ s with the AWG cleared (output off), so subtraction is clean.
- **90 stray keys** with $l\in[0,19]$ (short-probe session, $l$ in rounded ns) excluded.

Per-key extraction: `medfilt(spec,15)`; line peak = argmax of the filtered trace over the full span; adaptive window = contiguous filtered-trace excursion above $0.25\times$(peak $-$ median); integrated net power $P = \sum_{\rm window}(\text{raw spec}) - \langle B\rangle\,N_{\rm win}$; FWHM at half-height; broadband background $B$ = median of filtered trace (reported separately); noise $\sigma$ = MAD of raw trace; integrated SNR $=P/(\sigma\sqrt{N_{\rm win}})$.

### 1.2 Threshold $a_{\rm th}(l)$

Detection criterion: net window power $>5\,\sigma\sqrt{N_{\rm win}}$ (5$\sigma$ integrated) with persistence ($\ge 3$ consecutive amplitudes, removes all spurious low-$a$ picks). Result (mV; $*$ = not detected even at $a=1000$):

| $\tau_p$ (ns) | $a_{\rm th}$ | $\tau_p$ | $a_{\rm th}$ | $\tau_p$ | $a_{\rm th}$ | $\tau_p$ | $a_{\rm th}$ |
|---|---|---|---|---|---|---|---|
| 0.3125 | 590 | 2.500 | 420 | 4.688 | 360 | 6.875 | 340 |
| 0.4688 | 590 | 2.656 | 420 | 4.844 | 370 | 7.031 | 330 |
| 0.6250 | 580 | 2.812 | 410 | 5.000 | 350 | 7.188 | 330 |
| 0.7812 | 560 | 2.969 | 390 | 5.156 | 360 | 7.344 | 340 |
| 0.9375 | 530 | 3.125 | 400 | 5.312 | 340 | 7.500 | 330 |
| 1.0938 | 520 | 3.281 | 390 | 5.469 | 340 | 7.656 | 330 |
| 1.2500 | 480 | 3.438 | 390 | 5.625 | 350 | 7.812 | 330 |
| 1.4062 | 480 | 3.594 | 390 | 5.781 | 350 | 7.969 | 330 |
| 1.5625 | 460 | 3.750 | 390 | 5.938 | 350 | 8.125 | 320 |
| 1.7188 | 460 | 3.906 | 370 | 6.094 | 350 | 8.281 | 320 |
| 1.8750 | 440 | 4.062 | 370 | 6.250 | 350 | 8.438 | 340 |
| 2.0312 | 440 | 4.219 | 360 | 6.406 | 350 | 8.594 | 330 |
| 2.1875 | 420 | 4.375 | 370 | 6.562 | 340 | 8.750 | 330 |
| 2.3438 | 420 | 4.531 | 360 | 6.719 | 330 | 8.906 | 330 |
|  |  |  |  |  |  | 9.062 | 330 |
|  |  |  |  |  |  | 9.219 | 320 |

- Global saturated threshold $a_\infty = 320$ mV (reached for $\tau_p\gtrsim 8$ ns; $340\pm10$ mV level already at $\tau_p\approx 5$–$7$ ns).
- Single-parameter build-up model $a_{\rm th}(l)=a_\infty/\sqrt{1-e^{-l/\tau_b}}$: **$\tau_b = 2.15\pm0.08$ ns** (rms 28 mV over 57 points, $l\ge0.47$ ns). The model reproduces the whole curve including the fast saturation below $\tau_p\sim4$ ns.

### 1.3 Ring-up fits $P(l) = c_0 + P_\infty(1-e^{-l/\tau_r})$

Only amplitudes $a=420$–$520$ mV are trusted (clean sustained line detection and monotone rising $P(l)$):

| $a$ (mV) | $\tau_r$ full range [ns] | $\tau_r$ from onset [ns] | $\chi^2_{\rm red}$ | SNR at 9.22 ns |
|---|---|---|---|---|
| 420 | $14.4\pm1.1$ | 2.67 | 38 | 63 |
| 430 | $9.3\pm0.4$ | 2.07 | 55 | 68 |
| 440 | $5.9\pm0.2$ | 1.65 | 61 | 67 |
| 450 | $5.4\pm0.1$ | 1.84 | 59 | 92 |
| 460 | $4.8\pm0.1$ | 2.38 | 76 | 97 |
| 470 | $4.1\pm0.1$ | 2.11 | 59 | 114 |
| 480 | $2.70\pm0.04$ | 1.19 | 59 | 106 |
| 490 | $2.14\pm0.03$ | 0.78 | 60 | 104 |
| 500 | $2.15\pm0.03$ | 0.65 | 93 | 99 |
| 510 | $1.83\pm0.03$ | 0.37 | 96 | 105 |
| 520 | $1.27\pm0.02$ | 0.24 | 80 | 89 |

Caveats: $\chi^2_{\rm red}\approx 40$–$100$ everywhere (the MAD-based noise bar underestimates the true point scatter, and the response contains a fast onset followed by slow growth, which a single exponential cannot capture — hence "full-range" $\tau_r$ is inflated at low $a$ and collapses at $a\ge530$ where $P(l)$ is flat/decaying). The roll-modulation term was not needed: pump detuning from $2f_{\rm line}$ is only $\delta\approx0.2$–$0.5$ MHz and the total burst ($\le 0.184\ \mu$s) is far shorter than the roll period ($\gtrsim2\ \mu$s).

**Small amplitudes ($a\le 410$ mV): no trustworthy $\tau_r$.** Below $a_{\rm th}$ the line is undetected at any $l$ (SNR $<5$); just above threshold the fits are degenerate ($\tau_r$ hits the 50 ns bound with $100$–$200$ ns errors, or collapses to the 0.02 ns bound). Raw numbers ($l=60$ samples $=0.938$ ns vs $l=590=9.219$ ns; net window power, W):

| $a$ | $P(0.938\,{\rm ns})$ | SNR | $P(9.219\,{\rm ns})$ | SNR |
|---|---|---|---|---|
| 300 | $1.2\times10^{-7}$ (off-peak) | 2.7 | $8.6\times10^{-8}$ | 2.0 |
| 320 | $8.2\times10^{-8}$ | 1.4 | $5.2\times10^{-7}$ | 6.1 |
| 330 | $4.9\times10^{-8}$ | 1.0 | $4.9\times10^{-7}$ | 6.0 |
| 350 | $2.2\times10^{-7}$ | 3.5 | $7.2\times10^{-7}$ | 9.1 |
| 400 | $2.7\times10^{-7}$ | 2.8 | $3.9\times10^{-6}$ | 42 |

### 1.4 Implication for the mode damping scale

Converting every dynamic time scale extracted from A into the small-signal rate $1/(2\tau)$:

- threshold build-up: $\tau_b=2.15$ ns $\Rightarrow 1/(2\tau_b)=233$ MHz;
- onset-restricted ring-up at the weakest trusted drives ($a=420$–$430$ mV): $\tau_r=2.1$–$2.7$ ns $\Rightarrow 185$–$238$ MHz;
- the most conservative fits just above threshold ($a=380$–$410$, poor $\chi^2$): $\tau_r\approx4$–$10$ ns $\Rightarrow 50$–$125$ MHz;
- total-band (dilution-free) power at $a\ge450$ already saturates by $\tau_p\approx5$–$7$ ns ($P_{\rm tot}(9.2)/P_{\rm tot}(4.7)=1.06$–$1.15$).

**Conclusion: $\Gamma_{\rm eff\text{-}eq}=1/(2\tau_r)\sim 10^2$ MHz (working range 80–250 MHz, all central values $\ge 185$ MHz).** This is consistent with the two-tone delay scale ($\tau\approx2.4$ ns $\leftrightarrow\sim200$ MHz) and **inconsistent with the naive small-signal 8–10 MHz** ($\tau=50$–$60$ ns): with $\tau=50$ ns one would predict (i) no saturation of $a_{\rm th}$ across the whole $l$ range ($a_{\rm th}\propto(1-e^{-l/50{\rm ns}})^{-1/2}$ would keep rising beyond $a=1000$ mV for $\tau_p<20$ ns), (ii) $P(9.2)/P(3.4{\rm ns})\approx2.6$ at threshold-adjacent drive — observed only at the very weakest drives and already $1.1$–$1.9$ at $a=420$–$500$ mV, and (iii) a line FWHM broadened by the $2\tau_p$-period gate modulation everywhere — whereas the FWHM narrows from $\sim$16 MHz ($\tau_p\approx1$ ns) to $\sim2$ MHz ($\tau_p\ge6$ ns, limited by RBW/mode linewidth).

### 1.5 Pulse-train geometry signature

Because the 10-pulse train repeats back-to-back during the sweep, the mode sees a **continuous 50 %-duty gate**, for which the rate-equation threshold ($\oint G\,dt$ per unit time vs. loss) is independent of the pulse granularity — so the observed $a_{\rm th}(l)$ rise at short $\tau_p$ is **not** a build-up threshold effect. Two measured facts separate the effects: (i) the line FWHM at fixed $a$ scales with the emission comb envelope (16 MHz at $\tau_p=0.94$ ns $\leftrightarrow 1/T_{\rm burst}=53$ MHz; 2 MHz at $\tau_p\ge8$ ns), i.e. short pulses spread idler power over many comb lines and only a fraction lands in any fixed detection window — this produces most of the short-$\tau_p$ threshold elevation; (ii) nevertheless the *total*-band power (immune to in-span spreading at $\tau_p\gtrsim0.5$ ns, where the comb spacing $<75$ MHz) still grows with $\tau_p$ and saturates on a few-ns scale ($P_{\rm tot}(9.2{\rm ns})/P_{\rm tot}(4.7{\rm ns})=1.06$ at $a=500$ vs $1.78$ at $a=400$), so a genuine diminishing-returns buildup on the few-ns scale coexists with the dilution. No signature of per-10-pulse (burst-level) memory effects appeared beyond this — consistent with $\tau_{\rm decay}\lesssim$ few ns $\ll$ burst durations.

Line position when detected: $2036.0\pm0.4$ MHz ($=f/2$ within the RBW grid); sub-threshold "peaks" are noise spikes scattered over the whole span ($\pm20$ MHz, FWHM $<1.5$ bins), which is what the persistence + 5$\sigma$ criteria remove.

---

## 2. Dataset B — `phase gating/delay_f1_JULY_varyf1a1.h5`: $a_1\times f_1$ sweep, two tones, tone-2-first

### 2.1 Geometry (from the embedded `script_content`, computed from samples)

Record $=96000$ samples $=1500$ ns at 64 GSa/s:

- **Tone 2** at $f_2=4150$ MHz, $a_2=588$ mV: `window2[35200:57600]` $= 550\,{\rm ns}\rightarrow900\,$ns, length $22400$ samples $=350$ ns (tukey $\alpha=64/22400$, effectively rectangular). **It starts first** and **ends exactly at sample 57600**.
- **Tone 1** at $f_1\in[4020,4119]$ MHz, $a_1\in[256,328]$ mV: `window1[57600:64000]` $=900\,{\rm ns}\rightarrow915.625$ ns, length $1000$ samples $=15.625$ ns. It begins exactly where tone 2 ends.
- **Overlap = 0 samples; gap = 0 samples.** The two tones are exactly adjacent (edge-to-edge meeting point at $t=900$ ns); tone 2 comes first.
- Both carriers use the common time origin $t_0=550$ ns $=(64000-16000-12800)/64\,$GHz (start of tone 2) with phases $[-\varphi/2,\,+\varphi/2]$ on [tone 2, tone 1]. Therefore the stored `phase` $\varphi$ is, to within the drift term $2\pi(f_1-f_2)\cdot350\,{\rm ns}\le16.4^\circ$ over the sweep, **the relative carrier phase at the meeting point $t=900$ ns** — physically the relative parametric phase at the hand-off from the strong tone to the probe tone (which inherits the gated idler phase).
- SA: center $= f_2/2 = 2075$ MHz exactly (the idler line sits at band center by construction), span 15 MHz, 1001 pts (0.015 MHz/bin), RBW 1 MHz, avg 1 (single repeat), sweep 0.25 s.

Key inventory: 31399 keys, $a_1$: 31 values (256, 264, 268…328 step 2 — irregular low end), $f_1$: 95 values ($4020\ldots4119$, mixed 1/2 MHz steps), phases on a uniform grid: 794 blocks with the full 32-point grid ($\pi/16$) and 364 with the 16-point grid ($\pi/8$), 10 partial blocks (9–29 phases) — 167 keys in non-uniform/near-empty blocks were skipped. Grid is **sparse**: 1158 of $31\times95=2945$ $(a_1,f_1)$ cells populated (71–80 % for the 32-phase cell set). $\Delta=4150-f_1$ coverage $31$–$130$ MHz.

### 2.2 Harmonic decomposition and maps

Per $(a_1,f_1)$ block: $I(\varphi_j)=$ bin-summed baseline-subtracted power in $|f-2075|\le1.75$ MHz (single repeat); then
$c_k=\frac{2}{N}\sum_j I(\varphi_j)e^{-ik\varphi_j}$, $A_0=\langle I\rangle$, $A_1=|c_1|$, $A_2=|c_2|$, $\varphi_0=\arg c_1$, $\pi$-jump indicator $\arg c_2$ (jumped $\equiv |\arg c_2|>\pi/2$).

Saved npy maps (rows $a_1$ 256…328, cols $f_1$ 4020…4119, NaN = no data/corrupt):

`B_A1.npy, B_A2.npy, B_R21.npy (A2/A1), B_A0.npy, B_ARG2.npy, B_PHI0.npy, B_CORRUPT.npy, B_NPH.npy, B_a1s.npy, B_f1s.npy`; per-block $I(\varphi)$ and timestamps in `B_blocks.pkl`; figures `fig_B_maps.png`, `fig_B_delta.png`, `fig_A_summary.png`.

Summary statistics (clean blocks):

| map | max | location | min | median |
|---|---|---|---|---|
| $A_1$ [W·bin-sum] | $7.49\times10^{-5}$ | $a_1{=}328$, $f_1{=}4042$ ($\Delta=108$) | $2.5\times10^{-7}$ | $1.4\times10^{-5}$ |
| $A_2$ | $4.88\times10^{-5}$ | $a_1{=}328$, $f_1{=}4108$ ($\Delta=42$) | $9.0\times10^{-8}$ | $9.0\times10^{-6}$ |
| $A_2/A_1$ | 23.8 (A1-near-zero point, $a_1{=}326$, $\Delta=122$)  | — | $7.3\times10^{-3}$ | 0.55 |
| $A_0$ | $6.04\times10^{-4}$ | $\Delta=130$, $a_1=256$ | $2.18\times10^{-4}$ | $3.7\times10^{-4}$ |
| $\arg c_2$ | — | bimodal, median $-17^\circ$; 55.8 % of points have $|\arg c_2|>\pi/2$ | | |
| $\varphi_0=\arg c_1$ | median $+53^\circ$; wanders over the whole circle between cells | | | |

### 2.3 $\Delta$-dependence and the $c_2$ $\pi$-jump line

$\Delta$-binned averages over $a_1$ (clean blocks):

| $\Delta$ (MHz) | $n$ | $\langle A_1\rangle$ | $\mathrm{med}\,A_2/A_1$ | frac $\pi$-jumped | med $\arg c_2$ |
|---|---|---|---|---|---|
| 30–40 | 57 | $2.6\times10^{-5}$ | 0.71 | 0.79 | $+89^\circ$ |
| 40–50 | 71 | $2.9\times10^{-5}$ | 0.86 | 0.96 | $-145^\circ$ |
| 50–60 | 71 | $2.2\times10^{-5}$ | 0.64 | 0.92 | $-113^\circ$ |
| 60–70 | 186 | $2.1\times10^{-5}$ | 0.51 | 0.75 | $+69^\circ$ |
| 70–80 | 180 | $3.3\times10^{-5}$ | 0.39 | 0.60 | $-78^\circ$ |
| 80–90 | 179 | $2.4\times10^{-5}$ | 0.44 | 0.11 | $-16^\circ$ |
| 90–100 | 180 | $1.3\times10^{-5}$ | 0.68 | 0.24 | $+3^\circ$ |
| 100–110 | 75 | $2.3\times10^{-5}$ | 0.47 | 0.51 | $-8^\circ$ |
| 110–130 | 157 | $1.0\times10^{-5}$ | 0.76 | 0.76 | $-106^\circ$ |

- **$A_1$ does not peak at $\Delta=50$–$58$** here; its maximum band is $\Delta=70$–$80$ MHz and the minimum is at $\Delta=90$–$100$ (and 110–130). At the highest probe amplitude ($a_1=328$) $A_1$ peaks at $\Delta=108$; at $a_1=320/300$ at $\Delta\simeq73$–$78$; at low $a_1$ (256–288) the $A_1(\Delta)$ maximum moves to $\Delta\lesssim50$. No single low-$\Delta$ memory-maximum band is reproduced in this (different-geometry) dataset.
- **$c_2$ $\pi$-jump structure: reproduced qualitatively, shifted slightly.** At $a_1\ge304$ the wrapped $\arg c_2(\Delta)$ rows are clean and show exactly two sign flips of $\cos(\arg c_2)\,A_2$: at $\Delta\approx76\pm2$ MHz (individual rows: 74, 76, 76, 76, 78, 78, 78, 78) and a second at $\Delta\approx103\pm2$ MHz; between them lies an "unjumped island" ($\arg c_2\approx0$, $\Delta\approx80$–$100$). The lower edge at $f_1\simeq4072$–$4076$ brackets the known $\Delta=74.5$ MHz line (there measured at $\Sigma\approx4106$; here $\Sigma=(f_1+f_2)/2\approx4113$ at $f_1=4075$), i.e. reproduces within $\sim2$–$4$ MHz. At $f_1=4075/4076$ specifically, 72–78 % of points are in the jumped state vs 44–56 % at $f_1=4073/4074$ — consistent with a flip immediately below/at these $f_1$, but the lower-$a_1$ rows are too noisy to sharpen it (single repeat, RBW 1 MHz).
- The $\varphi_0=\arg c_1$ map is smooth along $a_1$ within a column but jumps between $f_1$ columns — consistent with the tone-2-gated idler phase plus an $f_1$-dependent propagation phase.

### 2.4 Data hygiene (baseline recalculation steps)

The baseline cache is invalidated every $>90$ s; the executed loop resets it only at session boundaries, so recalcs can land inside a 32-phase block (block $\approx18$ s; run span 12.6 h over 2 days $\Rightarrow O(500)$ recalcs, $\sim20$ % of blocks in principle exposed). Screen: per-phase off-resonance median power ($|f-2075|>2.5$ MHz) ordered by timestamp, best contiguous-prefix two-level split; flag if $|\Delta{\rm level}|>6\sigma_{\rm pooled}$.

- **2 of 1158 blocks flagged and excluded**: $(a_1,f_1)=(256,4062)$ and $(272,4088)$; 3 at $>5\sigma$, 7 at $>4\sigma$ ($0.6$ %). The two excluded blocks are removed from all harmonic maps above; no harmonic-fit results in this report use them. The floor is quiet enough that the large number of recalcs leaves the off-line level unchanged to within the per-phase scatter, so residual corruption is limited to these blocks.
- Phase sub-sampling: blocks with $<9$ usable phases or non-uniform grids (167 keys, 10 blocks) were excluded rather than fit.

---

## 3. Files

- `analysis_3_and_july.md` (this report), `fig_A_summary.png`, `fig_B_maps.png`, `fig_B_delta.png`
- A: `A_Pnet.npy` (primary $P_{f/2}(l,a)$ map), `A_Pwin.npy`, `A_Nseg.npy`, `A_Pkpos.npy`, `A_fwhm.npy`, `A_Bg.npy`, `A_noise.npy`, `A_Snr.npy`, `A_Ptot.npy`, `A_ath_naive.npy`, `A_ath_persist.npy`, `A_amps.npy`, `A_ls.npy`, `A_freqs.npy`, `A_ringup_fits.pkl`
- B: `B_A1.npy`, `B_A2.npy`, `B_R21.npy`, `B_A0.npy`, `B_ARG2.npy`, `B_PHI0.npy`, `B_CORRUPT.npy`, `B_NPH.npy`, `B_a1s.npy`, `B_f1s.npy`, `B_blocks.pkl`
- scripts: `an_a_*.py`, `an_b_*.py`, `an_figs.py`, `an_tables.py`
