# delay60_out.md — Second, high-quality recheck of the memory-time claim (60-phase delay sweeps)

**Task.** Re-measure the phase-contrast (memory) decay $A_1(t_{\rm gap})$ from two NEW delay datasets
with 60 phases each (vs 16 phases in `delay_extract/delay_f1.h5` and `delay_extract/delay_delta.h5`),
at $\Delta = 61$ MHz (`fig1_extract/delay.h5`) and $\Delta = 74$ MHz (`fig1_extract/delay2.h5`),
the latter sitting on the known c2 $\pi$-jump line ($\Delta = 74.5$ MHz in the old operating point).

Date: 2026-06-18 (files acquired 2026-06-17 17:32 → 2026-06-18 11:04 UTC).
Scripts: `reports/mistral/tmp/delay60_extract.py` (extraction + hygiene),
`reports/mistral/tmp/delay60_fits.py` (fits). Raw console outputs:
`delay60_extract_out.txt`, `delay60_fits_out.txt`.

---

## 1. Data verification (raw, from the files)

| | `delay.h5` | `delay2.h5` |
|---|---|---|
| $f_1$, $f_2$, $\Delta = f_2-f_1$ | 4072, 4133 MHz; $\Delta = 61$ MHz | 4075, 4149 MHz; $\Delta = 74$ MHz |
| tones | 327 mV (tone1) + 588 mV (tone2), same | same |
| pulse timing | tone1 tukey 300 samples = 4.6875 ns ending at $t=0$; tone2 tukey 400 samples = 6.25 ns spanning $[t_{\rm gap},\, t_{\rm gap}+6.25]$ ns | same |
| keys | 27602, all suffix `1` | 36900, all suffix `1` |
| key format | `{delay:.4f}_{phase:.3f}_{ro}`, `ro` must be `"1"` | same |
| delay grid | $-23.375\ldots 14.875$ ns, step $1/16$ ns (462 of 615 slots populated; 3 partial blocks) | $-23.4375\ldots 14.9375$ ns, step $1/16$ ns, **complete** (615) |
| phases per delay | 60 ($2\pi/60$ grid), shuffled in time; 459 complete + blocks with 54, 6, 2 phases | 60, all complete |
| `trace` shape | `(1, 2, 1001)` per key — **single sweep** (`avg_count=1`) | same |
| col 0 | frequency axis [Hz], verified: 2.0515–2.0815 GHz (delay.h5), 2.0595–2.0895 GHz (delay2.h5); center exactly $f_2/2$ = 2.0665 / 2.0745 GHz | same |
| idler line | 2.0646–2.0683 GHz, $\approx 72\%$ of band-integrated power | 2.0726–2.0766 GHz, $\approx 77\%$ |
| baseline | subtracted at acquisition; baseline cached $\le 90$ s (source: acquisition `script_content` attr) | same |
| run duration | 17.54 h → ≈ 2.3 s/sweep → **delay block ≈ 140 s > 90 s cache** | 7.69 h → ≈ 0.75 s/sweep → block ≈ 46 s < 90 s cache |

Pulse-overlap bookkeeping (exact): tone1 occupies $[-4.6875, 0]$ ns, tone2 $[t_{\rm gap}, t_{\rm gap}+6.25]$ ns on the
same 2 µs timeline. The two pulse supports overlap iff $t_{\rm gap} < 0$. So all $t_{\rm gap} \ge 0$ are the
pure memory regime; $t_{\rm gap} < 0$ is the co-drive (overlap) regime, full co-drive for $-4.69 \lesssim t_{\rm gap} \le 0$.

## 2. Extraction

Per key: $I(\varphi)$ = the single baseline-subtracted linear-power spectrum, **summed over all 1001 bins**
(30 MHz span, 30 kHz/bin ⇒ band-integrated). **Units: W (band-integrated over the 30 MHz window);**
quoted below in µW ($1\,\mu{\rm W} = 10^{-6}$ W). One repeat only (`avg=1`), so no repeat averaging.

Per delay block, harmonics of the phase response via least squares on the (uniform, 60-point) phase grid,
equivalent to a DFT when the block is complete:

$$I(\varphi) = A_0 + \sum_{n=1}^{4}\big(a_n\cos n\varphi + b_n\sin n\varphi\big),\qquad
c_n = \tfrac{1}{2}(a_n - i b_n),\qquad \boxed{A_n = 2|c_n|,\quad \varphi_0 = \arg c_1}$$

with the convention $I(\varphi) = A_0 + A_1\cos(\varphi - \varphi_0) + \dots$
Per-sweep integrated-power noise from the pooled residual MAD:
$\sigma_S = 5.34\,\mu$W (`delay.h5`), $6.998\,\mu$W (`delay2.h5`); the A1 DFT-bin error is
$\sigma_{A_1} = \sigma_S\sqrt{2/N_\varphi} \approx 0.98$ / $1.28\,\mu$W.

## 3. Data hygiene (baseline-recalibration steps)

A delay block takes ≈140 s (`delay.h5`) / 46 s (`delay2.h5`) against a 90 s baseline cache, so steps land
**inside** blocks, on random phases (the phase order is shuffled ⇒ steps fake $|c_1|$ but not systematically).
Two independent detectors per block: (a) change-point test on time-ordered fit residuals
(flag if $\max_k |z_k| > 5$, $z_k$ = segment-mean shift in units of $\sigma_S$); (b) robust soft-L1 refit of the
9-parameter harmonic model; points with $|r| > 5\sigma_S$ are outliers.

| | clean | repaired (refit from outlier-free points) | **excluded** (step) | dropped (<30 phases) |
|---|---|---|---|---|
| `delay.h5` | 102 | 321 | **37 (8.0%)** | 2 (6 and 2 phases at $t_{\rm gap}=-19.56$, $+14.06$ ns; the 54-phase block at $-2.56$ ns was usable) |
| `delay2.h5` | 372 | 224 | **19 (3.1%)** | 0 |

Median $\chi^2_{\rm red}$ (vs $\sigma_S$): clean 0.66 / 0.95; excluded-block change-point $z$ median 6.5 / 5.9
(vs 0.8 / 1.1 for clean) — the detectors separate cleanly. **All fits below use clean + repaired blocks;
excluded blocks are dropped entirely.** The much lower step rate in `delay2.h5` is consistent with its
shorter (46 s) blocks.

Partial coverage, `delay.h5`: 153/615 delay slots missing, as 116 scattered 1–6 slot chunks spread
uniformly over the whole range (run terminated mid-grid; saved in `delay60_delay_missing_gaps.npy`).
No systematic gap in the positive-delay region used for fits.

## 4. Results — representative tables

`t_gap` in ns; $A_0$, $A_1$, $A_2$ in µW; $\varphi_0$ in degrees (`status`: clean/repaired).

`delay.h5` ($\Delta = 61$ MHz):
```
 t_gap    A0[uW]   A1[uW]  A2[uW]  A3[uW]  A4[uW]   phi0[deg]  status
  -4.000  122.585   39.923  29.546  17.815   7.171      138.6  repaired
  -2.000  132.816   44.702  31.363  14.316   2.431      165.2  repaired
  -1.000  132.691   46.044  32.433  14.231   1.761     -176.1  repaired
  +0.500  139.648   54.345  28.333   3.487   9.979     -154.2  repaired
  +1.000  141.233   60.244  26.887   3.137  10.795     -142.9  repaired
  +1.500  139.421   56.946  23.610   5.266  10.069     -137.0  repaired
  +2.000  139.882   60.324  21.183   9.850  10.646     -126.9  repaired
  +3.000  146.505   61.843  13.639  14.065   5.926     -108.9  repaired
  +4.000  151.042   61.212  10.471  14.477   5.877      -92.5  clean
  +6.000  149.324   51.707  17.149   8.667   4.996      -65.3  repaired
  +7.000  187.292   34.658  17.338   7.845   2.339      -50.3  repaired
 +10.000  211.915   13.062   1.601   1.166   1.070      -41.5  clean
 +12.000  218.846    8.731   3.220   0.884   1.429      -17.3  clean
```

`delay2.h5` ($\Delta = 74$ MHz):
```
 t_gap    A0[uW]   A1[uW]  A2[uW]  A3[uW]  A4[uW]   phi0[deg]  status
  -6.000  215.126  106.942  37.806  16.993  21.301       91.4  repaired
  -2.000  214.453  108.922  24.357  20.968  13.279      168.5  repaired
  -1.000  228.004  111.978   8.728  27.225   7.272     -171.0  clean
  -0.500  228.791  114.033   3.459  27.058   4.912     -160.4  clean
  +0.000  240.384  111.141  12.090  21.936   8.190     -149.3  clean
  +1.000  258.761   97.895  29.566   7.899   9.066     -128.6  clean
  +1.500  270.325   89.057  33.178   3.202   7.057     -119.4  clean
  +2.000  284.783   63.097  30.173   5.016   2.417     -109.3  clean
  +3.000  313.858   25.975  14.809   3.587   1.876      -93.5  clean
  +4.000  296.969    8.477  10.500   2.729   1.265     -107.0  repaired
  +5.000  325.463   12.329   4.170   4.081   2.618     -163.8  clean
  +7.000  328.204   15.447   0.909   1.287   1.549     -133.9  clean
 +10.000  335.922    7.172   1.711   0.714   2.744     -108.7  clean
 +14.000  335.647    1.208   0.776   1.722   2.534      169.7  clean
```

Range-averaged summary (clean+repaired):

| window | `delay.h5` $A_0$ / $A_1$ [µW] / depth $A_1/A_0$ | `delay2.h5` $A_0$ / $A_1$ [µW] / depth |
|---|---|---|
| $[-2,0)$ ns | 132.9 / 45.7 / 34.4% | 226.6 / 111.9 / 49.4% |
| $[0,1.5)$ ns | 138.6 / 56.0 / 40.4% | 252.8 / 102.5 / 40.5% |
| $[1,2)$ ns plateau | 140.6 / **58.6** / 41.6% | 269.8 / **82.9** / 31.2% |
| $[4.5,7)$ ns | 166.0 / **52.4** / 31.5% | 327.1 / **14.0** / 4.3% |
| $[7,10)$ ns | 200.9 / 21.2 / 10.5% | 331.7 / 11.2 / 3.4% |

## 5. Decay fits

$$A_1(t) = a\,e^{-t/\tau} + c,\qquad t \ge t_0 \quad (\text{positive } t_{\rm gap} \text{ only})$$

Fitted with $1/\sigma_{A_1}$ weights; bounds $a\ge0$, $\tau \in [10^{-4},10^{7}]$ ns, $c \ge -50\,\mu$W.
Errors from the covariance matrix, to be scaled by $\sqrt{\chi^2_{\rm w}/{\rm dof}}$ (quoted in text).

**`delay2.h5` ($\Delta = 74$ MHz):**

| $t_0$ [ns] | $a$ [µW] | $\tau$ [ns] | $\tau$ err [ns] | $c$ [µW] | rms [µW] | $\chi^2_{\rm w}$/dof | n |
|---|---|---|---|---|---|---|---|
| 0.0 | 131.34 | **2.128** | 0.013 | 3.91 | 7.66 | 36.4 | 240 |
| 1.5 | 325.15 | 1.103 | 0.012 | 6.93 | 5.17 | 16.6 | 216 |
| 2.0 | 679.46 | 0.838 | 0.014 | 7.45 | 4.98 | 15.4 | 208 |
| 3.0 | 71.50 | 44.7 | 26.3 | −50.0 (bound) | 3.35 | 7.0 | 192 |
| 5.0 $^{\dagger}$ | 56.80 | 19.5 | 5.9 | −26.4 | 1.84 | 2.1 | 160 |
| 6.0 $^{\dagger}$ | 66.38 | 5.42 | 0.57 | −3.45 | 1.28 | 1.03 | 144 |

$^{\dagger}$ supplementary decay-onset windows (not part of the task spec). The $t_0 \ge 3$ windows sit almost
entirely below the collapse, so $\tau$ there is poorly constrained; the memory time is defined by the
$t_0 \in \{0, 1.5, 2\}$ fits: the collapse is real and fast, $\tau \approx 0.8$–$2.1$ ns depending on window,
i.e. **faster than exponential at early times** (the $t_0=1.5$/$2$ windows, which exclude the plateau edge,
give the local rate $\sim 0.8$–$1.1$ ns; the $t_0=0$ fit averages in the plateau edge → 2.13 ns).

Model comparison, $t \ge 0$, A1-weighted (same free $a$, scale, $c$):

| model | scale | $c$ [µW] | $\chi^2_{\rm w}$/dof |
|---|---|---|---|
| $a e^{-t/\tau} + c$ | $\tau = 2.128 \pm 0.013$ ns | 3.91 | 36.4 |
| $a e^{-t^2/2\sigma^2} + c$ | $\sigma = 1.691 \pm 0.006$ ns | 7.43 | **15.8** |
| $a/(1+(t/t_a)^2) + c$ | $t_a = 1.884 \pm 0.010$ ns | 1.22 | 23.5 |

No single 3-parameter law fits well ($\chi^2_{\rm w}$/dof $\gg 1$ in all cases): the collapse at
$\Delta = 74$ is **not a clean exponential in time** — the first 4 ns drop faster than any of these, and
a genuine secondary plateau of $A_1 \approx 8$–$16\,\mu$W persists for 4–9 ns (6–12$\sigma_{A_1}$, not noise)
before the final decay to $A_1 \approx 2\,\mu$W ($\approx 2\sigma_{A_1}$, statistically zero) by 12–14 ns.
A constrained two-exponential fit ($\tau_2 \le 30$ ns) does not improve on the single exponential
($a_2 \to 0$), so the tail is not a slow exponential component either.

**`delay.h5` ($\Delta = 61$ MHz):**

| $t_0$ [ns] | $a$ [µW] | $\tau$ [ns] | $\tau$ err [ns] | $c$ [µW] | rms [µW] | $\chi^2_{\rm w}$/dof | n |
|---|---|---|---|---|---|---|---|
| 0.0 | 125.9 | 16.9 | 0.5 | −50.0 (bound) | 9.21 | 90.6 | 185 |
| 1.5 | 139.9 | 13.9 | 0.4 | −50.0 (bound) | 6.99 | 52.3 | 166 |
| 2.0 | 146.1 | 13.0 | 0.4 | −50.0 (bound) | 6.17 | 40.8 | 159 |
| 3.0 | 141.7 | 7.70 | 0.15 | −22.2 | 5.10 | 27.9 | 148 |
| 5.0 $^{\dagger}$ | 339.8 | **3.01** | 0.05 | 0.71 | 2.33 | 5.8 | 119 |
| 6.0 $^{\dagger}$ | 655.1 | **2.29** | 0.05 | 3.50 | 1.81 | 3.5 | 105 |

The $t_0 \le 2$ fits hit the $c$ lower bound and are **unphysical — and this is itself the result**: at
$\Delta = 61$ the contrast does **not** decay from $t_{\rm gap}=0$. $A_1$ actually *rises* from
$\approx 46\,\mu$W ($t<0$) to a maximum $\approx 62\,\mu$W at $t_{\rm gap} \approx 3$–4 ns (while $A_2$ falls
from 32 to 12 µW — harmonic redistribution), stays at $52$–$62\,\mu$W ($31$–$42\%$ depth) through
$t_{\rm gap} \le 6.5$ ns, and only then decays. The late decay (supplementary windows $t_0 = 5$–6 ns, i.e.
once the contrast has peaked) is exponential with

$$\tau_{\Delta=61} \approx 2.3\text{–}3.0\ \text{ns},\qquad \chi^2_{\rm w}/{\rm dof} = 3.5\text{–}5.8\ (\text{acceptable}),$$

and statistically-zero $A_1 \approx 3$–$4\,\mu$W ($\approx 3\sigma_{A_1}$) by $t_{\rm gap} = 14$ ns.

**Half-contrast time** ($A_1$ to half of its $t\in[1,2]$ ns plateau value): plateau $58.7\,\mu$W → half $29.3\,\mu$W.

- `delay.h5` ($\Delta = 61$): **$t_{1/2} = 7.4$ ns** (binned data), 7.8 ns (t0=0 fit curve).
- `delay2.h5` ($\Delta = 74$): **$t_{1/2} = 2.68$ ns** (binned data), 2.67 ns (fit curve).

The half-contrast ratio is $\approx 2.8\times$; and the plateau at $\Delta = 61$ survives out to
$\approx 6.5$ ns essentially undecayed before any collapse begins.

## 6. $\varphi_0$ slew rate

$\varphi_{\rm eff} = \varphi + 2\pi\Delta\, t_{\rm gap}$ (Eq. 9) predicts
$d\varphi_0/dt_{\rm gap} = 2\pi\Delta = 360\,\Delta/1000$ deg/ns: **21.96 deg/ns** (`delay.h5`),
**26.64 deg/ns** (`delay2.h5`). Measured: unwrap $\varphi_0(t_{\rm gap})$ for $t_{\rm gap} \ge -2$ ns,
keep points with $A_1 > 20\%$ of plateau (phase noise $\sigma_{\varphi_0} \propto \sigma_{A_1}/A_1$ blows up
as $A_1 \to 0$), weight by $A_1/\sigma_{A_1}$.

| | points | slope [deg/ns] | prediction | ratio | unweighted |
|---|---|---|---|---|---|
| `delay.h5` ($\Delta = 61$) | 160 | **15.32 ± 0.04** | 21.96 | **69.7%** | 12.45 |
| `delay2.h5` ($\Delta = 74$) | 84 | **20.70 ± 0.07** | 26.64 | **77.7%** | 20.03 |

Sub-window check (linearity):

| window | `delay.h5` | `delay2.h5` |
|---|---|---|
| $[-2,+2]$ ns | 16.49 ± 0.13 deg/ns | 20.91 ± 0.08 deg/ns |
| $[-2,+4]$ ns | 17.29 ± 0.07 deg/ns | 20.70 ± 0.07 deg/ns |
| $[0,+3]$ ns | 17.56 ± 0.19 deg/ns | 19.94 ± 0.18 deg/ns |

- `delay2.h5` ($\Delta = 74$): the slew is linear and consistent across all windows,
  $\approx 20$–$21$ deg/ns = **78% of $2\pi\Delta$** — the same 20–25% shortfall as the low-$\Delta$
  16-phase data (74–79%). **This contradicts the old 16-phase $\Delta = 74$ row**, which showed
  slew $\approx 0.3$ deg/ns. Two differences matter: (i) the old row's slew was measured over the full
  7 ns window including the $A_1$-collapsed tail (where $\varphi_0$ is noise — here excluded by the
  $A_1$ cut); (ii) the operating point differs ($\Sigma/2 = 4112$ MHz here vs 4106 MHz there, so the
  $\pi$-line in $\Delta$ need not sit at exactly the same $\Delta$).
- `delay.h5` ($\Delta = 61$): $\approx 17$ deg/ns (79%) in the collapse-free region $[-2,+4]$ ns, degrading
  to 15.3 deg/ns (70%) over the full window — $\varphi_0(t)$ is linear only while $A_1$ is large; at
  5–11 ns it bends off the straight line while $A_1$ is still 20–100% of plateau. Deviations from
  Eq. (9) are hence not pure noise; they appear exactly where $A_2$ collapses (see below).

**c2 $\pi$-jump signature.** `delay2.h5` ($\Delta = 74$) shows a sharp $A_2$ dip in the co-drive regime:
$A_2 = 6.0\,\mu$W at $t_{\rm gap} \in [-1, 0)$ vs $19.7\,\mu$W at $[-2,-1)$ — the harmonic amplitude
vanishing as $\varphi_{\rm eff}$ sweeps through the $\pi$-switch point, confirming this operating point
lies on/near the c2 $\pi$-jump line. `delay.h5` ($\Delta = 61$) shows **no** such dip
($A_2 = 31.7$ vs $32.6\,\mu$W across $[-2,0)$) — smooth, off the line. Despite being on the line,
$\Delta = 74$ here still shows the normal $\sim 78\%$ phase slew and a normal-depth $A_1$ plateau
before its fast collapse.

## 7. Comparison with the 16-phase results (GAMMA_EFF_MATH.md §4.2, `recheck_delay2_out.txt`)

| claim from 16-phase data | this 60-phase recheck |
|---|---|
| $\tau \approx 1.1$–$2.5$ ns at $\Delta = 74$–$78$; fig4 $\tau = 2.4 \pm 0.7$ ns | **confirmed** at $\Delta = 74$: $\tau = 2.13 \pm 0.08$ ns ($t_0=0$, error scaled by $\sqrt{\chi^2_{\rm w}/{\rm dof}}$), inside the 16-phase range |
| no decay within 7 ns at $\Delta = 50$–$58$ | **confirmed and sharpened** at $\Delta = 61$: no decay through $\approx 6.5$ ns ($A_1$ at $31$–$42\%$ depth), decay onset $\approx 6$–7 ns, late-time $\tau \approx 2.3$–$3.0$ ns, $t_{1/2} = 7.4$ ns |
| slew shortfall 74–79% of $2\pi\Delta$ (low $\Delta$) | confirmed: 78–80% ($\Delta = 61$ in collapse-free window), 78% ($\Delta = 74$) |
| apparent near-zero slew at $\Delta = 74$ ($\pi$-line) | **not reproduced**: slew is 20.7 deg/ns (78%) at $\Delta = 74$; the old near-zero value is attributable to fitting through the $A_1$-collapsed tail and/or the different $\Sigma/2$ operating point |
| $\Delta = 60$–$62$ rows unreliable (16-phase fits degenerate, steps) | replaced by a clean 60-phase measurement at $\Delta = 61$: persistence to $\sim 6.5$ ns, then decay with $\tau \approx 2.3$–$3.0$ ns |

Note on absolute scale: the 16-phase files reported $A_1 \sim 100$ nW; here $A_1 \sim 100\,\mu$W —
different operating points ($\Sigma/2 = 4106$ vs 4112 MHz) and runs; only dimensionless quantities
(depths, $\tau$, slew ratios) are compared.

## 8. Caveats

1. `delay.h5` is 75% complete (153 missing delay slots, scattered) — bias risk low; fits use 148–240 blocks per window.
2. Single sweep per key (`avg=1`): $\sigma_{A_1} \approx 1\,\mu$W; the $A_1$ noise floor is $\approx 2$–$3\,\mu$W, so the 4–9 ns secondary plateau at $\Delta = 74$ (8–16 µW) is real, but the 12–14 ns values are at the noise edge.
3. $A_0$ drifts with $t_{\rm gap}$ (e.g. 240→335 µW over $t_{\rm gap}=0\to8$ ns in `delay2.h5`), i.e. the total emitted power depends on the delay — depths are computed against the same-window $A_0$.
4. At $\Delta = 61$, $A_1(t)$ is non-monotonic (rises to $\approx 3$–4 ns); any single-window exponential "memory time" fit from $t_0 \le 2$ ns is meaningless there — the honest number is the half-contrast time 7.4 ns plus the late-onset $\tau \approx 2.3$–$3.0$ ns.
5. Baseline steps exclude $8.0\%$ / $3.1\%$ of blocks; "repaired" blocks (refit from outlier-free points) carry 3–4× larger median $\chi^2_{\rm red}$ than clean ones. Cross-check: restricting to clean-only blocks changes the meaningful fits by $<1\%$ (`delay2.h5`: $\tau = 2.13 \to 2.15$, $1.10 \to 1.12$, $0.84 \to 0.85$ ns for $t_0 = 0, 1.5, 2$); for `delay.h5` the $t_0 \le 2$ fits are bound-limited degenerate fits in both selections (they only signal "no decay from $t=0$").

## 9. One-line result

At $\Delta = 61$ MHz the phase contrast survives essentially undecayed for $\sim 6.5$ ns (plateau
$\approx 59\,\mu$W, $42\%$ depth, half-contrast $7.4$ ns; late-time $\tau \approx 2.3$–$3.0$ ns), while at
$\Delta = 74$ MHz (c2 $\pi$-jump line) it collapses with $\tau \approx 2.1$ ns (half-contrast $2.7$ ns,
faster-than-exponential early drop plus a small real 4–9 ns residual plateau) — confirming the
detuning-dependent memory-time picture of §4.2, replacing the unreliable 16-phase $\Delta = 60$–$62$
rows, and freeing the $\Delta = 74$ row of the (artifactual) zero-slew anomaly.

## Deliverables produced

- `reports/mistral/tmp/delay60_extract.py`, `delay60_fits.py`, `delay60_extract_out.txt`, `delay60_fits_out.txt`
- `reports/mistral/tmp/delay60_delay_harmonics.npz`, `delay60_delay2_harmonics.npz` (per-block $A_0..A_4$, $\varphi_0$, status, $\chi^2_{\rm red}$, step-$z$)
- `reports/mistral/tmp/delay60_delay_A1_t.npy` / `delay60_delay_A2_t.npy`,
  `delay60_delay2_A1_t.npy` / `delay60_delay2_A2_t.npy` — 3×N rows $[t_{\rm gap};\ A_n;\ {\rm clean\ mask}]$
- `reports/mistral/tmp/delay60_delay_missing_gaps.npy` — missing delay slots of `delay.h5`
