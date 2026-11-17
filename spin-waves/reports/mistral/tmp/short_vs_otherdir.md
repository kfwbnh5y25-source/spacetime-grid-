# A (broad, pump1 fixed) vs B (other_dir, reciprocal) — CSI2 / Fourier / harmonics comparison

## Which plot is which
- **A** = `A_old_params/short/interaction_plots_broad/` — pump1 **fixed** at f1=4069 MHz, a1=324 mV; scanning pump2 (Sobol, a2≈220–588 mV, f2≈3960–4200 MHz), 4 phases {0, π/2, π, 3π/2}. Horizontal axis of all maps = scanned **pump2** frequency.
- **B** = `A_old_params/short/other_dir/interaction_plots_broad/` — reciprocal: pump (originally pump2) **fixed** at f=4130 MHz, a=266 mV; scanning the other tone (a≈280–350 mV, f≈4050–4090 MHz). Horizontal axis = scanned **pump1** frequency.
- Both read at f/2 on a swept spectrum analyzer, normalized to the only-f2 (only-tone) reference ("Relative response").

## 1. CSI2.png (absolute + relative response)

**A** (axis = scanned pump2 frequency):
- Absolute: response rises steeply with a2; yellow band (≈0 dB plateau) from ≈4090→4200 MHz for a2≳400 mV, signal floor ≈−40 dB.
- Relative (colorbar −18…+20 dB, both ends saturated at grid-cell scale):
  - **Hot bands (clipped at +20 dB):** (i) a narrow vertical band ≈**4058–4068 MHz, peak ≈4063–4065 MHz**, present for a2≈220–350 mV and fading above; (ii) a broad V/bowl band ≈**4090–4180 MHz** below a2≈350 mV (saturated red up to 4200 mV at lowest a2).
  - **Suppression (blue):** −8 to −18 dB lobes at ≈3980–4060 MHz (all amplitudes) and ≈4095–4130 MHz above ≈400 mV.
  - **Peak contrast:** ≈**+20 dB enhancement vs −16 to −18 dB suppression ⇒ ≈36–38 dB full pixel-to-pixel range** (colorbar-limited).
- **Diagonal/antidiagonal structure:** strong positive-slope diagonal ridge in the right panel (ptp) from (≈4090 MHz, 300 mV) to (≈4170 MHz, 450 mV); no clear antidiagonal. The hot-band pair (4065 vertical + 4130s bowl) is the "anti-crossing-like" structure around the fixed pump1=4069.

**B** (axis = scanned pump1 frequency):
- Absolute: nearly flat (−5…−15 dB) over the whole 4050–4090 / 280–350 mV window; no strong frequency structure.
- Relative: **+20 dB saturated (red) over almost the entire map above a≈325 mV**, and a soft diagonal threshold boundary (red extends down to ≈290 mV at 4080–4090 MHz, stays at 0 dB = white below ≈320 mV for 4050–4060 MHz). No blue (suppression) pixels at all in range.
- **Peak contrast:** ≈+20 dB (saturated) vs 0 dB below threshold; no measurable suppression lobe, so effective on/off contrast >20 dB but one-sided.
- **Diagonal/antidiagonal structure:** only the soft threshold boundary (positive diagonal: higher scanned f → enhancement starts at lower amplitude). No sharp diagonal/antidiagonal features.

## 2. Fourier_4f1_f2.png (4-phase components)

**A** (colorbar ±0.88):
- Dominant order: **cos(φ), c1**. ℜ(c1) shows a diagonal stripe (alternating negative ≈−0.3 and positive lobes) along ≈4100–4160 MHz, a2≈340–520 mV; ℑ(c1) similar, amplitude up to ≈±0.5.
- **c2 (cos 2φ)** is nonzero but weaker: ℜ(c2) up to ≈±0.3 along the same diagonal; ℑ(c2) ≈ 0 (blank).
- arg(c1)/π is structured along the diagonal (phase flips sign across the ridge); arg(c2)/π is ≈ uniformly +π (c2 mostly real-negative) — consistent with a mostly in-phase 2nd-harmonic component.
- Spatial structure of dominant c1: **diagonal ridge sweeping up in frequency with a2**, not vertical.

**B** (colorbar ±0.88):
- Dominant order: **stronger low-order pocket in c1 AND c2 concentrated at a≳320 mV, f≈4063–4080 MHz**: ℜ(c1) and ℑ(c1) both ≈−0.5…−0.7 (deep red/orange pocket), max |c1| = 1.000; ℜ(c2) has a similar pocket ≈−0.5 (peak |c2| = 1.000); ℑ(c2) ≈ 0.
- Spatial structure: a **localized pocket** (not a sweeping diagonal) near the top of the amplitude range, i.e. fixed-frequency, threshold-in-amplitude.
- arg maps: below threshold both phases are noise; above threshold arg(c1) ≈ −π*0.6…−0.8 pocket, arg(c2) ≈ −π (dark blue = −π).

## 3. phase_harmonics.csv (located in each interaction_plots_broad/ dir; columns amplitude, frequency, c0, c1, c2)

| | A (f1 fixed 4069, scan pump2) | B (fixed 4130, scan pump1) |
|---|---|---|
| rows | 10268 | 868 |
| mean \|c0\| | 72.34 | 7.53 |
| mean \|c1\| | 0.0936 | 0.0937 |
| mean \|c2\| | 0.0673 | 0.1227 |
| **\|c1\|/\|c2\| (mean)** | **1.39 (c1-dominant)** | **0.76 (c2-dominant)** |
| rows where \|c2\|>\|c1\| | 32% | 69% |
| max \|c1\| | 1.000 @ a2=424 mV, f2=4121 MHz | 1.000 @ a=349 mV, f=4073 MHz |
| max \|c2\| | 1.000 @ a2=357 mV, f2=4069 MHz | 1.000 @ a=349 mV, f=4067 MHz |
| mean \|c1\|/\|c0\| | 0.13% | 1.25% |
| mean \|c2\|/\|c0\| | 0.09% | 1.63% |
| max 2(\|c1\|+\|c2\|) (CSV units) | 2.94 @ a2=415 mV, f2=4118 MHz | 2.91 @ a=347 mV, f=4068 MHz |

So: **A is harmonic-order c1 (cos φ)-dominant; B is c2 (cos 2φ)-dominant** on average. Absolute maxima of c1/c2 are nearly identical (both clip at 1.000), but B reaches that clipping at only ≈345–350 mV whereas A needs ≈360–425 mV.

## 4. Vertical stripe near 4063 MHz in A (artifact / fixed-pump1 reference?)

- Visible in **A's ptp map**: faint vertical dashed stripe at ≈4063 MHz spanning the full amplitude range, values ≈3–10 interaction units (vs ridge peak ≈45–47 units), i.e. ~5–20× weaker than the real ridge.
- Also visible in **A's CSI2 third panel** (faint vertical streak ≈4060–4065 MHz) and as a faint dotted vertical line in **A's Fourier ℜ(c1) and ℜ(c2)** at ≈4063 MHz. It appears in *every* A map at a fixed frequency independent of the scanned pump2 axis → consistent with leakage/reference at the **fixed pump1 = 4069 MHz** (small offset from the plot pixel grid), not a spin-wave feature.
- **Not present in B** — B's axis is the scanned pump1 frequency, so the fixed tone (4130 MHz) lies outside B's frequency window; nothing analogous is seen. This supports the artifact interpretation.

## 5. A vs B: reciprocity / asymmetry

- **Landscape not reciprocal.** A (scan pump2 around fixed pump1=4069) shows a two-hot-band structure (narrow 4065 MHz band + broad 4130s bowl), a −16…−18 dB suppression lobe, and a diagonal c1 ridge tracking pump2 = f(K) with amplitude. B (scan pump1 around fixed 4130) shows one flat, one-sided, +20 dB-saturated enhancement above a≈325 mV with no suppression lobe and only a localized Fourier pocket.
- **Fixed-pump1 (A) gives the richer landscape; fixed-tone-at-4130 (B) gives the stronger *relative* gating at low amplitude**: B saturates >+20 dB by a≈325–330 mV and reaches max phase-modulation 2(|c1|+|c2|) ≈ 2.91 already at 347 mV; A reaches its max 2.94 only at a2=415 mV. In fractional terms B's phase modulation is ≈10× larger relative to |c0| (1.3–1.6% vs 0.1%).
- **Phase-order flip is the clearest asymmetry:** A → cos(φ) (c1) dominant (|c1|/|c2| = 1.39); B → cos(2φ) (c2) dominant (|c2|/|c1| ≈ 1.9 median, c2>c1 in 69% of rows). The phase-dependent gating is therefore qualitatively different depending on which tone is fixed, not a simple mirror image.
- Peak contrast numbers: A ≈ 36–38 dB pixel range (+20/−18, colorbar-limited) with max ptp ≈ 45–47 interaction units at ≈4135 MHz / 430 mV; B >20 dB one-sided (no blue side) with no resolvable ptp structure above noise in the right panel.
- Caveat: B's window (4050–4090 MHz, 280–350 mV) is much smaller than A's, so B samples only the lower-left corner of A's landscape; part of the asymmetry may be sampling-range, not physics.
