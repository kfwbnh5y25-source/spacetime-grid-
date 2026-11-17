# measure/ — acquisition scripts for setup PC (Paper A)

Transfer this whole folder to the setup PC and run there (needs `kNET`). All scripts
import the shared chassis `reports/mistral/hw_core_spectro.py` — transfer it too
(scripts reference it via the parent directory; or drop a copy next to them and
remove one `sys.path` level).

## Pre-flight (1 min, no hardware)

```
python meas1_longgap_memory.py --dry-run --out /tmp/t1.h5
python meas2_threshold_vs_phase.py --dry-run --out /tmp/t2.h5
python meas3_depth_repeat.py      --dry-run --out /tmp/t3.h5
```

Each writes a `.dryrun.json`; if that works, the code path is valid.

## What each measures, and which claim it serves

| Script | Measurement | Claim it verifies | Runtime (est) |
|---|---|---|---|
| `meas1_longgap_memory.py` | Long-gap delay scans, Δ = 50 / 58 MHz, t_gap 0→200 ns, 60 phases, 2 tones 327/588 mV, delay.h5 geometry | memory time τ_mem(Δ): theory predicts 51–61 ns; contrast still alive at 7 ns at low Δ. Also φ_eff winding at long gaps | ~3.5 h (`--trials 1` halves; `--config 0/1` runs one detuning, ~1.8 h) |
| `meas2_threshold_vs_phase.py` | a2 swept 268→312 mV (2 mV) across threshold × 16 phases, f1=4069/a1=324, f2=4111, long overlap burst | threshold formula h_th(Ψ) = [√(Γ²+δ²−W²sin²Ψ) − W cosΨ]/V — linear-in-cosΨ threshold shift; gives W from spectroscopy | ~25 min |
| `meas3_depth_repeat.py` | 60 phases × 8 trials, fresh baseline EVERY trial, two overlap operating points (4111/271 money point, 4135/430 strong zone) | ≥50% phase-gating depth with error bars; kills the baseline-step |c1| artifact (~9% of old blocks) | ~25–35 min |

## Phase conventions (IMPORTANT for analysis)

- meas1: phases referenced to a FIXED t0_ref (waveform start), phase on tone2 only —
  this reproduces the delay.h5 convention so φ₀ winds by 2πΔ·t_gap across the scan.
  The per-key meta carries `phi_eff_turn_on` (wrapped tone2 phase at its turn-on).
- meas2/meas3: split-phase (+φ/2, −φ/2), both referenced to the overlap start — the
  f1_4069_5 / fig3 convention; relative phase inside the overlap = φ exactly.
- Every h5 carries root attrs `phase_convention`, `waveform_digests` (full pulse
  specs + digests), `sa_settings`. Raw signal (`__SIG`) and raw baseline (`__BL`)
  are stored separately — never subtract in hardware.

## After the runs

Copy the three h5 files back. Analysis entries will be added under
`reports/mistral/tmp/` (see AGENT_STATE.md section "pending analysis hooks" for the
ready-to-run parser entry points).
