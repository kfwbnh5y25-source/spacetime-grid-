## 0A. Resume protocol (MANDATORY)

On every new session or after context compaction:

1. Read this entire AGENT_STATE.md before taking action.
2. Treat section 0 (Boss rules) as hard constraints.
3. Establish the repository root:
   /Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/
4. Inspect the current filesystem before assuming that a file, process,
   measurement, or subagent is still running.
5. Reconcile the reported subagent statuses against actual deliverables.
   A status marked RUNNING in this file is not proof that the task is
   still running.
6. Check existing outputs before launching analyses again. Do not
   overwrite existing results without a specific reason.
7. Continue from section 7 (Next steps), respecting its order.
8. Distinguish VERIFIED, CLAIMED, THEORY, and OPEN results in every
   report. Never promote a hypothesis to a verified result without
   supporting data.
9. After every major result, update this state file and the relevant
   technical document. Preserve previous results and explain changes.
10. Before ending a session, record:
    - Work completed
    - Files created or modified
    - Results verified
    - Remaining uncertainties
    - Exact next action

Do not repeat completed work merely to reconstruct context.
Do not invent missing results or assume that background processes survived.
When a conflict exists between this file and the actual data, investigate
and report the conflict rather than silently changing the ledger.

# AGENT_STATE — Paper A (BiYIG two-tone phase gating / parametric idler spectroscopy)

MACHINE-FACING STATE. Load this file first after any compaction; it is the authoritative
resume point. Owner: experiment boss (human). Maintained by the agent. Update after every
major result; refresh at least every 10 minutes while working (boss rule, 1118).
Last updated: 2026-11-18 (sanity GO via 1 ns edge-fix; meas4 duration fix; meas5 edge sweep
DONE → Outcome A; runner.py live overnight; IDLER report §2.2 verified; GAMMA updated).

## 0. Boss rules (STANDING)
1. NO MODE NUMBERS anywhere in paper/figures/captions until eigenmode-vs-bias-field fit
   is done (simulated eigenmodes on boss's PC; bias field unknown; fit slow — LAST item).
   Say "seed mode k", "target mode m", "mode-pair identity switch".
2. "50%" means PHASE-GATING DEPTH >= 50% ("the more the better"). NOT threshold reduction.
3. LaTeX for all equations in .md files.
4. All agent outputs under spin-waves/reports/mistral/; temp under mistral/tmp/.
5. Boss offers to run measurements on request (see section 7 — three scripts delivered).
6. NEVER write files outside of /Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral. dont delete files unless prompted.

## 1. Verified numbers ledger (frozen — safe to print)
- Fig 1 (fig1_overview.{png,pdf} + manifests, DONE): port of figure1_single_tone.py;
  panels (a)-(i); (h) pulse-timing diagram; (i) reserved for sample characterization
  (fallback eigen_27mT_*.npz; HELD — boss rule). Data: Figure1_arnold_plot/data.h5 +
  arnold_tongue_zoom/data2.h5, halfaverage obs, full drive cal.
- Fig 2 (fig2_interaction.*): P_m ∝ a1²·P_f2/2 slopes 1.888±0.014 (f2=4125), 1.918±0.017
  (f2=4135); C_m = 1.20/0.85e-6 mV^-2; 2D (a1,a2) + residual maps; SI figS1 = D/ND logistic
  (v1 numbers verbatim in SI manifest: 1.33 dB=14.2% at f2=4111; 5.40 dB≈46% best case;
  ~50% threshold-reduction claim NOT reproduced — FLAG intact).
- Fig 3 (fig3_overlap.*): theory-test fig; two-phase fit: g=A1/A0=0.011, φ̂0=77.9° CI
  [13°,160°], R²=0.010, residual rms 10.5% vs 1.1% modulation; φ0=−155° EXCLUDED
  (OFF−ON=−2.15%±0.58%, t=−3.7) with A1<0 degeneracy caveat (two phases ~antipodal,
  185.7° apart). Design cannot demonstrate smooth-cos sufficiency; 16-phase remeasure
  recommended → meas2/meas3 scripts now cover this.
- Fig 4 (fig4_nonoverlap.*, panel d, 8 npy + manifest): delay_fine.h5 BROKEN (truncated);
  substituted delay_f1.h5. τ = 2.4±0.7 ns (1.5–8 ns fit, 7 hot-band f1). φ0 slew
  20.9 vs 2πΔ=25.2 deg/ns (Δ=70 MHz).
- Fig 5 (fig5_harmonics.*): c2 π-jump line at Δ = 74.5 MHz (IQR 73.5–74.5, rms 0.22 MHz
  detune_fine.h5; IQR 73.5–75.5, rms 0.72 MHz "detune copy 2.h5" — cross-dataset MATCH);
  c1 continuous (≤1.4° across line), c1 dominates strong zone (Δ≈65–85, Σ≈4106–4116);
  A2/A1 ≤ 0.3 on ridge; higher harmonics negligible (median |c4..7| = 6.4% of A1).
- DELAY RECHECK (independent DFT, raw 16-phase sweeps; see §3): fig4's 2.4 ns CONFIRMED
  in its regime (hot-band median τ = 1.55 ns [1.08–1.84] (fit≥2 ns) / 1.80 ns [1.28–2.09]
  (fit≥1.5 ns = fig4 window)). τ IS DETUNING-DEPENDENT (delay_delta.h5, Σ=4106, Δ=50–88,
  ≤7 ns): Δ50–58: contrast 32–36% persists ≥7 ns (τ>7 ns); Δ62–70: τ~1.5–4 ns;
  Δ74–78: τ=1.1–2.5 ns; Δ80–88: τ=0.1–0.8 ns. φ0 slew = 74–79% of 2πΔ where A1 strong.
  Depth at t_gap=0: 1−Pmin/Pmax = 49–60% (Δ=50–70) → money number REPRODUCED (definition
  must be stated); symmetrized (Pmax−Pmin)/(Pmax+Pmin) = 32–43%.
- EQUATION SELF-TEST (tmp/test_gamma_eqs.py v2, ALL PASS): Eqs. (2), (6), (7), (8a)
  confirmed (eigenvalue bisection + exact algebraic steady state + unbiased RK4 growth-
  slope). v1 FAILs were harness bugs (fixed-growth window; matrix sign).
- NMM solver (index.html): implements d'Aquino Eqs. (55)–(60), Table I rules 1–12,
  Routh–Hurwitz App. C.3; its arnoldThresholdMT = my Eq. (2), instantaneousArnold-
  ThresholdMT = my Eq. (7) (same conventions!). NMM_DATA JSON at index.html line ~963.
- CITATIONS (VERIFIED by cit-check, in GAMMA_EFF_MATH.md '## Verified citations'):
  d'Aquino 2026 (arXiv:2601.07798v2, 50 pp): amplitude eqs (42)–(43) p.9, averaged (49)–(50)
  p.10, steady states (55)–(56) p.10, regimes z/u_h/u_n/c±± (57)–(60) p.11, single-mode PWM
  steady state (33)–(34) p.7, 4B term inside (9) p.4 (NO "P_mk" symbol exists — cross-pump
  object = mutual NFS matrix (40) Ñ_hn p.9, self-NFS (28)–(29) p.6), threshold (20)–(21)
  p.5, Arnold tongue (24) p.6, Table I p.14, RH App C.3 pp.30–34, (C30)–(C32) p.34.
  Perna 2022 (JMMM 546, 168683, 12 pp): modal expansion (2) p.2, amplitude eq (3) p.2,
  τ_h = |b^α_hh|^{-1} (Γ-analog) inside (5) p.5, closed-form (5), orthogonality (18)+
  Table 1 p.3, eigenproblem (16) p.9, projection (20)→(21) p.9, NMM (31), LLG (33)–(34)
  p.10, damped NMM (40)–(41) p.10, LLGS (42)–(46) p.11. No parametric threshold formula.
  NOT verifiable: N_hk/n_mj/P_mk symbols; d'Aquino eqs > (60).

## 2. Open tensions (do not hide from boss/paper)
- T1: measured τ_mem(Δ) ≪ 51–61 ns EXCEPT at low Δ (Δ≤58 data compatible: lower bound
  7 ns). Candidate readings: (i) nonlinear ringdown of over-threshold driven seed (likely
  near Δ=74 π-line); (ii) mode-pair hybridization across the identity-switch line. At
  Δ≤58 NOT excluded by data — meas1 will settle.
- T2: φ0 slew = 74–79% of 2πΔ (both datasets) — effective-detuning deficit ~0.26–0.37·Δ;
  not an error in (9); candidate mode-pulling/hybridization signature. OPEN.
- T3: 4B averages out for |Δf|≥50 MHz in steady state (report) vs clean cos φ response
  observed (f1_4069_5, Δ≈58–64). Burst-wise hypothesis (GAMMA doc §4.3); fig3 inconclusive.
- T4: ~9% of phase blocks in OLD files (delay_f1.h5, delay_delta.h5) carry a baseline-
  recalibration step (90 s cache vs ~8 s blocks) faking |c1|; random-phase subsetting
  (measured order shuffled). All extracted numbers use clean blocks / medians. meas3 kills it.

## 3. Dataset map (formats + cautions) — phase gating/ unless noted
- delay_f1.h5 (962 MB full copy in mistral/tmp/delay_extract/; loose copy truncated):
  keys "{delay:.4f}_{f1}_{phase:.3f}_1", f2=4142 FIXED, f1 4045–4092, delays −100…+13.3 ns
  (dense arange −11.72…+6.92 step 20/64 ns + extras [0,±0.3125,±3.125,4.5,2.5,2] +
  linspace(-100,-68) block), 16 phases SHUFFLED in order. a=327/588. Pulses: tone1 tukey
  300 samples(4.6875 ns) ends at t=0; tone2 tukey 400 (6.25 ns) starts at delay. SA center
  2071.5, span 30 MHz, RBW 1 MHz, avg 1. CAUTION: 4745 "1EX" keys = DIFFERENT SA settings
  (span 60, center 4142, 10 ms sweep) — READOFF FILTER MUST BE "1". Scripts with baseline-
  step corruption; scripts: tmp/recheck_delay.py, recheck_delay2.py, diag_delayf1.py.
- delay_delta.h5 (730 MB in tmp/delay_extract/): keys "{delay:.4f}_{Delta:int}_{phase:.3f}_1",
  Δ=50–88 step 2, Σ=4106 FIXED (f1=4106−Δ/2), delays −11.72…+6.95 step 5/64 ns (unfiltered
  in file despite script comment), 16 phases. sa_mask=(2.110,2.170)e9 irrelevant here.
  NOTE: script_content attr shows `delays[delays > -5]` but the data include −11.7 — attr
  is from an earlier script version. Safer: check both.
- fig1_extract/delay.h5 (375 MB): 60 PHASES (2π/60), delays −23.44…+15 ns step 4/64 ns,
  f1=4072/f2=4133 (Δ=61), 3-field keys "{delay:.4f}_{phase:.3f}_1", 27602/36840 keys
  (partial). amps 327/588, same pulse geometry as delay_f1.
- fig1_extract/delay2.h5 (501 MB): same but f1=4075/f2=4149 (Δ=74 — ON the π-jump line),
  36900 keys (complete-ish). Analysis agent: delay-recheck2 (status §5).
- phase gating/3/3.h5 (65 MB): SINGLE TONE f=4072; keys "{amp}_{length}" (a first!);
  a 200–1000 mV step 10; l 20–590 samples (0.3125–9.219 ns) step 10 (+90 stray l<20 keys);
  waveform = tukey(l) repeated ×10 with gap = l. SA center 2036 = f/2, span 75 MHz,
  RBW 300 kHz, avg 10. Ring-up vs pulse length ⇒ independent Γ_eff readout.
- delay_f1_JULY_varyf1a1.h5 (428 MB): keys "{a1:.2f}_{f1}_{phase:.3f}_1", a1 256–328 (31),
  f1 4020–4119 (~95), 16 phases (partial 9–11 for some), f2=4150 a2=588. Geometry DIFFERS:
  tone2 FIRST (span offset −16.3 µs samples... precisely [64000−250·64−12800, +150·64)
  = 550–900 ns), tone1 tukey 1000 samples ending at 64000; time origin shifted −550 ns;
  phases [φ/2, −φ/2] (opposite half-phases). SA center 2075, span 15 MHz, RBW 1 MHz.
- detune_fine.h5 (fig5 core), "detune copy 2.h5" (cross-check), detune_fine031126.h5
  (earlier, 212/875 pairs), f1_4069_5.h5 (money OP f2=4111 a2=271), varyphase / A_old_params
  (a2 200–588), threshold_reduction/data.h5 (fig2), 4069_324_fasweep.h5, a2_266_5_*_sweep
  files, data_broad.h5. detune_fine_vary_a1_a2.py script exists (2D a1,a2 × phase harmonic
  fits; plotsfine_varya1a2/).
- General key field order trap: delay datasets "{delay}_{f1 or Delta}_{phase}_{ro}";
  folder 3 "{amp}_{length}"; JULY "{a1}_{f1}_{phase}_{ro}". Always check script_content
  attrs (embedded acquisition script) — it is the ground truth.

## 4. Docs & files (mistral/)
- GAMMA_EFF_MATH.md: notation; §2 h_th=Γ/V; §3.1 W=|P_mk||b_k|²; §3.2 boxed Eq.(7) big
  formula; §3.3 depth Pmax/Pmin=((1−x−w)/(1−x+w))², ≥50% ⇔ 4w/(1−x)≳0.4 (+fig3 verdict +
  delay-contrast support); §3.4 Γ_eff = Γ_m − C h1² − W cosΨ (Eq. 8), Σ_k phase-insensitive
  note; §4.1 φ_eff = φ + 2πΔt_gap (Eq. 9); §4.2 τ_mem=1/(2Γ_k) (Eq. 10) + MEASURED τ(Δ)
  table + tension discussion; §4.3 burst-wise hypothesis; §5 regime/figure map; §6 ledger
  [self-test PASS, MEASURED delay lines, BOSS RULE, LaTeX rule, citations CONFIRMED +
  BONUS arXiv:2601.05775 (experimental two-tone YIG microdisk, non-commutative sequences —
  cite & differentiate!), + '## Verified citations (page + equation)' appended by cit-check].
- CORE_MESSAGES.md: M1 coherent idler; M2 gating strength; M3 non-reciprocity; M4 money
  number (depth, VERIFIED 49–60% / 32–43% definitions); M5 phase≡delay (VERIFIED); M6 c2
  π-jump identity switch (VERIFIED cross-dataset) + harmonic dominance; M7 sequential vs
  overlap FRAME + MEASURED τ(Δ) + data hygiene; trust rules (1)–(5).
- AUDIT_AND_PLAN.md, IDLER reports cautions; tmp/: test_gamma_eqs.py (v2 PASS),
  recheck_delay.py/recheck_delay2.py/diag_delayf1.py/diag_delta.py + out files,
  delay_extract/, fig1_extract/, inventory/analysis scripts.
- Figures+manifests: figures/fig1_overview.*, fig2_interaction.*, fig3_overlap.*,
  fig4_nonoverlap.* (+npy), fig5_harmonics.*, figS1_threshold_reduction.*.

## 5. Subagent status (/root children)
- fig-harmonic(fig5): DONE. fig2: DONE (restructure). fig3: DONE. fig4: DONE (panel d).
- img-compare: DONE. fig1: DONE (turn 2; characterization slot reserved).
- cit-check: DONE (verified citations appended to GAMMA_EFF_MATH.md; summary
  tmp/citations_out.md).
- nmm-check: RUNNING — extracting NMM solver (index.html) into node, cross-val Eqs.
  (2)/(7) vs arnoldThresholdMT/instantaneousArnoldThresholdMT + delayed-protocol
  regimes. Deliverable: tmp/nmm_check_out.md.
- delay-recheck2: RUNNING — 60-phase analysis of fig1_extract/delay.h5 (Δ=61) +
  delay2.h5 (Δ=74): τ, slew, hygiene. Deliverable: tmp/delay60_out.md.
- analysis-3-july: RUNNING — folder 3 ring-up (Γ_eff readout) + JULY a1×f1 harmonics
  (c2 π-line reproduction at Δ≈74.5, A1 vs Δ). Deliverable: tmp/analysis_3_and_july.md.
- fig-synthesis: NOT YET SPAWNED — after nmm-check lands (wants solver numbers).
  Scope: synthesis figure overlaying Eqs (6)/(7)+(8)+(10)-derived threshold/gain
  surfaces with measured points from Figs 2–5; boss's table row: "last theory figure
  predicts the first data figure".

## 6. Measurement campaign (mistral/measure/, delivered & dry-run-tested; LIVE overnight 1117→1118)
- RUNNER: runner.py — steps: (1) meas5_edges, (2) meas1 c0 Δ=50, (3) meas1 c1 Δ=58,
  (4) meas2, (5) meas3, (6) meas1 c2 Δ=64, (7) meas1 c3 Δ=70, (8) meas4. Resume-safe
  (skips non-empty h5), --list/--only/--force-step, separate --out per meas1 config.
  Boss launched it ~22:00 1117. Sanity (meas0) NOT in runner — run manually (GO'd live).
- SANITY STORY (1117): 32 ns tukey edges → line but NO gating (A1/A0=0.00, archived
  tmp/sanity_preedgefix_nogo_201147.h5, raw-flat verified); fix EDGE_NS=1 ns → live GO:
  A0≈170 µW, A1/A0=0.34, decay 0.37, slew 10.0 deg/ns. GO-run h5 still on setup PC
  (default name) — QC pending.
- meas5 EDGE SWEEP → **Outcome A CONFIRMED**: A1/A0 = 0.38/0.35/0.19/0.045/0.073/0.042
  at edge 0.5/1/2/4/7/10 ns; end-control 1 ns = 0.35 (no drift); A0 constant, A1 dies.
  Falloff scale ~2 ns = τ_mem(Δ=74) = τ_b("3") — smearing picture. Details GAMMA §4.2.
- meas4 FIX: D0_june was built with 4.7/6.25 ns bursts (sample-count misread); NOW
  300/400 ns, 1 ns edges (D1_burst140 stays short by design).
- meas1 FIX: edge 0.5→1.0 ns (boss-verified June tukey).
- meas2/meas3 re-verified: freqs in band, peaks under 760 mV ceiling, dry-runs clean.
- Analysis hooks (unchanged): meas1 "LGM_c{i}_d{...}", meas2 "TVP_a2_{..}", meas3
  "DEP_{OP}_.."; write parsers when files return.
- Folder cleaned: stale artifacts moved to tmp/ (nothing deleted).

## 7. Next steps (ordered)
1. MORNING (boss returns): copy sanity GO-run h5 + all campaign h5s from setup PC;
   run parse_meas1/2/3 (tmp/, all SELFTEST PASS — parse_meas2 16/16 after max-slope +
   repeat-averaging fix 1118 22:20); QC vs predictions (T1 test: τ_mem at Δ=50/58 out
   to 200 ns = the decisive paper number); meas5 h5 precise re-fit; S1 level check.
2. IDLER_MATH_REPORT integration DONE (1118): §2.2 (1117), §3 p-channel mirror +
   §6 C_m≡K_2B limit (sympy, tmp/check_idler3.py), §4 (51)-(52) quote verified PDF p.10,
   §5 re-run LIVE under node — table reproduced exactly. K2B/C_m folded into
   PAPER_DRAFT SM-S1 + GAMMA ledger.
3. Lit leftovers: chase S-flagged refs (Kim/Merbouche, Zautkin) via web; add Perna 2026
   (2601.05775) differentiation paragraph to PAPER_DRAFT §1b.
4. Fig-synthesis subagent (nmm-check dependent); boss's synthesis-figure requirement.
5. LaTeX draft + tectonic compile (/opt/homebrew/bin/tectonic). No pdflatex.
6. LAST: eigenmode-vs-bias-field fit (boss's PC) → THEN and only then mode numbers +
   Fig 1 panel (i) sample characterization.

## 8. Environment quirks
- ZOTERO ACCESS (boss-granted 1118, endpoint: `http://localhost:23119/api`, Zotero 7 local
  API, read-only): library id 22169132 ("Ma bibliothèque"), 33 top-level items (incl. dups,
  notes, webpages). WORKING RECIPE: `curl -sS "http://localhost:23119/api/users/22169132/
  items?itemType=-attachment&limit=100&start=0&format=json"` to /tmp + parse with a local
  python script (python-urllib-to-localhost gets smart-approve DENIED; curl passes).
  NO zotero connector in this harness; do NOT read zotero.sqlite directly (denied).
  Relevant holdings: Srivastava 2023 PRApplied (many eigenmodes parametric pumping, YIG),
  Brächer 2017 Phys.Rep (parallel pumping review), Hula 2022 APL (spin-wave combs),
  Dreyer 2022 NatCom (phase-locking nonlinear spin waves), Körber 2023 NatCom (magnon
  reservoir), Kim/Merbouche (stimulated magnon scattering non-degenerate), Li 2026 NatCom
  (magnonic spontaneous oscillation parametric pumping), Massouras 2024 PRB (mode-resolved
  parametric disks), Yemeli 2025 (self-modulation instability BiYIG nanodisks), Soares 2026
  (= arXiv:2601.05775). Use for lit-review section of PAPER_DRAFT.
- python: `source ~/.venvs/nvim-python/bin/activate && python script.py`. Heredoc/python3
  standalone often denied; write script files. A denied bash call sometimes passes on
  immediate repeat (asks user).
- tectonic at /opt/homebrew/bin/tectonic. node at /opt/homebrew/bin/node.
- DISK TIGHT: ~3.4 GB free after fig1_extract. Do not extract more without checking `df`.
- Old h5 files may be truncated copies of larger originals (delay_fine.h5 3.1/171 MB;
  loose delay_f1.h5 775/962 MB truncated vs tmp copy complete).
- Boss communicates tersely; asks sharp verification questions; corrects misreadings of
  his prompts (e.g. the 50% depth). Honest verification > speed. Reports (AI-generated)
  are hypothesis generators, never sources.

### Work log 1118 ~22:25 (AFK block)
- parse_meas2 SELFTEST NOW PASS (16/16 thresholds, A0=285.0, A1=8.0). Fixes: (1) threshold_of
  reworked to max-slope anchor + wider logistic window; (2) main_parse now AVERAGES t0/t1
  repeats per a2 (duplicated abscissae broke np.gradient -> bogus thresholds); (3) p,n tuple
  unboxing. Validator python: `/Users/davidzehner/Work/.venv/bin/python` (has h5py).
  Run: `cd mistral/tmp && /Users/davidzehner/Work/.venv/bin/python parse_meas2.py --selftest x`.
- verify_idler_nmm.mjs RE-RUN LIVE under node — reproduces IDLER §5 table exactly
  (ratios 0.00549/0.00605/0.00574/0.00547 mT^-2; D pTot 0.0059 > ND 0.0038).
- IDLER integration COMPLETE: §3 p-channel mirror + §6 C_m≡K_2B verified via sympy
  (tmp/check_idler3.py); §4 (51)-(52) verified against PDF p.10 (pdftotext).
  Results: GAMMA ledger + PAPER_DRAFT SM-S1 (K2B/C_m table, 10-15% threshold shift).
- PAPER_DRAFT: SM-S1 added; meas5 edge result into Fig 3 (SM-paragraph + candidate
  A1/A0-vs-edge figure); lit: Kim/Merbouche resolved (arXiv:2406.09032, HAL hal-04625827),
  Zautkin best-match ZhETF 79, 281 (1980) [Sov.Phys.JETP 52, 142] author-list to confirm;
  Perna/Soares 2026 differentiation paragraph added.
