# Citations verification — one-page summary (2026 verification run)

Source PDFs (repo root), verified page-by-page with `pdftotext -f N -l N -layout`:
- d'Aquino, Perna, Merbouche, de Loubens, arXiv:2601.07798v2 (2026), 50 pp — `daquino.txt`
- Perna, Bruckner, Serpico, Suess, d'Aquino, J. Magn. Magn. Mater. 546 (2022) 168683, 12 pp — `perna.txt`

Both: printed page number == PDF page number (Perna page 1 unnumbered). Nothing was extrapolated;
anything not found is listed at bottom.

## d'Aquino 2026 (arXiv:2601.07798v2) — paper uses mode indices h, n and coefficients $N_h$, $N_{hn}$ (NOT $N_{hk}$/ $n_{mj}$ / $P_{mk}$)

| Item | Eq. | Page | What it is |
|---|---|---|---|
| Amplitude EOM (full NMM) | (9)+(10) | 4 | $\dot a_h = j\omega_h(b_h + \sum b_{hi}a_i + \sum c_{hij}a_ia_j + \sum \frac{d_{hijk}}{2}a_ia_ja_k)$; "third and fourth terms resemble those arising from three- and four-magnon scattering interactions" |
| Single-mode EOM w/ pump | (31) | 6 | $\dot a_h = j\omega_h[(1+N_h\|a_h\|^2-\delta h_a\|\varphi_h\|^2)a_h$ $-\delta h_a(\varphi_h,\varphi_h^*)a_h^*]-\lambda_h a_h$ |
| Two-modes EOM w/ both pumps | (42)-(43) | 9 | $\dot a_h = j\omega_h[1+N_{hh}\|a_h\|^2+N_{hn}\|a_n\|^2-\dots]a_h + \frac{j\omega_h}{2}\delta h_a(t)(\varphi_h,\varphi_h^*)a_h^*-\lambda_h a_h$ |
| Averaged two-modes | (49)-(50) | 10 | slow-variable system $(A_h,\Phi_h,A_n,\Phi_n)$ (this is what Table I/figs. 4-6 analyze) |
| Threshold, zero detuning | (20), (21) | 5 | $\delta h^{\rm thres,h} = 2\lambda_h/(\omega_h\|(\varphi_h,\varphi_h^*)\|)$; ellipticity form $=2\alpha\omega_h(1+e_h^2)/(1-e_h^2)$ |
| Threshold vs detuning (Arnold tongue) | (24), (25)-(26) | 6 | $\delta_h > \sqrt{(\omega_{\rm rf}-2\omega_h)^2/(\omega_h^2\|(\varphi_h,\varphi_h^*)\|^2) + (\delta h^{\rm thres})^2}$ — same structure as Eq. (2) in GAMMA_EFF_MATH.md, no factor-2 mismatch |
| Steady state, single-mode PWM | (33), (34) | 7 | $\|a_h\|^2_\pm = (\epsilon_h \pm s_h^\pm\epsilon_{\rm crit,h})/(2\omega_h N_h)$; under PWM only $\|a_h\|^2_+$ (34) admissible — index.html claim CONFIRMED |
| Four-magnon self coefficient | (27)-(29) | 6 | $N_h = \Re\tilde N_h$, $\tilde N_h = d_{h,h,-h,h} + d_{h,-h,h,h}/2$ |
| Four-magnon cross coefficient (≈ index.html "$P_{mk}$"/"$N_{hn}$") | (40)-(41) | 9 | $\tilde N_{hn} = d_{h,n,h,-n}+d_{h,h,n,-n}+d_{h,-n,h,n}$; $N_{hn}=\Re$, $\lambda_{hn}=\omega_h\Im$; (41): $\lambda_{hn}=\lambda_{nh}=0$, $N_{hn}=N_{nh}$ |
| Coupled steady-state amplitudes ("Eqs. (55)-(60)" claim) | (55)-(56) | 10 | $\|a_h\|^2 = L_{hn}[\|a_h\|^2_{\pm,\rm unc} - (N_{hn}/N_{hh})\|a_n\|^2_{\pm,\rm unc}]$, $L_{hn} = N_{hh}N_{nn}/(N_{hh}N_{nn}-N_{hn}^2)$ — NOT equations of motion |
| Regime definitions | (57)-(60) | 11 | $z$, $u_h$, $u_n$, $c_{\pm\pm}$ steady states |
| Table I rules 1-12 | Table I | 14 | rules 1-4 existence, 5-12 stability; rule 6: $c^{++}$ stable in $T^{++}$ near boundaries iff $\mathrm{sign}(N_{hh})\mathrm{sign}(N_{nn})(N_{hh}N_{nn}-N_{hn}^2)>0$ — index.html claim CONFIRMED |
| Routh-Hurwitz, Appendix C.3 | C 3, (C30)-(C32) | 30 (heading), 33 (RH text), 34 (C30-C32) | "Among various technical conditions, the Routh-Hurwitz criterion states that a necessary condition for the stability of P-mode $c^{++}$ is that all the coefficients of the polynomial are positive." Necessary-only caveat stated. — index.html claim CONFIRMED |
| Jacobian of two-modes P-modes | (C19)-(C23) | 30-31 | 4x4 linearized stability matrix |

## Perna 2022 (JMMM 546, 168683)

| Item | Eq. | Page | What it is |
|---|---|---|---|
| Modal expansion (projection) | (2) | 2 | $\delta\mathbf{m}_\perp = \sum_h a_h(t)\varphi_h(\mathbf{r})$; $a_h$ = "degrees of freedom" of the NMM |
| NMM amplitude equation | (3) | 2 | $\dot a_h = \frac{j\omega_h}{1+\alpha_G^2}(b_h + \sum b_{hi}a_i + \sum c_{hij}a_ia_j + \sum\frac{d}{2}a^3 + \sum\frac{e}{4}a^4 + \sum\frac{f}{4}a^5+\dots)$ |
| Orthogonality/normalization | (18), (19); Table 1 | 9 (3) | $(\varphi_h,\varphi_k)_{A^\perp_0} = \frac{1}{V_\Omega}\int_\Omega \varphi_h^*\cdot A^\perp_0\varphi_k dV = \delta_{hk}$; eigenproblem (16) p. 9 |
| Projection identity & free ODE | (20), (21) | 9 | $(\varphi_h,\partial_t\delta m_\perp) = \dot a_h$; $\dot a_h = j\omega_h a_h$ |
| Linear driven response + damping time | (4), (5) | 4, 5 | $\tau_h = \|b^\alpha_{hh}\|^{-1} = (\alpha_G\omega_h^2\|\varphi_h\|^2)^{-1}$ — the $\Gamma_m$-analog; diagonal approx $b_{hi}\approx\delta_{hi}(1+b^\alpha_{hi})$ |
| Conservative nonlinear NMM (Methods) | (31)-(32) | 10 | same structure as (3), coefficients $c^0, d^0, e^0$ |
| LLG with damping | (33)-(34) | 10 | LL eq ↔ Gilbert form with $\alpha_G$ |
| Damped nonlinear NMM | (40)-(41) | 10 | $(1+\alpha_G^2)\dot a_h = j\omega_h(\dots)$; $b_{hi}=\delta_{hi}+\alpha_G b^\alpha_{hi}$, $c_{hij}=c^0+\alpha_G c^\alpha$, $d$, $e$ likewise |
| LLGS (rf field + spin torque) | (42)-(43) | 11 | forcing decomposed via $\mathbf{f}_c = \mathbf{h}_{\rm rf}+\alpha_G\beta\mathbf{p}$, $\mathbf{f}_d = \alpha_G\mathbf{h}_{\rm rf}-\beta\mathbf{p}$ |
| Forcing projections + full NMM | (44)-(46) | 11 | mode-basis projections of $\mathbf{f}_c$, $\mathbf{f}_d$; full forced NMM |

## NOT verifiable / NOT present (do not cite)

- d'Aquino 2026: symbols "$P_{mk}$", "$N_{hk}$", "$n_{mj}$" — appear nowhere in the PDF. Closest
  real objects: $d_{hijk}$ (10), $N_h$ (28)-(29), $\tilde N_{hn}/N_{hn}$ (40)-(41).
- d'Aquino 2026: "Eqs. (55)-(60)" are steady states/regimes, NOT the amplitude equations of motion
  (those are (42)-(43) p. 9 and (49)-(50) p. 10).
- Perna 2022: no parametric threshold formula (use d'Aquino (20)/(24)); no $\Gamma_m$/Kp symbol —
  damping is $\alpha_G$-based, modal rate $\lambda_h = \alpha_G\omega_h^2\|\varphi_h\|^2$ (from Eq. (5), p. 5).
- d'Aquino 2026: no main-text equation beyond (60); appendices end at (A10)/(B6)/(C~32).

Full details with verbatim quotes: `reports/mistral/GAMMA_EFF_MATH.md`, section
"## Verified citations (page + equation)".
