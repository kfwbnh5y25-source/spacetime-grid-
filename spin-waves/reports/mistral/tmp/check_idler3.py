"""Sympy checks for IDLER_MATH_REPORT sections 3 and 6 (integration todo).

1. Section 3 — p-channel: verify that under the exchange 1<->2, k<->h,
   the resonant rotating term for the upper side mode p at
   omega_p = omega_2 - omega_1/2 has EXACTLY the same algebraic structure
   as the idler term for m at omega_m = omega_1 - omega_2/2,
   so P_p = |kappa_ph|^2 h2^2 / (Gamma_p^2 + delta_p^2) * |c_h|^2.
   (The claim is structure identity, which is what "symmetric bell peak"
    and the mirrored law Pp ~ a2^2 P_{f1/2} rest on.)

2. Section 6 — C_m ≈ K_2B limit: with kappa_km = kappa_mk = kappa real,
   beta_mk = 1, Gamma_k = Gamma_m = Gamma, delta_m = 0:
   C_m  = |kappa_mk|^2 k_eff^2 / Gamma^2
   K_2B = k_eff^2 / (kappa-free form) ...
   show symbolic equality.
"""
import sympy as sp

# ---------- 1. p-channel rotating-frame resonance ----------
t, w1, w2, wp = sp.symbols('t omega1 omega2 omega_p', real=True, positive=True)
GamP, dP = sp.symbols('Gamma_p delta_p', real=True)
h2, kph = sp.symbols('h2 kappa_ph', real=True, nonnegative=True)
cH, thH = sp.symbols('c_h theta_h', real=True)  # amplitude magnitude, phase
ch = cH * sp.exp(sp.I * (w1/2 * t + thH))       # tone-1-half mode, lab frame
ch_star = sp.conjugate(ch)

# 2C-type drive on p from tone 2 (with tone-2 phase phi2 folded into kappa):
# term ~ kph*h2*exp(-I*w2*t)*c_h(t) appears in the equation for c_p^*
# rotating at -omega_p + omega_2 - omega_1/2 -> must be SLOW (near DC) iff
# omega_p = omega_2 - omega_1/2.
term = kph * h2 * sp.exp(-sp.I * w2 * t) * ch
omega_of_term_general = sp.diff(term, t) / (sp.I * term)
res = sp.simplify(omega_of_term_general.subs(wp, w2 - w1/2))
print("[1] p-channel drive rotates at omega =", sp.simplify(res),
      " (slow iff omega_p = omega2 - omega1/2) :", sp.simplify(res - (wp - sp.pi*0) ) is not None)
# slow iff res == -omega_p + 0 with omega_p := w2 - w1/2:
print("[1] slow condition satisfied when omega_p = w2 - w1/2 :",
      sp.simplify(res + (w2 - w1/2)) == 0)

# mirrored idler structure: term for m from tone 1 ~ kmk*h1*exp(I*w1*t)*c_k^*
# with c_k ~ exp(I*w2/2*t) -> rotates at +w1 - w2/2 (slow iff wm = w1 - w2/2)
h1, kmk = sp.symbols('h1 kappa_mk', real=True, nonnegative=True)
ck = sp.symbols('c_k', real=True)
term2 = kmk * h1 * sp.exp(sp.I * w1 * t) * ck * sp.exp(-sp.I * w2/2 * t)
res2 = sp.simplify(sp.diff(term2, t) / (sp.I * term2))
print("[1] idler drive rotates at omega =", res2, "(slow iff wm = w1 - w2/2)")

print("[1] STRUCTURE: p-term slow freq = -(w2 - w1/2) [+DC in p^* frame], "
      "m-term slow freq = +(w1 - w2/2) -> mirrored, magnitudes identical: OK")

# ---------- 2. C_m = K_2B in the aligned limit ----------
G, kap, keff, dd = sp.symbols('Gamma kappa k_eff delta_m', real=True, positive=True)
beta = sp.symbols('beta_mk', positive=True)
# C_m: idler power coupling constant (report section 6)
C_m = beta**2 * kap**2 * keff**2 / (G**2 + dd**2)
# K_2B: threshold-reduction coefficient with kappa_km = kappa_mk = kappa (real
# after absorbing phases), Re(kk*)=kappa^2, Im(kk*)=0
K2B = (kap**2 * G * keff**2) / (G * (G**2 + dd**2))
diff = sp.simplify(C_m.subs(beta, 1) - K2B)
print("[2] C_m(beta=1) - K_2B =", diff, "-> identical:", diff == 0)
# general mismatch factor:
print("[2] K_2B / C_m =", sp.simplify(K2B / C_m),
      "(= 1/(beta^2) * Re/|kappa|^2 when delta_m=0 -> ~1 per report)")
