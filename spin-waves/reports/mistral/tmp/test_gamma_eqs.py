"""Numerical self-test of GAMMA_EFF_MATH.md equations (2), (6)/(7), (8a). v2.

v1 BUG (fixed): the 5 us RK4 window could only detect growth >~3.7 MHz above threshold,
so the bisection measured a fixed-growth-rate contour ~1.47x the true threshold. v2:
  - (2)/(6)+(7): exact eigenvalues of the 2x2 linear system, bisection on max Re(eig) = 0
  - (8a): exact algebraic steady-state solve;
  - RK4 cross-check: bisection on the ASYMPTOTIC log-slope of |b| (threshold-unbiased).

Model (rate units 1/s, rotating frame):
    db/dt = -(G + j*delta) b + j*A conj(b) + s
"""
import numpy as np
from scipy.optimize import brentq

G = 8.0e6
W = 3.0e6


def eig_maxre(A, delta):
    M = np.array([[-G - 1j * delta, 1j * A],
                  [-1j * np.conj(A), -G + 1j * delta]], dtype=complex)
    return np.max(np.linalg.eigvals(M).real)


print("=== Eq (2): single drive threshold (eigenvalue bisection) ===")
ok2 = True
for delta in (0.0, 3.0e6, -3.0e6):
    pred = np.sqrt(G**2 + delta**2)
    sim = brentq(lambda a: eig_maxre(a, delta), 0.1 * G, 4.0 * G, xtol=1e-3)
    good = abs(sim - pred) / pred < 1e-6
    ok2 &= good
    print(f"delta={delta/1e6:+.0f} MHz: sim = {sim/1e6:.6f} MHz, pred = {pred/1e6:.6f} MHz  "
          f"{'PASS' if good else 'FAIL'}  (rel {abs(sim-pred)/pred:.2e})")

print("\n=== Eq (6)/(7): external + internal 4B pump (W = 3 MHz) ===")
ok6 = True
psi0 = 0.3
for delta in (0.0, 2.0e6):
    for psi_deg in (0, 45, 90, 135, 180, 210, 225, 270, 315):
        psi = np.deg2rad(psi_deg) - psi0
        def g(a, psi=psi, delta=delta):
            return eig_maxre(a + W * np.exp(-1j * psi), delta)
        a_sim = brentq(g, 0.1 * G, 4.0 * G, xtol=1e-3)
        pred = np.sqrt(G**2 + delta**2 - (W * np.sin(psi))**2) - W * np.cos(psi)
        good = abs(a_sim - pred) / max(abs(pred), 1e-12) < 1e-6
        ok6 &= good
        print(f"delta={delta/1e6:+.0f} MHz Psi={psi_deg:3d} deg: sim = {a_sim/1e6:.4f} MHz, "
              f"pred = {pred/1e6:.4f} MHz  {'PASS' if good else 'FAIL'}")

print("\n=== Eq (8a): sub-threshold coherent gain, delta = 0 ===")
ok8 = True


def steady_state(A, s):
    """Solve (-(G)) b + iA b* + s = 0 (delta=0) as a 2x2 real system, b = x + i y.

    iA b* = i(ar + i ai)(x - i y) = (-ai x + ar y) + i (ar x + ai y)
    -> real:  -(G + ai) x + ar y = -s_r
       imag:  ar x - (G - ai) y = -s_i
    """
    ar, ai = A.real, A.imag
    M = np.array([[-(G + ai), ar],
                  [ar, -(G - ai)]], dtype=float)
    x, y = np.linalg.solve(M, [-s.real, -s.imag])
    return complex(x, y)


sigma = 1e-6
# for real A>0: growing eigenvector direction (1,1)/sqrt2, decaying (1,-1)/sqrt2
forms = [("s = (1+1j)sigma/rt2 (amplified)", sigma * (1 + 1j) / np.sqrt(2)),
         ("s = (1-1j)sigma/rt2 (damped)", sigma * (1 - 1j) / np.sqrt(2)),
         ("s = sigma (real)", sigma + 0j)]
for label, s in forms:
    for Afac in (0.0, 0.3, 0.6, 0.85):
        A = Afac * G
        b = steady_state(A, s)
        p1 = (s.real + s.imag) / np.sqrt(2)  # amplified-quadrature projection
        p2 = (s.real - s.imag) / np.sqrt(2)  # deamplified-quadrature projection
        pred = p1**2 / (G - A) ** 2 + p2**2 / (G + A) ** 2
        good = abs(abs(b)**2 - pred) / pred < 1e-9
        ok8 &= good
        print(f"{label:32s} |A|/G={Afac:.2f}: sim |b|^2 = {abs(b)**2:.4e}, pred = {pred:.4e}  "
              f"{'PASS' if good else 'FAIL'}")

print("\n=== RK4 cross-check (asymptotic log-slope = 0 bisection) ===")
DT = 2.0e-9
NSTEP = 6000


def slope(A, delta):
    A = A if isinstance(A, complex) else A
    b = 1e-10 + 0j
    pts = []
    for n in range(NSTEP):
        f = lambda x: -(G + 1j * delta) * x + 1j * A * np.conj(x)
        k1 = f(b); k2 = f(b + 0.5 * DT * k1); k3 = f(b + 0.5 * DT * k2); k4 = f(b + DT * k3)
        b = b + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if abs(b) > 1e-2:
            pts.append((n * DT, np.log(abs(b))))
            break
        if n % 5 == 0:
            pts.append((n * DT, np.log(abs(b))))
    if len(pts) < 10:
        return None
    (t1, y1), (t2, y2) = pts[-10], pts[-1]
    return (y2 - y1) / (t2 - t1)


okr = True
for delta in (0.0, 3.0e6):
    pred = np.sqrt(G**2 + delta**2)
    a_sim = brentq(lambda a: slope(a, delta), 0.5 * G, min(3.99 * G, pred * 3.5), xtol=1e-3)
    good = abs(a_sim - pred) / pred < 1e-3
    okr &= good
    print(f"delta={delta/1e6:+.0f} MHz RK4: |A_th| = {a_sim/1e6:.4f} MHz vs pred {pred/1e6:.4f} MHz  "
          f"{'PASS' if good else 'FAIL'}")

print("\nSUMMARY:", "ALL PASS" if (ok2 and ok6 and ok8 and okr) else f"ok2={ok2} ok6={ok6} ok8={ok8} okRK4={okr}")
