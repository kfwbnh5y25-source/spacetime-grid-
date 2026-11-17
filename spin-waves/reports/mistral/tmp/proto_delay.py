import h5py, numpy as np
from collections import defaultdict
p = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/delay_extract/delay_f1.h5"
F1 = 4072.0
rows = []
with h5py.File(p, "r") as f:
    for k in f.keys():
        a, b, c, d = k.split("_")
        if float(b) != F1 or d != "1":
            continue
        tr = f[k]["trace"][:][0]
        rows.append((float(a), float(c), float(tr[1].astype(np.float64).mean())))
G = defaultdict(dict)
for delay, ph, I in rows:
    G[delay][ph] = I
def fit(ph, I):
    X = np.column_stack([np.ones_like(ph), np.cos(ph), np.sin(ph), np.cos(2*ph), np.sin(2*ph)])
    c, *_ = np.linalg.lstsq(X, I, rcond=None)
    pred = X @ c
    r2 = 1 - ((I-pred)**2).sum()/max(((I-I.mean())**2).sum(), 1e-30)
    return np.hypot(c[1], c[2]), np.degrees(np.arctan2(-c[2], c[1])), c[0], r2
ds = sorted(G)
print("delay_ns  A1_nW  phi0_deg  A0_nW  R2   n")
res = []
for d in ds:
    ph = np.array(sorted(G[d])); I = np.array([G[d][q] for q in ph])
    A1, phi0, A0, r2 = fit(ph, I)
    res.append((d, A1, phi0, A0, r2, len(ph)))
    print(f"{d:8.3f} {A1*1e9:8.2f} {phi0:8.1f} {A0*1e9:8.1f} {r2:6.3f} {len(ph)}")
# exponential decay fit on positive delays
res_p = [(d, A1, ph0) for (d, A1, ph0, A0, r2, n) in res if d >= 0 and n >= 12]
D = np.array([r[0] for r in res_p]); Y = np.array([r[1] for r in res_p]); P = np.array([r[2] for r in res_p])
sel = D <= 8.0
Df, Yf = D[sel], Y[sel]
# fit ln(A1) ~ -d/tau + c weighted, and also linear decay
from scipy.optimize import curve_fit
def expo(d, a, tau, c): return a*np.exp(-d/tau) + c
try:
    po, pc = curve_fit(expo, Df, Yf*1e9, p0=[Yf.max()*1e9, 3.0, 0.0], maxfev=20000)
    print("exp fit: A1_0=%.1f nW tau=%.3f ns offset=%.1f nW" % (po[0], po[1], po[2]))
except Exception as e:
    print("exp fit failed:", e)
sl, ic = np.polyfit(Df, P, 1)
print("phi0 slope: %.2f deg/ns -> %.2f MHz (2*pi*Delta). Expected Delta=%.1f MHz" % (sl, sl/(2*np.pi)*1000*1e0/1e0*360/360, 4142-F1))
print("phi0 slope (deg/ns):", sl, "-> Delta_equiv MHz =", sl/360*1e3)
# linear (non-exp) decay comparison
cl, il = np.polyfit(Df, Yf*1e9, 1)
print("linear decay slope: %.2f nW/ns, intercept %.1f nW -> zero crossing at %.2f ns" % (cl, il, -il/cl))
