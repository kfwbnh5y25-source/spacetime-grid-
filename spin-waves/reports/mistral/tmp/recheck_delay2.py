"""Recheck part 2: delay_delta.h5 (keys {delay:.4f}_{Delta}_{phase:.3f}_1, Delta = f2-f1 integer MHz,
delays -5..+7 ns, f1 = 4106 - Delta/2, f2 = 4106 + Delta/2) + cleaned delay_f1 tau fits.

Outputs:
  - per-f1 tau (fit window 2-13.3 ns and 1.5-8 ns) across hot band -> median +- IQR vs fig4's 2.4+-0.7 ns
  - per-delta (delay_delta): A1 vs t_gap, tau fit, and phi0 slew rate vs 2*pi*Delta prediction
"""
import numpy as np
import h5py
from collections import defaultdict
from scipy.optimize import curve_fit

BASE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/delay_extract"
HOT = (4066.0, 4078.0)


def integrated_power(g):
    tr = g["trace"][:]
    s = tr[:, 1, :].astype(np.float64)
    return float(s.sum() / s.shape[0])


def harmonics(path, keyfields):  # keyfields: ("delay","f1"|"delta")
    with h5py.File(path, "r") as f:
        rec = defaultdict(dict)
        for k in sorted(f.keys()):
            pr = k.split("_")
            if len(pr) != 4 or pr[3] != "1":
                continue
            a, b, ph = pr[0], pr[1], pr[2]
            rec[(float(a), float(b))][float(ph)] = k
        out = {}
        for (d, b), pdict in rec.items():
            phs = np.array(sorted(pdict))
            I = np.array([integrated_power(f[pdict[p]]) for p in phs])
            u = np.exp(-1j * phs)
            c1 = np.mean(I * u)
            out[(d, b)] = (float(np.mean(I)), float(2 * abs(c1)),
                           float(np.degrees(np.angle(c1))))
        return out, sorted(set(d for d, _ in out)), sorted(set(b for _, b in out))


def fit_tau(dv, av, t0):
    m = np.array(dv) >= t0
    t, y = np.array(dv)[m], np.array(av)[m]
    if t.size < 4:
        return None
    def model(t, a, tau, c):
        return a * np.exp(-t / tau) + c
    try:
        p, _ = curve_fit(model, t, y, p0=[max(y[0] - y[-1], 1e-13), 2.0, y[-1]], maxfev=20000)
        return float(p[1]), float(p[0]), float(p[2])
    except Exception:
        return None


print("=== delay_f1.h5 (ro=1 only): per-f1 tau ===")
H, delays, f1s = harmonics(f"{BASE}/delay_f1.h5", "f1")
taus, taus15 = [], []
for f1 in f1s:
    dv = [d for d in delays if (d, f1) in H]
    av = [H[(d, f1)][1] for d in dv]
    r8 = fit_tau(dv, av, 2.0)
    r15 = fit_tau(dv, av, 1.5)
    tag = "*" if HOT[0] <= f1 <= HOT[1] else " "
    s8 = f"tau={r8[0]:5.2f}ns a={r8[1]*1e6:5.1f} c={r8[2]*1e6:5.1f}" if r8 else "  no fit"
    s15 = f"tau15={r15[0]:5.2f}ns" if r15 else "no fit"
    print(f"f1={f1:7.1f}{tag} {s8}  {s15}")
    if r8 and HOT[0] <= f1 <= HOT[1]:
        taus.append(r8[0] * 1e9 / 1e9)  # delays already ns? values in ns units
        taus15.append(r15[0] if r15 else np.nan)
taus = np.array([t for t in taus if np.isfinite(t)])
print(f"hot-band tau (fit>=2 ns): median={np.median(taus):.2f} ns  IQR=[{np.percentile(taus,25):.2f},{np.percentile(taus,75):.2f}]  n={taus.size}")
taus15 = np.array([t for t in taus15 if np.isfinite(t)])
print(f"hot-band tau (fit>=1.5 ns, fig4 window): median={np.median(taus15):.2f} ns  IQR=[{np.percentile(taus15,25):.2f},{np.percentile(taus15,75):.2f}] n={taus15.size}")

print("\n=== delay_delta.h5 (ro=1 only): per-delta memory ===")
Hd, delays_d, deltas = harmonics(f"{BASE}/delay_delta.h5", "delta")
print("delays:", delays_d)
print("delta:", deltas)
for D in deltas:
    f1 = 4106 - D / 2
    dv = [d for d in delays_d if (d, D) in Hd]
    av = [Hd[(d, D)][1] for d in dv]
    a0s = [Hd[(d, D)][0] * 1e6 for d in dv]
    r = fit_tau(dv, av, 0.5)
    # phi0 slew: linear fit of wrapped phi0 vs t on |phichanges|<pi unwrapped sequence
    ph = np.array([Hd[(d, D)][2] for d in dv])
    o = np.argsort(dv)
    phu = np.unwrap(np.radians(ph[o])) if np.isfinite(ph[o]).all() else None
    slew = np.nan
    if phu is not None and len(dv) > 4:
        tt = np.array(dv)[o]
        m = tt >= -2.0
        if m.sum() > 3:
            slew = np.polyfit(tt[m], np.degrees(phu[m]), 1)[0]
    n_ph = len(dv)
    a_at0_1 = np.mean([a for d, a in zip(dv, av) if 0 <= d <= 1.5]) * 1e6 if any(0 <= d <= 1.5 for d in dv) else np.nan
    a_at4_7 = np.mean([a for d, a in zip(dv, av) if 4.5 <= d <= 7]) * 1e6 if any(4.5 <= d <= 7 for d in dv) else np.nan
    st = f"tau={r[0]:5.2f}ns" if r else "no-fit"
    print(f"Δ={D:3.0f}MHz (f1={f1:.0f}): A1[0-1.5ns]={a_at0_1:6.1f}nW A1[4.5-7]={a_at4_7:6.1f}nW {st} n_del={n_ph} slew={slew:7.1f}deg/ns (2πΔ={2*np.pi*D:6.1f})")
