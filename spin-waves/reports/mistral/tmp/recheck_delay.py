"""Independent recheck of the fig4 memory-time claim (tau ~ 2.4 ns) from RAW delay data.

delay_f1.h5: keys {delay:.4f}_{f1}_{phase:.3f}_1, f2 = 4142 MHz fixed, 16 phases.
    I(phi) = tone2-on minus baseline? (baseline-subtracted linear power, as everywhere).
    Harmonics via DFT: c_n = (1/N) sum_phi I(phi) e^{-j n phi}; A1 = 2|c1| (phase-referencing
    cannot fake |c1|), phi0 = angle(c1) relative to I = A0 + A1 cos(phi - phi0).
delay_delta.h5: inventory only (key structure) + same analysis if layout matches.
"""
import re
import numpy as np
import h5py
from collections import defaultdict
from scipy.optimize import curve_fit

BASE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp/delay_extract"

HOT = (4066.0, 4078.0)  # per fig4: hot band f1


def inventory(path, nprint=6):
    with h5py.File(path, "r") as f:
        keys = sorted(f.keys())
        print(f"\n### {path.split('/')[-1]}: {len(keys)} keys; root attrs: {dict(f.attrs)}")
        for k in keys[:nprint] + keys[-2:]:
            g = f[k]
            print("  ", k, "->", list(g.keys()), getattr(g["trace"], "shape", None) if "trace" in g else "")
        return keys


def integrated_power(g):
    """Mean-over-repeats bin-summed baseline-subtracted power [mW]."""
    tr = g["trace"][:]
    spec = tr[:, 1, :].astype(np.float64)
    return float(spec.sum() / spec.shape[0])


def harmonics_I(path, want_f1=None):
    """DFT harmonics of I(phi) for every (delay, f1)."""
    with h5py.File(path, "r") as f:
        keys = sorted(f.keys())
        rec = defaultdict(dict)
        for k in keys:
            pr = k.split("_")
            if len(pr) != 4:
                continue
            d, f1, ph, ro = pr
            rec[(float(d), float(f1))][float(ph)] = k
        out = {}
        for (d, f1), pdict in rec.items():
            if want_f1 is not None and f1 not in want_f1 and not (HOT[0] <= f1 <= HOT[1]):
                continue
            phs = np.array(sorted(pdict))
            I = np.array([integrated_power(f[pdict[p]]) for p in phs])
            n = phs.size
            # assumes uniform phase grid starting at 0: DFT bin 1
            u = np.exp(-1j * phs)
            c1 = np.sum(I * u) / n
            c2 = np.sum(I * np.exp(-2j * phs)) / n
            out[(d, f1)] = dict(A0=float(np.mean(I)), A1=float(2 * abs(c1)),
                                phi0=float(np.angle(c1)), A2=float(2 * abs(c2)),
                                phi2=float(np.angle(c2) / 2), n_ph=n,
                                ph_lo=float(phs.min()), ph_hi=float(phs.max()))
        return out, sorted(set(d for d, _ in out)), sorted(set(f1 for _, f1 in out))


def fit_decay(dvals, a1vals, t_min=0.0):
    m = np.asarray(dvals) > t_min
    t = np.asarray(dvals)[m]
    y = np.asarray(a1vals)[m]
    if t.size < 4:
        return None
    def model(t, a, tau, c):
        return a * np.exp(-t / tau) + c
    p0 = [max(y[0] - y[-1], 1e-12), 3.0, y[-1]]
    try:
        popt, _ = curve_fit(model, t, y, p0=p0, maxfev=20000)
        r = y - model(t, *popt)
        return dict(tau=float(popt[1]), a=float(popt[0]), c=float(popt[2]),
                    rms=float(np.sqrt(np.mean(r**2))), t=t, y=y, yfit=model(t, *popt))
    except Exception as e:
        return dict(err=str(e), t=t, y=y)


def main():
    inv_keys = inventory(f"{BASE}/delay_delta.h5")
    inventory(f"{BASE}/delay_f1.h5", nprint=3)

    print("\n### delay_f1.h5: harmonic extraction (all f1, positive+negative delays)")
    H, delays, f1s = harmonics_I(f"{BASE}/delay_f1.h5")
    print("f1 values:", f1s)
    print("delays:", delays)
    hot = [f1 for f1 in f1s if HOT[0] <= f1 <= HOT[1]]
    print("hot-band f1:", hot)

    d_arr = np.array(delays)
    print("\nper-f1 A1(t_gap) [mW x 1e6 -> nW]:")
    print("t_gap[ns]:", " ".join(f"{d:+6.2f}" for d in delays))
    for f1 in f1s:
        a1 = np.array([H[(d, f1)]["A1"] * 1e6 if (d, f1) in H else np.nan for d in delays])
        if np.all(np.isnan(a1)):
            continue
        tag = "*" if HOT[0] <= f1 <= HOT[1] else " "
        print(f"f1={f1:7.1f}{tag}:", " ".join(f"{v:6.2f}" if np.isfinite(v) else "   nan" for v in a1))

    print("\nphi0(t_gap) [rad*100; only hot band]:")
    for f1 in hot:
        ph = np.array([H[(d, f1)]["phi0"] if (d, f1) in H else np.nan for d in delays])
        print(f"f1={f1:7.1f} :", " ".join(f"{v*100:6.1f}" if np.isfinite(v) else "   nan" for v in ph))

    print("\ndecay fits A1 = a*exp(-t/tau) + c on positive delays:")
    for f1 in f1s:
        a1 = [H[(d, f1)]["A1"] for d in delays if (d, f1) in H]
        dv = [d for d in delays if (d, f1) in H]
        r = fit_decay(dv, a1)
        if r and "tau" in r:
            print(f"f1={f1:7.1f}: tau={r['tau']*1e9:.2f} ns  a={r['a']*1e6:.1f} nW  c={r['c']*1e6:.1f} nW  rms={r['rms']*1e6:.2f} nW  n={len(r['t'])}")
        elif r:
            print(f"f1={f1:7.1f}: FIT FAILED: {r['err']}  y={np.array(r['y'])*1e6}")

    # also A2 harmonic decay (cross-check of the boss's c2 interest)
    print("\nA2(t_gap) hot band:")
    for f1 in hot:
        a2 = [H[(d, f1)]["A2"] * 1e6 if (d, f1) in H else np.nan for d in delays]
        print(f"f1={f1:7.1f}:", " ".join(f"{v:6.2f}" for v in a2))


if __name__ == "__main__":
    main()
