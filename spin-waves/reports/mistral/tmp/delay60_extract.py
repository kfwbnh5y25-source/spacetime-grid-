"""Extraction of I(phi) harmonics from the 60-phase delay datasets (delay.h5, delay2.h5).

Key format: {delay:.4f}_{phase:.3f}_{readoff}; readoff must be '1'.
I(phi) = basis-subtracted linear power spectrum (single sweep, avg=1), summed
over all 1001 bins (30 MHz band centered at f2/2). Units: W (band-integrated).

Per (delay): harmonic fit over the 60 phases:
    I(phi) = A0 + Sum_n A_n cos(n phi - n*phi_n)  with  c_n = (1/N) Sum I e^{-i n phi},
    A_n = 2|c_n|, phi0 = angle(c1).

Hygiene (baseline-recal steps every ~90 s inside a ~140 s delay block):
  - change-point test on time-ordered residuals (steps land on random phases
    because the phase order is shuffled);
  - robust soft-L1 refit -> outlier count.
  status: clean / repaired (A1 refit from clean points) / excluded_step.
Blocks with < 30 phases: dropped.
Outputs: reports/mistral/tmp/delay60_{name}_harmonics.npz  +  A1/A2 .npy rows.
"""
import re
import numpy as np
import h5py
from datetime import datetime
from scipy.optimize import least_squares
from collections import defaultdict

BASE = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp"
FILES = [("delay", f"{BASE}/fig1_extract/delay.h5", 61.0),
         ("delay2", f"{BASE}/fig1_extract/delay2.h5", 74.0)]

HARM = [1, 2, 3, 4]


def parse_ts(x):
    if isinstance(x, (bytes, np.bytes_)):
        x = x.decode()
    return datetime.fromisoformat(str(x))


def design(phs, harm):
    cols = [np.ones_like(phs)]
    for n in harm:
        cols += [np.cos(n * phs), np.sin(n * phs)]
    return np.vstack(cols).T


def harmonic_fit(I, phs, harm=HARM):
    """OLS harmonic fit; returns coefs (A0, a1..a4, b1..b4), residual."""
    X = design(phs, harm)
    beta, *_ = np.linalg.lstsq(X, I, rcond=None)
    resid = I - X @ beta
    return beta, resid


def main():
    rng = np.random.default_rng(0)
    for name, path, delta_MHz in FILES:
        print(f"\n===== {name} =====")
        f = h5py.File(path, "r")
        all_keys = list(f.keys())
        rec = defaultdict(dict)  # delay -> phase -> (key, timestamp)
        for k in all_keys:
            pr = k.split("_")
            if len(pr) != 3 or pr[2] != "1":
                continue
            rec[round(float(pr[0]), 4)][float(pr[1])] = k
        delays = sorted(rec)
        nfull = 615
        print(f"delays present: {len(delays)}/{nfull} (missing {nfull - len(delays)}), "
              f"range [{delays[0]}, {delays[-1]}]")
        missing = np.array([d for d in (np.round(np.arange(-1500 / 64, 15, 4 / 64), 4))
                            if round(d, 4) not in [round(float(dd), 4) for dd in delays]])
        np.save(f"{BASE}/delay60_{name}_missing_gaps.npy", missing)
        for d in delays:
            if len(rec[d]) != 60:
                print(f"  partial block: t_gap={d:+.4f} ns has {len(rec[d])} phases")

        # --- pass 1: load all sweeps, OLS harmonics per block
        data = {}
        res_pool = []
        freq_axis = None
        for d in delays:
            phs_list, t_list, I_list, Sline_list = [], [], [], []
            for ph, k in rec[d].items():
                tr = f[k]["trace"][:].astype(np.float64)
                fr = tr[:, 0, :].ravel()
                spec = tr[:, 1, :].ravel()
                if freq_axis is None:
                    freq_axis = fr
                phs_list.append(ph)
                t_list.append(parse_ts(f[k]["timestamp"][0][0]))
                I_list.append(spec.sum())          # I(phi) [W], mean over 1 repeat
                # line-region sum (cross-check; line from mean spectrum)
            phs = np.array(phs_list)
            I = np.array(I_list)
            ts = np.array(t_list)
            order = np.argsort(ts)
            phs_t, ts_t = phs[order], ts[order]
            if len(phs) >= 30:
                beta, resid = harmonic_fit(I, phs)
                sigma_line = None
                data[d] = dict(phs=phs, phs_t=phs_t, ts_t=ts_t, I=I,
                               beta=beta, resid=resid, n=len(phs))
                res_pool.append(resid)
            else:
                data[d] = dict(phs=phs, phs_t=phs_t, ts_t=ts_t, I=I, beta=None,
                               resid=np.zeros(len(phs)), n=len(phs))
        # verify frequency axis
        c_pred = FILES[[x[0] for x in FILES].index(name)][2]
        fc = (freq_axis.min() + freq_axis.max()) / 2
        print(f"freq axis: {freq_axis.min()/1e9:.4f}-{freq_axis.max()/1e9:.4f} GHz, "
              f"center {fc/1e9:.4f} GHz (expect f2/2 = {(4072 if name=='delay' else 4075)/2/1e3:.4f})")

        # sigma_S: robust scatter of pooled residuals (per block OLS resid)
        rp = np.concatenate(res_pool)
        sigma_S = 1.4826 * np.median(np.abs(rp - np.median(rp)))
        print(f"sigma_S (per-sweep integrated-power noise, pooled resid MAD): {sigma_S*1e9:.2f} nW")

        # --- pass 2: hygiene per block
        results = {}
        n_excl = n_rep = 0
        for d in delays:
            dd = data[d]
            if dd["beta"] is None or dd["n"] < 30:
                results[d] = dict(status="low_n", beta=None, chi=None, step=None,
                                  n_out=dd["n"], n=dd["n"])
                continue
            beta, resid = dd["beta"], dd["resid"]
            # (a) change-point test on time-ordered residuals
            e = resid[np.argsort(dd["ts_t"])]
            n = len(e)
            cs = np.cumsum(e)
            tot = cs[-1]
            # split anywhere in [5, n-5]; z = shift of segment means in sigma units
            ks = np.arange(5, n - 4)
            s_after = (tot - cs[ks - 1]) / (n - ks)   # mean of e after split
            s_before = cs[ks - 1] / ks                # mean of e before split
            s = s_after - s_before                    # step size per split
            se = sigma_S * np.sqrt(1 / ks + 1 / (n - ks))
            z = np.abs(s) / se
            zmax = float(z.max())
            kbest = int(ks[np.argmax(z)])
            step_ns = (dd["ts_t"][kbest] - dd["ts_t"][0]).total_seconds()
            # (b) robust refit with soft-L1 loss
            X = design(dd["phs"], HARM)
            fs = max(sigma_S, 1e-30)
            sol = least_squares(lambda b: (dd["I"] - X @ b) / fs, beta,
                                loss="soft_l1", f_scale=1.0, max_nfev=2000)
            r_rob = dd["I"] - X @ sol.x
            n_out = int(np.sum(np.abs(r_rob) > 5 * fs))
            frac_out = n_out / dd["n"]
            if zmax > 5 or frac_out >= 0.25:
                status = "excluded_step"
            elif n_out > 0:
                status = "repaired"
            else:
                status = "clean"
            results[d] = dict(status=status, beta=(beta, sol.x), chi=None,
                              step_z=float(zmax), step_at=step_ns, n_out=n_out, n=dd["n"])

        for d in delays:
            st = results[d]["status"]
            if st == "excluded_step":
                n_excl += 1
            elif st == "repaired":
                n_rep += 1
        print(f"hygiene: clean={len(delays)-n_excl-n_rep - sum(1 for d in delays if results[d]['status']=='low_n')}, "
              f"repaired={n_rep}, excluded_step={n_excl}, low_n={sum(1 for d in delays if results[d]['status']=='low_n')}")

        # --- pass 3: final harmonics per block
        rows = []
        for d in delays:
            dd = data[d]
            r = results[d]
            beta_use = dd["beta"]
            if r["status"] == "repaired":
                beta_use = r["beta"][1]
            if beta_use is None:
                rows.append((d, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                             dd["n"], "low_n", np.nan, np.nan))
                continue
            A0 = beta_use[0]
            An = {}
            phin = {}
            for i, nh in enumerate(HARM):
                a = beta_use[1 + 2 * i]
                b = beta_use[2 + 2 * i]
                c = (a - 1j * b) / 2  # I = A0 + Re[(a - i b) e^{i n phi}] with c_n=(a-ib)/2
                An[nh] = 2 * abs(c)
                phin[nh] = np.angle(c) / nh
            chi2 = float(np.mean((dd["I"] - design(dd["phs"], HARM) @ beta_use) ** 2) / sigma_S**2) \
                if dd["beta"] is not None else np.nan
            rows.append((d, A0, An[1], An[2], phin[1], An[3], An[4], dd["n"], r["status"],
                         chi2, results[d].get("step_z", np.nan)))

        arr = np.array([(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]) for r in rows],
                       dtype=[("t", "f8"), ("A0", "f8"), ("A1", "f8"), ("A2", "f8"),
                              ("phi0", "f8"), ("A3", "f8"), ("A4", "f8"), ("nph", "f8")])
        status = np.array([r[8] for r in rows])
        chi2 = np.array([r[9] for r in rows])
        stepz = np.array([r[10] for r in rows])
        # residuals per block (time-ordered), for hygiene audit
        resids = np.full((len(delays), 60), np.nan)
        for i, d in enumerate(delays):
            dd = data[d]
            e = dd["resid"][np.argsort(dd["ts_t"])]
            resids[i, :len(e)] = e
        np.savez(f"{BASE}/delay60_{name}_harmonics.npz", **{k: arr[k] for k in arr.dtype.names},
                 status=status, chi2=chi2, step_z=stepz, sigma_S=sigma_S,
                 delta_MHz=delta_MHz)
        np.save(f"{BASE}/delay60_{name}_A1_t.npy", np.vstack([arr["t"], arr["A1"], status == "clean"]))
        np.save(f"{BASE}/delay60_{name}_A2_t.npy", np.vstack([arr["t"], arr["A2"], status == "clean"]))
        print(f"saved delay60_{name}_harmonics.npz, A1/A2 npy")

        # --- compact table at representative delays
        print(f"\nt_gap    A0[uW]   A1[uW]  A2[uW]  A3[uW]  A4[uW]   phi0[deg]  status")
        sel = [-6, -4, -2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14]
        tabd = {}
        for r in rows:
            tabd[round(r[0], 4)] = r
        for s in sel:
            cands = [k for k in tabd if abs(k - s) < 0.031]
            if not cands:
                continue
            r = tabd[cands[0]]
            if np.isnan(r[2]):
                print(f"{r[0]:+8.3f}  {'dropped':>7s}")
                continue
            print(f"{r[0]:+8.3f}  {r[1]*1e6:7.3f}  {r[2]*1e6:7.3f}  {r[3]*1e6:6.3f}  "
                  f"{r[5]*1e6:6.3f}  {r[6]*1e6:6.3f}  {np.degrees(r[4]):9.1f}  {r[8]}")
        f.close()


if __name__ == "__main__":
    main()
