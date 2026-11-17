"""
Figure 1 -- single-tone resonance overview (parametric f/2 response vs drive amplitude and frequency)

PORT of spin-waves/Figure!_antigraity/figure1_single_tone/figure1_single_tone.py into the
mistral pipeline framework (provenance manifest, figures + manifests under reports/mistral/),
per the boss directives:

* Fig 1 is SINGLE TONE ONLY: same data sources as the original script
  (Figure!/Figure1_arnold_plot/data.h5, Figure!/arnold_tongue_zoom/data2.h5, read-only),
  same channels (`halfaverage` f/2 observable, `full` drive calibration), same reduction
  (trim, Hampel despike with protected narrow line, normalised-Gaussian smoothing of the maps
  in calibrated drive power, broken-stick/hinge analysis, spectra, drive-chain control).
* Style conventions of the reference kept: rcParams (font 6, etc.), 7.2-in-wide inch-based
  manual layout, panel letters as bold fig.text, maps = pcolormesh(rasterized) on magma with
  cmap.set_bad("white") for missing cells, horizontal colorbar in its own cax ABOVE each panel
  with top ticks, LogNorm(vmin = 3 sigma, vmax = max)/linear 99.5th-pct zoom scale.
* Added panel (h): pulse-timing diagram for the two-tone experiments of Figs 3-4 (tone-2
  delayed onset / gap vs overlap; AWG 64 GSa/s, 1280-sample = 20 ns timing quanta).
* Added panel (i): RESERVED slot for the sample characterisation (FMR/mode-spectrum panel of
  the ~12 BiYIG microdisks); content lands in a later revision (fallback: eigen_27mT_*.npz).
* No dual-tone result content (the earlier detune_fine031126 pcolor maps are withdrawn to
  Figs 2-5).

Only h5py, numpy, scipy, matplotlib are used.  Source HDF5 files are opened READ-ONLY.

Data layout (derived from the files)
------------------------------------
Every HDF5 group is one acquisition batch named  "<a>_<f>_<channel>_<YYYYMMDD>_<HHMMSS>_<index>_0"
    a        AWG drive amplitude  [mV]  (single sine tone, waveform a/1000*sin(2 pi f t))
    f        drive frequency      [MHz]
    channel  full        : SA centred at f   (span 10 MHz, RBW 20 kHz)
             half        : SA centred at f/2 (span 5 MHz) -- parametric response, continuous drive
             halfaverage : SA centred at f/2 (span 5 MHz) -- same, drive pulsed via AWG sequence (50:3), SA averaged
Each group holds   trace  (n_repeat, 2, 1001) float32 :  [:,0,:] SA frequency axis [Hz],  [:,1,:] baseline-subtracted
                   linear power (10**(dBm/10) - baseline, i.e. mW per RBW bin; can be negative),  timestamp (n_repeat, 1)

Run:  python fig1_overview.py            # extract (cached) + figure + validation
      python fig1_overview.py --reextract
"""
import re
import sys
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/pipeline")
from provenance import Provenance

MISTRAL = Path("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral")
FIG_ROOT = MISTRAL / "figures"
CACHE = MISTRAL / "tmp" / "cache"
CACHE.mkdir(parents=True, exist_ok=True)

DATA_ROOT = Path("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/Figure!")  # read only
FILE_BROAD = DATA_ROOT / "Figure1_arnold_plot" / "data.h5"
FILE_ZOOM = DATA_ROOT / "arnold_tongue_zoom" / "data2.h5"

KEY_RE = re.compile(r"^(\d+)_(\d+)_(half|full|halfaverage)_(\d{8})_(\d{6})_(\d+)_(\d+)$")

p = Provenance("fig1_overview")
p.step("revision_1_pcolor",
       params={"boss_feedback": "panel (c) was to become dual-tone pcolor maps from detune_fine031126.h5"},
       note="SUPERSEDED by revision 2: dual-tone content withdrawn to Figs 2-5")
p.step("revision_2_single_tone_port",
       params={"boss_directive": "Fig 1 = port of figure1_single_tone.py (single tone only), "
                                 "provenance framework, style conventions of the reference, "
                                 "plus pulse-timing panel (h) and reserved sample-characterisation slot (i)"},
       note="current revision")
for _f in (FILE_BROAD, FILE_ZOOM):
    p.record_source(str(_f), note="single-tone parameter sweep, read-only source")

# ----------------------------------------------------------------------------------------------
# Extraction (selective reads only; repeat-weighted pooling across acquisition groups)
# ----------------------------------------------------------------------------------------------
TRIM_MIN_REPEATS = 5     # points with at least this many repeats lose their highest and their lowest repeat (by bin-summed power)


def _group_stats(key, tr):
    """Per-group reduction used by extract(): sums and extreme-repeat bookkeeping of one HDF5 group."""
    m = KEY_RE.match(key)
    a, f = float(m.group(1)), float(m.group(2))
    spec = tr[:, 1, :].astype(np.float64)          # (n_repeat, 1001) power slice only
    fax = tr[0, 0, :].astype(np.float64)           # SA axis [Hz]
    pint = spec.sum(1)                             # per-repeat bin-summed power
    imax, imin = int(np.argmax(pint)), int(np.argmin(pint))
    return dict(a=a, f=f, n=int(spec.shape[0]), fax=fax,
                spec_sum=spec.sum(0), pint_sum=float(pint.sum()), pint_sq=float((pint ** 2).sum()),
                pint_max=float(pint[imax]), spec_max=spec[imax].copy(),
                pint_min=float(pint[imin]), spec_min=spec[imin].copy())


def reduce_acc(acc, trim):
    n_trimmed = 0
    ent = {}
    for k, e in acc.items():
        n_, s1, s2 = e["n"], e["pint_sum"], e["pint_sq"]
        if trim and n_ >= TRIM_MIN_REPEATS:
            s1 = s1 - e["pint_max"] - e["pint_min"]
            s2 = s2 - e["pint_max"] ** 2 - e["pint_min"] ** 2
            n_ -= 2
            n_trimmed += 1
        ent[k] = (n_, s1, s2)
    keys_sorted = sorted(acc)
    n_ = np.array([ent[k][0] for k in keys_sorted], float)
    s1 = np.array([ent[k][1] for k in keys_sorted])
    s2 = np.array([ent[k][2] for k in keys_sorted])
    with np.errstate(invalid="ignore", divide="ignore"):
        var = np.where(n_ > 1, (s2 - s1 ** 2 / n_) / (n_ - 1), np.nan)
    return dict(
        pint_sem=np.sqrt(np.maximum(var, 0) / n_),
        a=np.array([k[0] for k in keys_sorted]),
        f=np.array([k[1] for k in keys_sorted]),
        n_rep=np.array([ent[k][0] for k in keys_sorted]),
        n_groups=np.array([acc[k]["n_g"] for k in keys_sorted]),
        spec_mean=np.array([acc[k]["spec_sum"] / acc[k]["n"] for k in keys_sorted], dtype=np.float32),
        sa_freq=np.array([acc[k]["fax"] for k in keys_sorted]),
        n_points_trimmed=np.array(n_trimmed),
    )


def extract_pair(path, channel, tag, reextract=False):
    """Read all groups of one channel selectively (provenance-wrapped), accumulate per (a, f) over *all*
    groups/repeats (several runs revisit the same point); returns the (untrimmed, trimmed) reductions.
    Only trace[:, 1, :] and trace[0, 0, :] are read -- no whole-file loads.  Cached in tmp/cache."""
    cfiles = (CACHE / f"{tag}.npz", CACHE / f"{tag}T.npz")
    if all(c.exists() for c in cfiles) and not reextract:
        out = []
        for cf in cfiles:
            z = np.load(cf, allow_pickle=False)
            out.append({k: z[k] for k in z.files})
        p.step("extract_cache_hit", params={"tag": tag, "cache": [c.name for c in cfiles]},
               result={"n_points": int(out[0]["a"].size)})
        return out
    acc = {}
    skipped = {"unmatched_or_other_channel": 0}

    def _sel(k):
        m = KEY_RE.match(k)
        if m is None or m.group(3) != channel:
            skipped["unmatched_or_other_channel"] += 1
            return False
        return True

    per_group = p.load_traces(str(path), select=_sel, loader=_group_stats, pass_agg="all")
    for d in per_group.values():
        e = acc.get((d["a"], d["f"]))
        if e is None:
            acc[(d["a"], d["f"])] = dict(spec_sum=d["spec_sum"], n=d["n"], n_g=1, fax=d["fax"],
                                         pint_sum=d["pint_sum"], pint_sq=d["pint_sq"],
                                         pint_max=d["pint_max"], spec_max=d["spec_max"],
                                         pint_min=d["pint_min"], spec_min=d["spec_min"])
        else:
            if not np.allclose(e["fax"], d["fax"], rtol=0, atol=1.0):
                raise RuntimeError(f"SA axis mismatch for repeated point {(d['a'], d['f'])}")
            e["spec_sum"] = e["spec_sum"] + d["spec_sum"]
            e["n"] += d["n"]; e["n_g"] += 1
            e["pint_sum"] += d["pint_sum"]; e["pint_sq"] += d["pint_sq"]
            if d["pint_max"] > e["pint_max"]:
                e["pint_max"], e["spec_max"] = d["pint_max"], d["spec_max"]
            if d["pint_min"] < e["pint_min"]:
                e["pint_min"], e["spec_min"] = d["pint_min"], d["spec_min"]
    n_keys = len(per_group)
    del per_group
    raw = reduce_acc(acc, trim=False)
    trm = reduce_acc(acc, trim=True)
    for red, cf in zip((raw, trm), cfiles):
        np.savez_compressed(cf, **red)
    p.step("extract_selective_read", params={"tag": tag, "channel": channel, "path": str(path)},
           note="trace[:,1,:] + trace[0,0,:] slices only; trim = highest/lowest repeat dropped when n >= 5",
           result={"n_hdf5_groups_read": int(n_keys), "skipped_other_keys": int(skipped["unmatched_or_other_channel"]),
                   "n_points": int(raw["a"].size), "n_points_trimmed": int(trm["n_points_trimmed"])})
    return raw, trm


# ----------------------------------------------------------------------------------------------
# Reduction / helper functions  (ported from the reference script)
# ----------------------------------------------------------------------------------------------
NW = 1e-6   # 1 nW = 1e-6 mW.  Raw stored unit is mW per SA bin (10**(dBm/10)), we report nW.


def hampel_despike(spec, protect=None, win=15, nsig=5.0, dil=3):
    from scipy.ndimage import median_filter, binary_dilation
    spec = np.asarray(spec, dtype=np.float64)
    med = median_filter(spec, size=(1, win), mode="nearest")
    r = spec - med
    sig = 1.4826 * np.maximum(median_filter(np.abs(r), size=(1, 101), mode="nearest"), 1e-30)
    mask = np.abs(r) > nsig * sig
    if protect is not None:
        mask &= ~protect
    if dil:
        mask = binary_dilation(mask, structure=np.ones((1, 2 * dil + 1), bool))
        if protect is not None:
            mask &= ~protect
    return np.where(mask, med, spec), mask


PROTECT_A_KHZ = 100.0     # the narrow (RBW 20 kHz) f/2 line of the overview data is never despiked within +-100 kHz


def despike_dataset(d, protect_khz=None):
    out = dict(d)
    protect = None
    if protect_khz is not None:
        off = (d["sa_freq"] - d["f"][:, None] / 2 * 1e6) / 1e3
        protect = np.abs(off) < protect_khz
    clean, mask = hampel_despike(d["spec_mean"], protect)
    out["spec_mean"] = clean.astype(np.float32)
    return out, dict(bins_replaced_fraction=float(mask.mean()),
                     points_with_replaced_bins=float(mask.any(1).mean()))


def pint_nW(d):
    return d["spec_mean"].astype(np.float64).sum(1) / NW


def noise_sigma(P):
    neg = P[P < 0]
    return float(np.median(np.abs(neg)) / 0.6745)


def cut_at_frequency(d, f0, quantity, nneigh=2, max_dist=3.0):
    a_all = np.unique(d["a"])
    vals, flags = [], []
    for a in a_all:
        sel = np.flatnonzero(d["a"] == a)
        fs = d["f"][sel]
        hit = sel[fs == f0]
        if hit.size:
            vals.append(quantity[hit[0]]); flags.append(0); continue
        left = sel[(fs < f0) & (fs >= f0 - max_dist)]
        right = sel[(fs > f0) & (fs <= f0 + max_dist)]
        if left.size == 0 or right.size == 0:
            vals.append(np.full(np.shape(quantity[0]), np.nan)); flags.append(2); continue
        left = left[np.argsort(d["f"][left])][-nneigh:]
        right = right[np.argsort(d["f"][right])][:nneigh]
        use = np.concatenate([left, right])
        x = d["f"][use]
        vals.append(PchipInterpolator(x, quantity[use], axis=0, extrapolate=False)(f0)); flags.append(1)
    return a_all, np.array(vals), np.array(flags)


def loo_cut(d, f0, quantity, nneigh=2, max_dist=3.0):
    keep = d["f"] != f0
    d2 = {k: (v[keep] if isinstance(v, np.ndarray) and v.shape[:1] == d["f"].shape else v) for k, v in d.items()}
    a_i, v_i, fl = cut_at_frequency(d2, f0, quantity[keep], nneigh, max_dist)
    a_m = d["a"][d["f"] == f0]; v_m = quantity[d["f"] == f0]
    order = np.argsort(a_m)
    a_m, v_m = a_m[order], v_m[order]
    common = np.isin(a_i, a_m)
    a_i, v_i, fl = a_i[common], v_i[common], fl[common]
    ok = fl == 1
    return a_i[ok], v_i[ok], v_m[np.isin(a_m, a_i[ok])]


def smooth_map(a, f, P, f_grid, a_grid, sig_f, sig_a, rel_cut=2.5):
    from scipy.ndimage import gaussian_filter
    df, da = f_grid[1] - f_grid[0], a_grid[1] - a_grid[0]
    ja = np.round((a - a_grid[0]) / da).astype(int); jf = np.round((f - f_grid[0]) / df).astype(int)
    ok = (ja >= 0) & (ja < a_grid.size) & (jf >= 0) & (jf < f_grid.size)
    gw = np.zeros((a_grid.size, f_grid.size)); gp = np.zeros_like(gw)
    np.add.at(gw, (ja[ok], jf[ok]), 1.0); np.add.at(gp, (ja[ok], jf[ok]), P[ok])
    sg = (sig_a / da, sig_f / df)
    num = gaussian_filter(gp, sg, mode="constant", truncate=4.0)
    den = gaussian_filter(gw, sg, mode="constant", truncate=4.0)
    delta = np.zeros((9, 9)); delta[4, 4] = 1.0
    peak = gaussian_filter(delta, sg, mode="constant", truncate=4.0).max()
    valid = den >= peak * np.exp(-0.5 * rel_cut ** 2)
    Z = np.where(valid, num / np.where(den > 0, den, 1), np.nan)
    return Z, float(valid.mean())


def smooth_cut(a, P, a_eval, sig_a):
    K = np.exp(-0.5 * ((a_eval[:, None] - a[None, :]) / sig_a) ** 2)
    return (K * P[None, :]).sum(1) / K.sum(1)


def threshold_amplitude(a, P, P_cut):
    o = np.argsort(a); a, P = a[o], P[o]
    below = np.flatnonzero(P < P_cut)
    if below.size == 0 or below[-1] == len(a) - 1:
        return np.nan
    i = below[-1]
    return a[i] + (P_cut - P[i]) * (a[i + 1] - a[i]) / (P[i + 1] - P[i])


def raw_point_pint(path, channel, a, f):
    num_, n_tot, axis = 0.0, 0, None
    with __import__("h5py").File(path, "r") as h:
        for k in h.keys():
            m = KEY_RE.match(k)
            if m and m.group(3) == channel and float(m.group(1)) == a and float(m.group(2)) == f:
                tr = h[k]["trace"][...]
                num_ += tr[:, 1, :].astype(np.float64).sum()
                n_tot += tr.shape[0]
                axis = tr[0, 0, :]
    return num_ / n_tot / NW, n_tot, axis


def hinge_fit(a, P, lo, hi, step=0.5, margin=6.0):
    m = (a >= lo) & (a <= hi)
    a, P = a[m], P[m]
    best = None
    for ak in np.arange(lo + margin, hi - margin + 1e-9, step):
        X = np.column_stack([np.ones_like(a), a - ak, np.maximum(a - ak, 0.0)])
        c = np.linalg.lstsq(X, P, rcond=None)[0]
        sse = float(((X @ c - P) ** 2).sum())
        if best is None or sse < best[0]:
            best = (sse, ak, c)
    X1 = np.column_stack([np.ones_like(a), a])
    c1 = np.linalg.lstsq(X1, P, rcond=None)[0]
    sse1 = float(((X1 @ c1 - P) ** 2).sum())
    sse, ak, c = best
    Fstat = (sse1 - sse) / 2.0 / (sse / (a.size - 4))
    return dict(ak=float(ak), s1=float(c[1]), s2=float(c[1] + c[2]), c0=float(c[0]), F=float(Fstat), n=int(a.size),
                lo=float(lo), hi=float(hi))


def hinge_model(h, aa):
    return h["c0"] + h["s1"] * (aa - h["ak"]) + (h["s2"] - h["s1"]) * np.maximum(aa - h["ak"], 0.0)


def hinge_bootstrap(a, P, lo, hi, nb=500, seed=1):
    rng = np.random.default_rng(seed)
    m = (a >= lo) & (a <= hi)
    a, P = a[m], P[m]
    r = []
    for _ in range(nb):
        i = rng.integers(0, a.size, a.size)
        h = hinge_fit(a[i], P[i], lo, hi)
        r.append((h["ak"], h["s1"], h["s2"]))
    return np.percentile(np.array(r), [16, 50, 84], axis=0)


def local_slope(a, P, half_width=5.5):
    sl, se = np.full(a.size, np.nan), np.full(a.size, np.nan)
    for i, ai in enumerate(a):
        m = np.abs(a - ai) <= half_width
        if m.sum() >= 6:
            c, cov = np.polyfit(a[m], P[m], 1, cov=True)
            sl[i], se[i] = c[0], np.sqrt(cov[0, 0])
    return sl, se


def fwhm_contiguous(off, s, halfwin=700.0):
    sel = np.flatnonzero(np.abs(off) < halfwin)
    i0 = sel[np.argmax(s[sel])]
    pk = s[i0]
    if pk <= 0:
        return np.nan, np.nan
    half = pk / 2
    j = i0
    while j > 0 and s[j] > half:
        j -= 1
    k = i0
    while k < s.size - 1 and s[k] > half:
        k += 1
    xl = off[j] + (half - s[j]) / (s[j + 1] - s[j]) * (off[j + 1] - off[j])
    xr = off[k - 1] + (half - s[k - 1]) / (s[k] - s[k - 1]) * (off[k] - off[k - 1])
    return float(xr - xl), float(0.5 * (xr + xl))


def synthetic_hinge_test():
    rng = np.random.default_rng(3)
    a = np.arange(285.0, 336.0)
    truth = dict(ak=315.0, s1=0.8, s2=3.0, c0=30.0)
    P = hinge_model(truth, a) + rng.normal(0, 1.5, a.size)
    h = hinge_fit(a, P, 285, 335)
    assert abs(h["ak"] - truth["ak"]) < 2.0 and abs(h["s2"] - truth["s2"]) < 0.4, h
    return h


RBW_FULL_HZ = 20e3
ENBW_FACTOR = 1.056


def fit_drive_model(d_full):
    ok = d_full["a"] > 0
    a, f = d_full["a"][ok], d_full["f"][ok]
    P = d_full["spec_mean"][ok].astype(np.float64).sum(1) / NW
    fs = np.unique(f); idx = np.searchsorted(fs, f)
    X = np.column_stack([np.log(a)] + [(idx == i).astype(float) for i in range(fs.size)])
    coef = np.linalg.lstsq(X, np.log(P), rcond=None)[0]
    res = np.log(P) - X @ coef
    return float(coef[0]), dict(zip(fs, np.exp(coef[1:]))), float(res.std()), (a, f, P)


def drive_power_nW(d_half, d_full, model):
    e, cf, _, (a_f, f_f, P_f) = model
    look = {(float(x), float(y)): z for x, y, z in zip(a_f, f_f, P_f)}
    dfb = (d_full["sa_freq"][0, -1] - d_full["sa_freq"][0, 0]) / (d_full["sa_freq"].shape[1] - 1)
    out, used_fit = np.empty(d_half["a"].size), 0
    for i, (a, f) in enumerate(zip(d_half["a"], d_half["f"])):
        v = look.get((float(a), float(f)))
        if v is None:
            v = cf[f] * a ** e if f in cf else np.nan
            used_fit += 1
        out[i] = v * dfb / (ENBW_FACTOR * RBW_FULL_HZ)
    return out, used_fit


def to_dBm(P_nW):
    return 10.0 * np.log10(np.maximum(P_nW, 1e-30) * 1e-6)


def group_scatter_factor(path, f0, a_min=250.0):
    groups = {}
    with __import__("h5py").File(path, "r") as h:
        for k in h.keys():
            m = KEY_RE.match(k)
            if not m or m.group(3) != "halfaverage" or float(m.group(2)) != f0 or float(m.group(1)) < a_min:
                continue
            pr = h[k]["trace"][:, 1, :].astype(np.float64).sum(1) / NW
            if pr.size >= 2:
                groups.setdefault(float(m.group(1)), []).append((pr.size, pr.mean(), pr.var(ddof=1)))
    z = []
    for L in groups.values():
        for i in range(len(L)):
            for j in range(i + 1, len(L)):
                se = np.sqrt(L[i][2] / L[i][0] + L[j][2] / L[j][0])
                if se > 0:
                    z.append((L[i][1] - L[j][1]) / se)
    z = np.array(z)
    return float(z.std()), int(z.size)


def cv_compare(a, P, lo, hi):
    m = (a >= lo) & (a <= hi)
    a, P = a[m], P[m]
    n = a.size
    x = (a - a.mean()) / 25.0
    res = {}
    for deg in (1, 2, 3, 4):
        e = [P[i] - np.polyval(np.polyfit(x[np.arange(n) != i], P[np.arange(n) != i], deg), x[i]) for i in range(n)]
        res[f"poly{deg}"] = float(np.sqrt(np.mean(np.square(e))))
    e = []
    for i in range(n):
        k = np.arange(n) != i
        e.append(P[i] - hinge_model(hinge_fit(a[k], P[k], lo, hi), a[i:i + 1])[0])
    res["hinge"] = float(np.sqrt(np.mean(np.square(e))))
    return res


def wild_bootstrap_p(a, P, lo, hi, deg=3, nb=1000, seed=7):
    rng = np.random.default_rng(seed)
    m = (a >= lo) & (a <= hi)
    aw, Pw = a[m], P[m]
    xs = (aw - aw.mean()) / 25.0
    sm = np.polyval(np.polyfit(xs, Pw, deg), xs)
    res = Pw - sm
    Fobs = hinge_fit(aw, Pw, lo, hi)["F"]
    Fnull = np.array([hinge_fit(aw, sm + res * rng.choice([-1.0, 1.0], res.size), lo, hi)["F"] for _ in range(nb)])
    return float(np.mean(Fnull >= Fobs)), float(np.percentile(Fnull, 95)), float(Fobs)


def raw_reduced_point(path, channel, a, f, trim, protect_khz):
    reps = []
    with __import__("h5py").File(path, "r") as h:
        for k in h.keys():
            m = KEY_RE.match(k)
            if m and m.group(3) == channel and float(m.group(1)) == a and float(m.group(2)) == f:
                tr = h[k]["trace"][...]
                reps.append(tr[:, 1, :].astype(np.float64))
                fax = tr[0, 0, :].astype(np.float64)
    R = np.concatenate(reps, axis=0)
    if trim and R.shape[0] >= TRIM_MIN_REPEATS:
        ps = R.sum(1)
        keep = np.ones(R.shape[0], bool); keep[np.argmax(ps)] = False; keep[np.argmin(ps)] = False
        R = R[keep]
    mean = R.mean(0)[None, :]
    protect = None
    if protect_khz is not None:
        protect = np.abs((fax - f / 2 * 1e6) / 1e3)[None, :] < protect_khz
    clean, _ = hampel_despike(mean, protect)
    return float(clean.sum() / NW), clean[0]


# ----------------------------------------------------------------------------------------------
# Figure constants  (ported)
# ----------------------------------------------------------------------------------------------
Y_STEP = 0.03
Y_STEP_B = 0.01
SIG_A_CANDIDATES = [(1.5, 0.10), (2.0, 0.13), (2.5, 0.17), (3.0, 0.20), (4.0, 0.27)]
SIG_ZOOM = (1.0, 0.03)
SIG_CUT_MV = 2.0
SIG_SPEC_BINS = 5
SIG_SPEC20_BINS = 2
A_SLICE_MV, SLICE_BAND_MV, SIG_SLICE_MHZ = 335.0, 2.0, 2.5
K1_WINDOW, K2_WINDOW = (285.0, 335.0), (322.0, 350.0)
SPEC_TARGETS = [300, 310, 325, 335, 348]
C_BELOW, C_ABOVE = ["#9ecae1", "#3182bd"], ["#fdae6b", "#e6550d", "#a63603"]
C_K1, C_K2 = "#1f78b4", "#238b45"
A_MIN_ZOOM = 250.0
F_PLOT_MAX_B = 4080.0
F_EXC = 4065.0


def support_neff(y, f, f_grid, y_grid, sig_f, sig_y):
    from scipy.ndimage import gaussian_filter
    df, dy = f_grid[1] - f_grid[0], y_grid[1] - y_grid[0]
    ja = np.round((y - y_grid[0]) / dy).astype(int); jf = np.round((f - f_grid[0]) / df).astype(int)
    ok = (ja >= 0) & (ja < y_grid.size) & (jf >= 0) & (jf < f_grid.size)
    gw = np.zeros((y_grid.size, f_grid.size)); np.add.at(gw, (ja[ok], jf[ok]), 1.0)
    sg = (sig_y / dy, sig_f / df)
    s1 = gaussian_filter(gw, sg, mode="constant", truncate=4.0)
    s2 = gaussian_filter(gw, (sg[0] / np.sqrt(2), sg[1] / np.sqrt(2)), mode="constant", truncate=4.0) / (4 * np.pi * sg[0] * sg[1])
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(s2 > 0, s1 ** 2 / s2, np.nan)


def threshold_curve(Z, y_grid, thr):
    out = np.full(Z.shape[1], np.nan)
    for j in range(Z.shape[1]):
        col = Z[:, j]
        bad = np.flatnonzero(~(col > thr))
        out[j] = y_grid[0] if bad.size == 0 else (y_grid[bad[-1] + 1] if bad[-1] + 1 < y_grid.size else np.nan)
    return out


def matplotlib_to_rgb(c):
    import matplotlib.colors as mcolors
    return mcolors.to_rgb(c)


def draw_sequence(ax, Wb, Hb):
    """Single-tone excitation + measurement sequence schematic (facts from run.py / run2.py):
    AWG 64 GS/s, segment 1 = a/1000 sin(2 pi f t) over 64000 samples (1 us), segment 2 = 64000 zeros (1 us);
    `halfaverage`: sequence ((1,2),(50,3)) = 50 x 1 us on, 3 x 1 us off, repeated; `half`: segment 1 continuously.
    Per repeat the SA is swept twice (0.25 s each, centred at f/2, 5 MHz span, 1001 points)."""
    from matplotlib.collections import LineCollection
    from matplotlib.patches import FancyBboxPatch, Circle, Rectangle
    BLUE, LBLUE, GREY, LGREY, RED = "#2b7bba", "#d4e6f4", "#6e6e6e", "#ececec", "#c0392b"
    ax.set_xlim(0, Wb); ax.set_ylim(0, Hb); ax.axis("off")
    fs_h, fs_t = 5.2, 4.7

    def heading(y, n, text):
        ax.add_patch(Circle((0.075, y + 0.017), 0.048, fc="#404040", ec="none", zorder=3))
        ax.text(0.075, y + 0.017, str(n), color="white", fontsize=4.6, fontweight="bold", ha="center", va="center", zorder=4)
        ax.text(0.16, y + 0.017, text, fontsize=fs_h, fontweight="bold", color="#202020", va="center", ha="left")

    def dim(x0, x1, y, text, ty, color="#404040", fs=fs_t):
        ax.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="<|-|>", lw=0.4, color=color, mutation_scale=3.2,
                                                                    shrinkA=0, shrinkB=0))
        ax.text(0.5 * (x0 + x1), ty, text, fontsize=fs, color=color, ha="center", va="center")

    heading(1.875, 1, "AWG waveform (segment of 1 µs)")
    x = np.linspace(0.11, 1.07, 400); y0, amp = 1.615, 0.090
    yy = y0 + amp * np.sin(2 * np.pi * (x - 0.11) / 0.16)
    segs = np.stack([np.column_stack([x[:-1], yy[:-1]]), np.column_stack([x[1:], yy[1:]])], axis=1)
    alpha = np.clip((1.07 - x[:-1]) / 0.30, 0.0, 1.0) ** 0.8
    lc = LineCollection(segs, colors=[(*matplotlib_to_rgb(BLUE), a_) for a_ in alpha], linewidths=0.9, zorder=3)
    ax.add_collection(lc)
    ax.plot([0.11, 1.12], [y0, y0], color="#bdbdbd", lw=0.35, zorder=1)
    ax.text(1.12, y0, "· · ·", fontsize=6, color=BLUE, va="center", ha="left")
    ax.annotate("", xy=(0.075, y0 + amp), xytext=(0.075, y0), arrowprops=dict(arrowstyle="<|-|>", lw=0.4, color="#404040", mutation_scale=3, shrinkA=0, shrinkB=0))
    ax.text(0.03, y0 + 0.5 * amp, "$a$", fontsize=5.2, color="#202020", ha="center", va="center", style="italic")
    dim(0.11 + 0.16, 0.11 + 0.32, y0 + amp + 0.035, "", 0)
    ax.text(0.11 + 0.24, y0 + amp + 0.085, "$T=1/f$", fontsize=fs_t, color="#404040", ha="center", va="center")
    ax.text(1.45, y0 + 0.035, "$a\\,\\sin(2\\pi f t)$", fontsize=5.4, color=BLUE, ha="left", va="center")
    ax.text(1.45, y0 - 0.060, "$f$ = 4065 MHz: 4065 periods\nper segment, 64 GS/s", fontsize=4.4, color=GREY, ha="left", va="center", linespacing=1.15)

    heading(1.445, 2, "Drive: pulsed vs continuous")
    per, on, off = 0.735, 0.690, 0.045
    xs0 = 0.11
    yb1, hb = 1.205, 0.125
    yb2 = 0.930

    def burst(x0, x1, y, h, face=LBLUE):
        ax.add_patch(Rectangle((x0, y), x1 - x0, h, fc=face, ec=BLUE, lw=0.5, zorder=2))
        for xx in np.arange(x0 + 0.01, x1 - 0.005, 0.0125):
            ax.plot([xx, xx], [y + 0.012, y + h - 0.012], color=BLUE, lw=0.12, alpha=0.55, zorder=2.5)

    xcur = xs0
    for k in range(3):
        burst(xcur, xcur + on, yb1, hb)
        xcur += on
        if k < 2:
            ax.plot([xcur, xcur + off], [yb1 + hb / 2] * 2, color="#9a9a9a", lw=0.5, zorder=2)
            xcur += off
    ax.annotate("", xy=(xcur + 0.17, yb1 + hb / 2), xytext=(xcur + 0.015, yb1 + hb / 2),
                arrowprops=dict(arrowstyle="-|>", lw=0.5, color=BLUE, mutation_scale=3.5, shrinkA=0, shrinkB=0))
    ax.text(xs0 + 0.5 * on, yb1 + hb / 2, "50 µs on", fontsize=fs_t, color="#12426b", ha="center", va="center", zorder=5, fontweight="bold")
    xg = xs0 + on + off / 2
    ax.annotate("", xy=(xg, yb1 - 0.004), xytext=(xg + 0.16, yb1 - 0.085),
                arrowprops=dict(arrowstyle="-|>", lw=0.4, color=RED, mutation_scale=3, shrinkA=0, shrinkB=0))
    ax.text(xg + 0.18, yb1 - 0.095, "3 µs off", fontsize=fs_t, color=RED, ha="left", va="center")
    dim(xs0, xs0 + per, yb1 - 0.038, "", 0)
    ax.text(xs0 + 0.5 * per, yb1 - 0.085, "53 µs period", fontsize=fs_t, color="#404040", ha="center", va="center")
    ax.text(xs0, yb1 + hb + 0.040, "halfaverage", fontsize=fs_t, color=BLUE, ha="left", va="center", fontweight="bold")
    ax.text(xs0 + 1.08, yb1 + hb + 0.040, "50 µs on / 3 µs off, 94 % duty", fontsize=4.4, color=GREY, ha="left", va="center")
    burst(xs0, xs0 + 3 * on + 2 * off, yb2, hb)
    ax.annotate("", xy=(xs0 + 3 * on + 2 * off + 0.17, yb2 + hb / 2), xytext=(xs0 + 3 * on + 2 * off + 0.015, yb2 + hb / 2),
                arrowprops=dict(arrowstyle="-|>", lw=0.5, color=BLUE, mutation_scale=3.5, shrinkA=0, shrinkB=0))
    ax.text(xs0 + 0.5 * (3 * on + 2 * off), yb2 + hb / 2, "half: continuous drive", fontsize=fs_t, color="#12426b", ha="center", va="center", zorder=5, fontweight="bold")

    heading(0.775, 3, "Measurement of one $(a, f)$ point (SA at $f/2$)")
    yb, hb3 = 0.40, 0.27

    def box(x0, x1, fc, ec, l1, l2):
        ax.add_patch(FancyBboxPatch((x0, yb), x1 - x0, hb3, boxstyle="round,pad=0,rounding_size=0.03", fc=fc, ec=ec, lw=0.5, zorder=2))
        ax.text(0.5 * (x0 + x1), yb + hb3 * 0.67, l1, fontsize=fs_t, fontweight="bold", color="#202020", ha="center", va="center", zorder=4)
        ax.text(0.5 * (x0 + x1), yb + hb3 * 0.27, l2, fontsize=4.3, color="#404040", ha="center", va="center", zorder=4)

    bx = [(0.06, 0.78), (0.90, 1.62), (1.74, Wb - 0.04)]
    box(*bx[0], LGREY, "#8a8a8a", "baseline sweep", "AWG idle · 0.25 s")
    box(*bx[1], LBLUE, BLUE, "signal sweep", "AWG on · 0.25 s")
    box(*bx[2], "#fdeee3", RED, "$P_{\\rm on}-P_{\\rm base}$", "$\\Sigma P$ spectrum")
    for (xa, xb) in ((bx[0][1], bx[1][0]), (bx[1][1], bx[2][0])):
        ax.annotate("", xy=(xb, yb + hb3 / 2), xytext=(xa, yb + hb3 / 2),
                    arrowprops=dict(arrowstyle="-|>", lw=0.5, color="#303030", mutation_scale=3.5, shrinkA=0, shrinkB=0))
    xl0, xl1, yl = 0.5 * (bx[0][0] + bx[0][1]), 0.5 * (bx[2][0] + bx[2][1]), 0.265
    ax.plot([xl1, xl1], [yb, yl], color="#303030", lw=0.5, zorder=1)
    ax.plot([xl1, xl0], [yl, yl], color="#303030", lw=0.5, zorder=1)
    ax.annotate("", xy=(xl0, yb - 0.002), xytext=(xl0, yl), arrowprops=dict(arrowstyle="-|>", lw=0.5, color="#303030", mutation_scale=3.5, shrinkA=0, shrinkB=0))
    ax.text(0.5 * (xl0 + xl1), yl + 0.045, "$\\times$ repeats (adaptive, up to 15)", fontsize=fs_t, color="#202020", ha="center", va="center",
            bbox=dict(fc="white", ec="none", pad=0.4))
    ax.text(0.5 * Wb, 0.155, "then channel full: SA at $f$ (10 MHz span) → drive power $P_{\\rm dr}$", fontsize=4.5, color=GREY, ha="center", va="center")
    ax.text(0.5 * Wb, 0.075, "$(a,f)$: scrambled Sobol points, 1 mV × 1 MHz grid; SA 5 MHz span, 1001 pts", fontsize=4.3, color=GREY, ha="center", va="center")


def draw_timing(ax, Wb, Hb):
    """Two-tone pulse-timing conventions of Figs 3-4 (adapted from the earlier fig1 draft; AWG facts from the
    phase-gating scripts).  For reading Figs 3-4 only -- no dual-tone data content in Fig 1."""
    from matplotlib.patches import Rectangle, Circle
    BLUE, ORANGE, GREY, RED = "#2b7bba", "#e08a2e", "#6e6e6e", "#c0392b"
    ax.set_xlim(0, Wb); ax.set_ylim(0, Hb); ax.axis("off")
    fs_h, fs_t, fs_s = 5.2, 4.7, 4.3

    def heading(y, n, text):
        ax.add_patch(Circle((0.075, y + 0.017), 0.048, fc="#404040", ec="none", zorder=3))
        ax.text(0.075, y + 0.017, str(n), color="white", fontsize=4.6, fontweight="bold", ha="center", va="center", zorder=4)
        ax.text(0.16, y + 0.017, text, fontsize=fs_h, fontweight="bold", color="#202020", va="center", ha="left")

    heading(1.32, 1, "Two-tone timing of Figs 3-4 (AWG 64 GSa/s; quantum = 1280 samples = 20 ns)")

    # row 1: sequential (non-overlap) pulses, 2-us record
    ax.text(0.11, 1.17, "sequential pulses ('varyphase'): tone-2 delayed onset, no overlap", fontsize=fs_t, color="#202020", va="center")
    y0, hh = 1.02, 0.085
    x0r = 0.62; scale = 1.80 / 2.0                       # inches per µs over the 2-µs record
    ax.add_patch(Rectangle((x0r + 0.70 * scale, y0), 0.30 * scale, hh, fc="#d4e6f4", ec=BLUE, lw=0.5, zorder=2))
    ax.add_patch(Rectangle((x0r + 1.00 * scale, y0), 0.40 * scale, hh, fc="#fbe3c9", ec=ORANGE, lw=0.5, zorder=2))
    ax.plot([x0r, x0r + 2.0 * scale], [y0, y0], color="#bdbdbd", lw=0.35, zorder=1)
    ax.plot([x0r + 1.00 * scale] * 2, [y0 - 0.032, y0 + hh + 0.032], color=RED, lw=0.6, ls=(0, (2, 1.5)), zorder=3)
    ax.text(x0r + 0.85 * scale, y0 + hh / 2, "tone 1\n300 ns", fontsize=fs_s, color="#12426b", ha="center", va="center", zorder=4, linespacing=1.1)
    ax.text(x0r + 1.20 * scale, y0 + hh / 2, "tone 2\n400 ns", fontsize=fs_s, color="#7a4a10", ha="center", va="center", zorder=4, linespacing=1.1)
    ax.text(x0r + 2.04 * scale, y0 + hh / 2, "record 128000 samples = 2 µs;\ntone phases $[0, \\varphi]$ at $t_{0,\\rm ref}$",
            fontsize=fs_s, color=GREY, ha="left", va="center", linespacing=1.2)
    ax.annotate("", xy=(x0r + 1.20 * scale, y0 - 0.010), xytext=(x0r + 1.20 * scale, y0 - 0.062),
                arrowprops=dict(arrowstyle="-|>", lw=0.4, color=ORANGE, mutation_scale=3, shrinkA=0, shrinkB=0))
    ax.text(x0r + 1.24 * scale, y0 - 0.062, "tone-2 onset carries the scanned 'delay' of Figs 3-4", fontsize=fs_s, color=ORANGE, ha="left", va="center")

    # row 2: overlap sequence, 51-us record (compressed)
    ax.text(0.11, 0.885, "overlap ('f1_4069_5' 3-segment AWG): tone 2 joins a continuous tone 1", fontsize=fs_t, color="#202020", va="center")
    y1_, h2_ = 0.70, 0.085
    x0s = 0.62; scon = 1.80 / 51.0
    ax.add_patch(Rectangle((x0s, y1_), 48.0 * scon, h2_, fc="#d4e6f4", ec=BLUE, lw=0.5, zorder=2))
    ax.add_patch(Rectangle((x0s + 2.999 * scon, y1_), (48.0 - 2.999) * scon, h2_ * 0.55, fc="#fbe3c9", ec=ORANGE, lw=0.4, zorder=3))
    ax.add_patch(Rectangle((x0s + 48.0 * scon, y1_), 3.0 * scon, h2_, fc="0.93", ec="none", zorder=1))
    ax.plot([x0s + 2.999 * scon] * 2, [y1_ - 0.032, y1_ + h2_ + 0.032], color=RED, lw=0.6, ls=(0, (2, 1.5)), zorder=4)
    ax.text(x0s + 25.5 * scon, y1_ + h2_ * 0.275, "both tones on, 45 µs (seg 2 × 15)", fontsize=fs_s, color="#7a4a10", ha="center", va="center")
    ax.text(x0s - 0.012, y1_ + h2_ / 2, "tone 1 only", fontsize=fs_s, color="#12426b", ha="right", va="center")
    ax.annotate("contact point $t_{0,\\rm ref}$ = 2.999 µs (phases $[\\varphi/2,-\\varphi/2]$ there)",
                xy=(x0s + 2.999 * scon, y1_ + h2_ + 0.006), xytext=(x0s + 13.0 * scon, y1_ + h2_ + 0.048),
                arrowprops=dict(arrowstyle="-|>", lw=0.4, color=RED, mutation_scale=3, shrinkA=0, shrinkB=0),
                fontsize=fs_s, color=RED, ha="left", va="center")
    for xb, txt in ((3.0, "seg 1"), (25.5, "seg 2 ×15"), (49.5, "seg 3")):
        xx = x0s + xb * scon
        ax.text(xx, y1_ - 0.038, txt, fontsize=fs_s, color="0.35", ha="center", va="top")
    ax.text(0.11, 0.47, "each segment 192000 samples = 3 µs; total 17 × 192000 samples = 51 µs; 64-sample (1 ns) Tukey ramps in/out",
            fontsize=fs_s, color=GREY, ha="left", va="center")
    ax.text(0.11, 0.37, "amplitudes: varyphase A1 = 327 / A2 = 588 mV (delay_vs_f1_fit.py); overlap A1 = 324 mV, A2 = 200-275 mV (f1_4069_5)",
            fontsize=fs_s, color=GREY, ha="left", va="center")
    ax.text(0.11, 0.27, "tone frequencies f1 ≈ 4060-4085 MHz, f2 = f1 + Δ; pulse envelopes on the 1280-sample (20 ns) timing grid; phases at $t_{0,\\rm ref}$",
            fontsize=fs_s, color=GREY, ha="left", va="center")


def draw_placeholder(ax, Wb, Hb):
    """RESERVED slot (boss addendum): FMR / mode-spectrum sample characterisation of the ~12 BiYIG microdisks.
    Content is filled in a later revision once the measured data land; fallback = computed eigenmode spectra
    at B_DC = 27 mT (spin-waves/eigen_27mT_35modes.npz / eigen_27mT_80modes.npz)."""
    from matplotlib.patches import FancyBboxPatch
    ax.set_xlim(0, Wb); ax.set_ylim(0, Hb); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.03, 0.05), 0.94, 0.90, boxstyle="round,pad=0,rounding_size=0.03",
                                fc="#f7f7f7", ec="0.55", lw=0.6, ls=(0, (4, 3))))
    ax.text(0.5, 0.66, "[ reserved: sample characterisation ]", fontsize=5.6, fontweight="bold", color="0.35", ha="center", va="center")
    ax.text(0.5, 0.48, "FMR / mode-spectrum panel of the ~12 BiYIG microdisks\n(where the driven modes of panels a-g come from)",
            fontsize=4.6, color="0.4", ha="center", va="center", linespacing=1.3)
    ax.text(0.5, 0.27, "panel slot per boss addendum; measured data pending,\nfallback: eigenmode spectra at $B_{\\rm DC}$ = 27 mT",
            fontsize=4.4, color="0.5", ha="center", va="center", linespacing=1.3)


def build(reextract=False):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle, ConnectionPatch
    from matplotlib.ticker import LogLocator, NullFormatter, MultipleLocator
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import argrelmin

    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 6, "axes.labelsize": 6.5, "xtick.labelsize": 5.8, "ytick.labelsize": 5.8, "legend.fontsize": 5.4,
        "axes.linewidth": 0.5, "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.minor.width": 0.35, "ytick.minor.width": 0.35,
        "xtick.major.size": 2.0, "ytick.major.size": 2.0, "xtick.minor.size": 1.1, "ytick.minor.size": 1.1,
        "xtick.major.pad": 1.6, "ytick.major.pad": 1.6, "axes.labelpad": 1.6,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "mathtext.default": "regular",
    })

    # =========================================================================================
    # data: observable `halfaverage` (f/2); `full` as the drive-power calibration / control
    # =========================================================================================
    A0, At = extract_pair(FILE_BROAD, "halfaverage", "broad", reextract=reextract)
    B0, Bt = extract_pair(FILE_ZOOM, "halfaverage", "zoom", reextract=reextract)
    assert np.array_equal(A0["a"], At["a"]) and np.array_equal(B0["f"], Bt["f"])
    A, stA = despike_dataset(At, PROTECT_A_KHZ)
    B, stB = despike_dataset(Bt, None)
    Af, _ = extract_pair(FILE_BROAD, "full", "broad_full", reextract=reextract)
    Bf, _ = extract_pair(FILE_ZOOM, "full", "zoom_full", reextract=reextract)
    p.step("analysis_despike", params={"protect_kHz_overview": PROTECT_A_KHZ, "win": 15, "nsig": 5.0, "dil": 3},
           result={"broad": stA, "zoom": stB},
           note="Hampel despike of the trimmed mean spectra; narrow RBW-20 kHz line protected in the overview")
    PA, PB = pint_nW(A), pint_nW(B)
    sigA, sigB = noise_sigma(PA), noise_sigma(PB)
    thrA, thrB = 3 * sigA, 3 * sigB
    num = {"channel": "halfaverage", "sigma_noise_A_nW": sigA, "sigma_noise_B_nW": sigB,
           "sigma_a0_reference_B_nW": float(np.std(PB[B["a"] == 0])),
           "n_points_A": int(A["a"].size), "n_points_B": int(B["a"].size),
           "n_repeats_A": int(A["n_rep"].sum()), "n_repeats_B": int(B["n_rep"].sum()),
           "native_coverage_A": float(A["a"].size / (301 * 491)),
           "mean_sumP_A_nW_around_fexc": float(np.mean(PA[(np.abs(A['f'] - F_EXC) <= 3) & (A['a'] >= 300)])),
           "mean_sumP_B_nW_same_window": float(np.mean(PB[(np.abs(B['f'] - F_EXC) <= 3) & (B['a'] >= 300)]))}
    p.step("analysis_noise", params={"estimator": "median(|P| | P<0)/0.6745 of bin-summed power"},
           result={"sigma_A_nW": sigA, "sigma_B_nW": sigB, "thresholds_3sigma_nW": [thrA, thrB]})

    # trimming / despiking effect on the bin-summed power
    P0B, PtB, P0A, PtA = pint_nW(B0), pint_nW(Bt), pint_nW(A0), pint_nW(At)
    sg0B, sg0A = noise_sigma(P0B), noise_sigma(P0A)
    sB_, sA_ = P0B > 10 * sg0B, P0A > 10 * sg0A
    num.update(
        trim_min_repeats=TRIM_MIN_REPEATS,
        trim_points_fraction_A=float(At["n_points_trimmed"] / At["a"].size), trim_points_fraction_B=float(Bt["n_points_trimmed"] / Bt["a"].size),
        trim_rel_change_median_B=float(np.median(np.abs(PtB[sB_] / P0B[sB_] - 1))), trim_rel_change_p95_B=float(np.percentile(np.abs(PtB[sB_] / P0B[sB_] - 1), 95)),
        trim_rel_change_median_A=float(np.median(np.abs(PtA[sA_] / P0A[sA_] - 1))), trim_rel_change_p95_A=float(np.percentile(np.abs(PtA[sA_] / P0A[sA_] - 1), 95)),
        despike_bins_fraction_B=stB["bins_replaced_fraction"], despike_points_fraction_B=stB["points_with_replaced_bins"],
        despike_bins_fraction_A=stA["bins_replaced_fraction"], despike_points_fraction_A=stA["points_with_replaced_bins"],
        despike_rel_change_median_B=float(np.median(np.abs(PB[sB_] / PtB[sB_] - 1))), despike_rel_change_p95_B=float(np.percentile(np.abs(PB[sB_] / PtB[sB_] - 1), 95)),
        despike_rel_change_median_A=float(np.median(np.abs(PA[sA_] / PtA[sA_] - 1))), despike_rel_change_p95_A=float(np.percentile(np.abs(PA[sA_] / PtA[sA_] - 1), 95)))

    # drive-power calibration from `full`
    modA, modB = fit_drive_model(Af), fit_drive_model(Bf)
    PdrA, fitA = drive_power_nW(A, Af, modA)
    PdrB, fitB = drive_power_nW(B, Bf, modB)
    yA, yB = to_dBm(PdrA), to_dBm(PdrB)
    eB, cfB = modB[0], modB[1]
    fac = ((Bf["sa_freq"][0, -1] - Bf["sa_freq"][0, 0]) / (Bf["sa_freq"].shape[1] - 1)) / (ENBW_FACTOR * RBW_FULL_HZ)

    def yofa(a):
        return to_dBm(cfB[F_EXC] * np.asarray(a, float) ** eB * fac)

    def aofy(y):
        return (10 ** (np.asarray(y, float) / 10) * 1e6 / fac / cfB[F_EXC]) ** (1.0 / eB)

    cA = np.array(list(modA[1].values())); cB = np.array(list(modB[1].values()))
    num.update(calib_exponent_A=modA[0], calib_exponent_B=modB[0], calib_rms_lnres_A=modA[2], calib_rms_lnres_B=modB[2],
               calib_c_range_A_nW_per_mV2=[float(cA.min()), float(cA.max())], calib_c_range_B_nW_per_mV2=[float(cB.min()), float(cB.max())],
               calib_c_ripple_A_dB=float(10 * np.log10(cA.max() / cA.min())), calib_c_ripple_B_dB=float(10 * np.log10(cB.max() / cB.min())),
               calib_points_from_fit_A=int(fitA), calib_points_from_fit_B=int(fitB),
               calib_enbw_factor=ENBW_FACTOR, calib_bin_over_enbw=float(fac),
               drive_power_range_A_dBm=[float(yA.min()), float(yA.max())], drive_power_range_B_dBm=[float(yB[B['a'] >= A_MIN_ZOOM].min()), float(yB.max())],
               drive_dBm_at_250_300_350mV_4065=[float(yofa(x)) for x in (250, 300, 350)])
    p.step("analysis_drive_calibration", params={"source": "channel `full`, ln SigmaP = ln c(f) + e ln a",
                                                 "ENBW_factor": ENBW_FACTOR, "RBW_full_Hz": RBW_FULL_HZ},
           result={"exponent_A": modA[0], "exponent_B": modB[0],
                   "points_from_fit_A": int(fitA), "points_from_fit_B": int(fitB)},
           note="every f/2 point gets the measured `full` drive power at the same (a, f); fit c(f) a^e for the rest")

    # maps in calibrated drive power
    fgA = np.arange(3760.0, 4251.0, 1.0)
    facA = ((Af["sa_freq"][0, -1] - Af["sa_freq"][0, 0]) / (Af["sa_freq"].shape[1] - 1)) / (ENBW_FACTOR * RBW_FULL_HZ)
    ylo_A = max(float(to_dBm(c * 200.0 ** modA[0] * facA)) for c in modA[1].values())
    yhi_A = min(float(to_dBm(c * 500.0 ** modA[0] * facA)) for c in modA[1].values())
    ygA = np.arange(ylo_A, yhi_A + Y_STEP, Y_STEP)
    sigA_f, sigA_y = SIG_A_CANDIDATES[-1]
    for sf, sy in SIG_A_CANDIDATES:
        ZA, vfrac = smooth_map(yA, A["f"], PA, fgA, ygA, sf, sy)
        num[f"valid_fraction_A_sig{sf}x{sy}"] = vfrac
        sigA_f, sigA_y = sf, sy
        if vfrac >= 0.99:
            break
    neffA = support_neff(yA, A["f"], fgA, ygA, sigA_f, sigA_y)
    num.update(sigma_smooth_A_MHz=sigA_f, sigma_smooth_A_dB=sigA_y, valid_fraction_A=vfrac,
               neff_A_median=float(np.nanmedian(neffA[np.isfinite(ZA)])), neff_A_frac_below2=float(np.nanmean(neffA[np.isfinite(ZA)] < 2)))
    selB = B["a"] >= A_MIN_ZOOM
    fgB = np.arange(B["f"].min(), B["f"].max() + 1, 1.0)
    ylo_B = max(float(to_dBm(c * A_MIN_ZOOM ** modB[0] * fac)) for c in modB[1].values())
    yhi_B = min(float(to_dBm(c * 350.0 ** modB[0] * fac)) for c in modB[1].values())
    ygB = np.arange(ylo_B, yhi_B + Y_STEP_B, Y_STEP_B)
    ZB, vfracB = smooth_map(yB[selB], B["f"][selB], PB[selB], fgB, ygB, *SIG_ZOOM)
    num.update(sigma_smooth_B_MHz=SIG_ZOOM[0], sigma_smooth_B_dB=SIG_ZOOM[1], valid_fraction_B=vfracB,
               native_coverage_B=float(selB.sum() / (fgB.size * (B["a"].max() - A_MIN_ZOOM + 1))))
    p.step("analysis_smooth_maps", params={"kernel": "normalised anisotropic Gaussian in (f, P_dr), linear power",
                                           "sigma_overview": [sigA_f, sigA_y], "sigma_zoom": list(SIG_ZOOM),
                                           "y_axis": "calibrated drive power P_dr (dBm)"},
           result={"grid_A": list(ZA.shape), "grid_B": list(ZB.shape), "valid_fraction_A": vfrac, "valid_fraction_B": vfracB},
           note="cells > 2.5 sigma away from every point stay NaN (white via cmap.set_bad)")
    p.array("map_overview_ZA_linear_mW", np.round(ZA, 12))
    p.array("map_zoom_ZB_linear_mW", np.round(ZB, 12))
    tcur = threshold_curve(ZA, ygA, thrA)
    ts = np.convolve(np.nan_to_num(tcur, nan=np.nanmax(tcur)), np.ones(15) / 15, mode="same")
    mins = [i for i in argrelmin(ts, order=12)[0] if 10 < i < fgA.size - 10]
    num["landscape_threshold_minima_f_MHz_dBm"] = [[float(fgA[i]), float(tcur[i])] for i in mins]
    num["landscape_threshold_at_4065_dBm_overview"] = float(tcur[int(F_EXC - fgA[0])])

    # frequency slice
    band = np.abs(A["a"] - A_SLICE_MV) <= SLICE_BAND_MV
    exact = A["a"][band] == A_SLICE_MV
    fS, PS = A["f"][band], PA[band]
    fg_s = np.arange(fgA[0], fgA[-1] + 1.0)
    Ksl = np.exp(-0.5 * ((fg_s[:, None] - fS[None, :]) / SIG_SLICE_MHZ) ** 2)
    ksum = Ksl.sum(1)
    sm_s = np.where(ksum >= 0.08, (Ksl * PS[None, :]).sum(1) / np.maximum(ksum, 1e-12), np.nan)
    ipk = int(np.nanargmax(sm_s)); pk_val = float(sm_s[ipk])
    o_ = np.argsort(fS)
    fs_o, ps_o = fS[o_], PS[o_]
    cand = [(ps_o[i + 1] - ps_o[i], fs_o[i], fs_o[i + 1], ps_o[i], ps_o[i + 1]) for i in range(fs_o.size - 1) if fs_o[i + 1] - fs_o[i] <= 8.0]
    cand = sorted(cand, reverse=True)
    jumps = []
    for c_ in cand:
        if all(abs(c_[1] - j[1]) > 5.0 for j in jumps):
            jumps.append(c_)
        if len(jumps) == 2:
            break
    jumps = sorted(jumps, key=lambda j: j[1])
    num.update(slice_amp_mV=A_SLICE_MV, slice_band_mV=SLICE_BAND_MV, slice_n_band=int(band.sum()), slice_n_exact=int(exact.sum()),
               slice_n_above_3sigma=int((PS > thrA).sum()), slice_peak_f_MHz=float(fg_s[np.nanargmax(sm_s)]), slice_peak_smoothed_nW=pk_val,
               slice_drive_dBm_median=float(np.median(yA[band])), slice_median_spacing_MHz=float(np.median(np.diff(fs_o))),
               slice_jumps=[dict(df_MHz=float(j[2] - j[1]), dP_nW=float(j[0]), f_lo=float(j[1]), f_hi=float(j[2]), P_lo=float(j[3]), P_hi=float(j[4])) for j in jumps],
               slice_plateau_between_jumps_nW=(float(np.mean([p_ for f_, p_ in zip(fs_o, ps_o) if jumps[0][2] < f_ < jumps[1][1]])) if len(jumps) == 2 and any(jumps[0][2] < f_ < jumps[1][1] for f_ in fs_o) else None))

    # cut at f_exc
    a_cut, P_cut_all, flag = cut_at_frequency(B, F_EXC, PB)
    sem_lookup = {float(a): float(sm) / NW for a, sm, ff in zip(B["a"], B["pint_sem"], B["f"]) if ff == F_EXC}
    sem_cut = np.array([sem_lookup.get(float(a), np.nan) for a in a_cut])
    kfac, npairs = group_scatter_factor(FILE_ZOOM, F_EXC)
    num.update(sem_inflation_factor=kfac, sem_inflation_n_pairs=npairs)
    okmask = (flag == 0) & (a_cut >= A_MIN_ZOOM)
    a_c, P_c = a_cut[okmask], P_cut_all[okmask]
    s_c = np.where(np.isfinite(sem_cut[okmask]), kfac * sem_cut[okmask], np.nan)
    num["n_zoom_points_below_A_MIN_excluded"] = int(((B["a"] > 0) & (B["a"] < A_MIN_ZOOM)).sum())
    a_loo, P_loo, P_meas_loo = loo_cut(B, F_EXC, PB)
    k_ = a_loo >= A_MIN_ZOOM
    a_loo, P_loo, P_meas_loo = a_loo[k_], P_loo[k_], P_meas_loo[k_]
    a_th = threshold_amplitude(a_c, P_c, thrB)
    a_fine = np.arange(a_c.min(), a_c.max() + 0.5, 0.5)
    P_sm = smooth_cut(a_c, P_c, a_fine, SIG_CUT_MV)
    strong = P_meas_loo > 10 * sigB
    rel = np.abs(P_loo[strong] / P_meas_loo[strong] - 1)
    num.update(n_cut_points=int(a_c.size), n_interpolated_at_fexc=int((flag == 1).sum()),
               n_missing_at_fexc=int((flag == 2).sum()), a_threshold_3sigma_mV=float(a_th),
               a_threshold_3sigma_dBm=float(yofa(a_th)), n_loo=int(a_loo.size),
               loo_median_rel_err=float(np.median(rel)), loo_p90_rel_err=float(np.percentile(rel, 90)),
               loo_max_rel_err=float(rel.max()), loo_n_strong=int(strong.sum()),
               loo_rms_abs_err_all_nW=float(np.sqrt(np.mean((P_loo - P_meas_loo) ** 2))))
    num["loo_rms_abs_err_over_sigma"] = num["loo_rms_abs_err_all_nW"] / sigB
    p.array("cut_fexc_a_mV", a_c)
    p.array("cut_fexc_SigmaP_nW", P_c)

    # slope-change analysis
    synthetic_hinge_test()
    h1 = hinge_fit(a_c, P_c, *K1_WINDOW); h2 = hinge_fit(a_c, P_c, *K2_WINDOW)
    ho = hinge_fit(a_c, P_c, 250.0, 285.0, margin=4.0)
    p.step("analysis_hinge", params={"windows_mV": [K1_WINDOW, K2_WINDOW], "model": "two straight segments, hinge scanned on 0.5 mV grid"},
           note="broken-stick fits, bootstrap CIs, LOO cross-validation and wild-bootstrap p-values follow")
    b1 = hinge_bootstrap(a_c, P_c, *K1_WINDOW); b2 = hinge_bootstrap(a_c, P_c, *K2_WINDOW)
    cv1 = cv_compare(a_c, P_c, *K1_WINDOW)
    cv1b = cv_compare(a_c, P_c, 295.0, 340.0); cv1c = cv_compare(a_c, P_c, 280.0, 350.0)
    p3, f95_3, fobs = wild_bootstrap_p(a_c, P_c, *K1_WINDOW, deg=3)
    p4, f95_4, _ = wild_bootstrap_p(a_c, P_c, *K1_WINDOW, deg=4)
    p5, f95_5, _ = wild_bootstrap_p(a_c, P_c, *K1_WINDOW, deg=5)
    num.update(
        onset_hinge_mV=ho["ak"],
        k1_mV=h1["ak"], k1_ci68_mV=[float(b1[0, 0]), float(b1[2, 0])], k1_dBm=float(yofa(h1["ak"])),
        k1_slope_below=h1["s1"], k1_slope_above=h1["s2"],
        k1_slope_below_ci68=[float(b1[0, 1]), float(b1[2, 1])], k1_slope_above_ci68=[float(b1[0, 2]), float(b1[2, 2])],
        k1_F=h1["F"], k1_slope_ratio=h1["s2"] / h1["s1"],
        k1_cv_rmse_285_335=cv1, k1_cv_rmse_295_340=cv1b, k1_cv_rmse_280_350=cv1c,
        k1_wild_p_cubic=p3, k1_wild_p_quartic=p4, k1_wild_p_quintic=p5,
        k2_mV=h2["ak"], k2_ci68_mV=[float(b2[0, 0]), float(b2[2, 0])], k2_dBm=float(yofa(h2["ak"])),
        k2_slope_below=h2["s1"], k2_slope_above=h2["s2"], k2_F=h2["F"], k2_slope_ratio=h2["s2"] / h2["s1"])
    slope, slope_se = local_slope(a_c, P_c)
    for lab, d_, P_ in (("raw", B0, P0B), ("trim", Bt, PtB)):
        _, Pcut_, _ = cut_at_frequency(d_, F_EXC, P_)
        hh_ = hinge_fit(a_c, Pcut_[okmask], *K1_WINDOW)
        num[f"k1_{lab}_mV"] = hh_["ak"]; num[f"k1_{lab}_F"] = hh_["F"]; num[f"k1_{lab}_slopes"] = [hh_["s1"], hh_["s2"]]
    lnP_c = np.log(np.maximum(P_c, 1.0))
    slope_ln, slope_ln_se = local_slope(a_c, lnP_c)
    h1l = hinge_fit(a_c, lnP_c, *K1_WINDOW)
    cv1l = cv_compare(a_c, lnP_c, *K1_WINDOW)
    pl3, _, _ = wild_bootstrap_p(a_c, lnP_c, *K1_WINDOW, deg=3)
    pl4, _, _ = wild_bootstrap_p(a_c, lnP_c, *K1_WINDOW, deg=4)
    w_ln = (a_c >= 290) & (a_c <= 345)
    base = (a_c >= 295) & (a_c <= 312)
    segm = (a_c >= 300) & (a_c <= 345)
    ipk = np.flatnonzero(segm)[np.nanargmax(slope_ln[segm])]
    sl_base, sl_base_sd = float(np.nanmean(slope_ln[base])), float(np.nanstd(slope_ln[base]))
    num.update(ln_slope_peak_a_mV=float(a_c[ipk]), ln_slope_peak=float(slope_ln[ipk]), ln_slope_peak_se=float(slope_ln_se[ipk]),
               ln_slope_baseline_295_312=sl_base, ln_slope_baseline_sd=sl_base_sd,
               ln_slope_peak_z=float((slope_ln[ipk] - sl_base) / max(slope_ln_se[ipk], 1e-9)))
    num.update(ln_k1_mV=h1l["ak"], ln_k1_F=h1l["F"], ln_wild_p_cubic=pl3, ln_wild_p_quartic=pl4, ln_cv_rmse_285_335=cv1l,
               ln_slope_mean_290_345_per_mV=float(np.nanmean(slope_ln[w_ln])), ln_slope_min_max_290_345=[float(np.nanmin(slope_ln[w_ln])), float(np.nanmax(slope_ln[w_ln]))],
               ln_slope_below_k1_mean=float(np.nanmean(slope_ln[(a_c > 295) & (a_c < h1["ak"] - 3)])),
               ln_slope_above_k1_mean=float(np.nanmean(slope_ln[(a_c > h1["ak"] + 3) & (a_c < 335)])))

    # control: drive chain smoothness
    mfull = (Bf["f"] == F_EXC) & (Bf["a"] >= A_MIN_ZOOM)
    o_ = np.argsort(Bf["a"][mfull]); a_fl = Bf["a"][mfull][o_]
    Pfull = (Bf["spec_mean"][mfull].astype(float).sum(1) / NW)[o_]
    expo = np.full(a_fl.size, np.nan)
    for i, ai in enumerate(a_fl):
        w_ = np.abs(a_fl - ai) <= 5.5
        if w_.sum() >= 6:
            expo[i] = np.polyfit(np.log(a_fl[w_]), np.log(Pfull[w_]), 1)[0]
    hfull = hinge_fit(a_fl, np.log(Pfull), *K1_WINDOW)
    in1 = (a_fl > h1["ak"] - 8) & (a_fl < h1["ak"] + 8)
    c_sq = float(np.sum(Pfull * a_fl ** 2) / np.sum(a_fl ** 4))
    dev_pct = 100.0 * (Pfull / (c_sq * a_fl ** 2) - 1.0)
    num.update(control_full_dev_from_a2_rms_pct=float(dev_pct.std()), control_full_dev_from_a2_maxabs_pct=float(np.abs(dev_pct).max()),
               control_full_dev_near_k1_rms_pct=float(dev_pct[in1].std()))
    num.update(control_full_exponent_4065=float(np.polyfit(np.log(a_fl), np.log(Pfull), 1)[0]),
               control_full_rms_lnres_4065=float(np.std(np.log(Pfull) - np.polyval(np.polyfit(np.log(a_fl), np.log(Pfull), 1), np.log(a_fl)))),
               control_full_local_exponent_around_k1=[float(np.nanmin(expo[in1])), float(np.nanmax(expo[in1]))],
               control_full_hinge_F_lnP=hfull["F"], control_full_hinge_ak_lnP=hfull["ak"],
               control_full_exp_min_290_345=float(np.nanmin(expo[(a_fl >= 290) & (a_fl <= 345)])), control_full_exp_max_290_345=float(np.nanmax(expo[(a_fl >= 290) & (a_fl <= 345)])))

    # slope-change line a_k(f)
    rows = []
    for f0 in fgB:
        sel = (B["f"] == f0) & (B["a"] >= A_MIN_ZOOM)
        if sel.sum() < 30:
            continue
        o = np.argsort(B["a"][sel])
        af, Pf_ = B["a"][sel][o], PB[sel][o]
        h = hinge_fit(af, Pf_, *K1_WINDOW)
        good = (h["F"] >= 30) and (h["s2"] >= 2.5 * max(h["s1"], 0.3)) and (K1_WINDOW[0] + 8 < h["ak"] < K1_WINDOW[1] - 8)
        cvw = False
        cvwl = False
        if good:
            cv = cv_compare(af, Pf_, *K1_WINDOW)
            cvw = cv["hinge"] < min(cv["poly2"], cv["poly3"], cv["poly4"])
            cvl_ = cv_compare(af, np.log(np.maximum(Pf_, 1.0)), *K1_WINDOW)
            cvwl = cvl_["hinge"] < min(cvl_["poly2"], cvl_["poly3"], cvl_["poly4"])
        rows.append((f0, h["ak"], h["s1"], h["s2"], h["F"], good, cvw, cvwl))
    rows = np.array(rows, dtype=float)
    gk = rows[:, 5] > 0
    fk, ak_f, cvwin = rows[gk, 0], rows[gk, 1], rows[gk, 6] > 0
    kfit = np.polyfit(fk, ak_f, 1)
    num.update(kink_line_f_range_MHz=[float(fk.min()), float(fk.max())], kink_line_n_rows=int(gk.sum()),
               kink_line_n_rows_cv_supported=int(cvwin.sum()), kink_line_n_rows_cv_supported_ln=int((rows[gk, 7] > 0).sum()),
               kink_line_ak_range_mV=[float(ak_f.min()), float(ak_f.max())], kink_line_slope_mV_per_MHz=float(kfit[0]),
               kink_line_rms_resid_mV=float(np.sqrt(np.mean((np.polyval(kfit, fk) - ak_f) ** 2))))
    ak_lookup = dict(zip(fk, ak_f))

    # spectra at f_exc
    m4065 = B["f"] == F_EXC
    a4065 = B["a"][m4065]; idx4065 = np.flatnonzero(m4065)
    sel_a = [float(a4065[a4065 > 0][np.argmin(np.abs(a4065[a4065 > 0] - t))]) for t in SPEC_TARGETS]
    spec510 = {}
    for ai in sel_a:
        i = idx4065[np.flatnonzero(a4065 == ai)[0]]
        off = (B["sa_freq"][i] - F_EXC / 2 * 1e6) / 1e3
        raw = B["spec_mean"][i].astype(np.float64)
        dens = 1.0 / (raw.sum() * (off[1] - off[0]) / 1e3)
        sm = gaussian_filter1d(raw, SIG_SPEC_BINS, mode="nearest")
        w, c = fwhm_contiguous(off, sm)
        spec510[ai] = (off, raw * dens, sm * dens, w, c, float(PB[i]))
    fw_a = []
    for k in np.flatnonzero(m4065 & (B["a"] >= 280)):
        off = (B["sa_freq"][k] - F_EXC / 2 * 1e6) / 1e3
        w, c = fwhm_contiguous(off, gaussian_filter1d(B["spec_mean"][k].astype(float), SIG_SPEC_BINS, mode="nearest"))
        fw_a.append((B["a"][k], w, c))
    fw_a = np.array(sorted(fw_a))
    bkf, akf = fw_a[:, 0] < h1["ak"] - 3, fw_a[:, 0] > h1["ak"] + 3
    num.update(fwhm510_below_k1_median_kHz=float(np.nanmedian(fw_a[bkf, 1])), fwhm510_above_k1_median_kHz=float(np.nanmedian(fw_a[akf, 1])),
               centroid510_below_k1_median_kHz=float(np.nanmedian(fw_a[bkf, 2])), centroid510_above_k1_median_kHz=float(np.nanmedian(fw_a[akf, 2])),
               spectra_amplitudes_mV=sel_a)
    grp = {"below": [], "above": []}
    for i in np.flatnonzero(np.abs(A["f"] - F_EXC) <= 3.0):
        akk = ak_lookup.get(float(A["f"][i]))
        if akk is None:
            continue
        if akk - 30 <= A["a"][i] <= akk - 5:
            grp["below"].append(i)
        elif akk + 5 <= A["a"][i] <= akk + 30:
            grp["above"].append(i)
    spec20 = {}
    for name, ii in grp.items():
        off0 = (A["sa_freq"][ii[0]] - A["f"][ii[0]] / 2 * 1e6) / 1e3
        for i in ii:
            assert np.allclose((A["sa_freq"][i] - A["f"][i] / 2 * 1e6) / 1e3, off0, atol=1.0)
        mean = np.mean([A["spec_mean"][i].astype(np.float64) for i in ii], axis=0)
        dens = 1.0 / (mean.sum() * (off0[1] - off0[0]) / 1e3)
        sm = gaussian_filter1d(mean, SIG_SPEC20_BINS, mode="nearest")
        w, c = fwhm_contiguous(off0, sm, halfwin=300.0)
        spec20[name] = (off0, mean * dens, sm * dens, w, c, len(ii))
    num.update(rbw20_n_points_below=len(grp["below"]), rbw20_n_points_above=len(grp["above"]),
               rbw20_fwhm_below_kHz=spec20["below"][3], rbw20_fwhm_above_kHz=spec20["above"][3],
               rbw20_centre_below_kHz=spec20["below"][4], rbw20_centre_above_kHz=spec20["above"][4])
    p.step("analysis_complete", params={"note": "all reductions done; figure next"})

    # =========================================================================================
    # Figure: 7.2 in wide (Nature double column), manual inch layout; two new panels at the bottom
    # =========================================================================================
    W, H = 7.2, 10.40
    YOFF = 3.35                      # vertical shift applied to the original reference layout
    fig = plt.figure(figsize=(W, H))

    def axi(x0, y0, w, h):
        return fig.add_axes([x0 / W, y0 / H, w / W, h / H])

    XL, WL = 0.52, 3.90
    XR, WR = 4.80, 2.35
    axA = axi(XL, 4.50 + YOFF, WL, 1.95); axB = axi(XR, 4.50 + YOFF, WR, 1.95)
    axS = axi(XL, 3.40 + YOFF, WL, 1.00)
    axC = axi(XL, 1.22 + YOFF, WL, 1.30); axD = axi(XL, 0.48 + YOFF, WL, 0.70)
    axE = axi(XR, 2.95 + YOFF, WR, 1.35); axF = axi(XR, 0.93 + YOFF, WR, 1.35)
    axH = axi(XL, 0.58, WL, 1.42)             # new: pulse-timing panel (h)
    axI_ = axi(XR, 0.58, WR, 1.42)            # new: reserved sample-characterisation slot (i)
    cmap = plt.get_cmap("magma").copy(); cmap.set_bad("white")

    def draw_map(ax, fg, yg, Z, thr, info, vmax=None, right=False):
        p_ = ax.get_position()
        cax = fig.add_axes([p_.x0, p_.y1 + 0.085 / H, p_.width, 0.07 / H])
        Zp = np.where(np.isfinite(Z), np.maximum(Z, thr), np.nan)
        im = ax.pcolormesh(fg, yg, np.ma.masked_invalid(Zp), cmap=cmap, norm=LogNorm(vmin=thr, vmax=vmax or np.nanmax(Z)),
                           shading="nearest", rasterized=True)
        cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend="min")
        cax.xaxis.set_major_locator(LogLocator(base=10, numticks=8)); cax.xaxis.set_minor_formatter(NullFormatter())
        cax.xaxis.set_ticks_position("top"); cax.xaxis.set_label_position("top")
        cb.set_label(r"Integrated $f/2$ power $\Sigma P$ (nW)", labelpad=3)
        cax.tick_params(length=1.6, width=0.4, pad=1.0)
        ax.text(0.99 if right else 0.015, 0.03, info, transform=ax.transAxes, color="white", fontsize=5.0,
                ha="right" if right else "left", va="bottom")
        ax.set_xlabel("Drive frequency $f$ (MHz)"); ax.set_ylabel("Drive power $P_{\\rm dr}$ (dBm)")

    # ---- a (overview map) with inset zoom
    draw_map(axA, fgA, ygA, ZA, thrA, f"Gaussian smoothing $\\sigma$ = {sigA_f:g} MHz, {sigA_y:g} dB\nRBW 20 kHz")
    axA.set_xlim(fgA[0] - 0.5, fgA[-1] + 0.5); axA.set_ylim(ygA[0] + 0.25, ygA[-1])
    zf = (fgB.min() - 0.5, F_PLOT_MAX_B + 0.5); zy = (ygB[0], ygB[-1])
    axA.add_patch(Rectangle((zf[0], zy[0]), zf[1] - zf[0], zy[1] - zy[0], fill=False, ec="#7fdcff", lw=0.7))
    axA.plot([F_EXC, F_EXC], [zy[0], zy[1]], color="#7fdcff", lw=0.4, ls=(0, (2, 2)))
    axA.xaxis.set_major_locator(MultipleLocator(100)); axA.xaxis.set_minor_locator(MultipleLocator(25))
    axA.set_xlabel(""); axA.tick_params(labelbottom=False)
    axA.yaxis.set_major_locator(MultipleLocator(2)); axA.yaxis.set_minor_locator(MultipleLocator(1))

    IW, IH, IX, IY = 0.215, 0.425, 0.080, 0.165
    axI = axA.inset_axes([IX, IY, IW, IH])
    vmx = float(np.nanpercentile(ZB[:, fgB <= F_PLOT_MAX_B], 99.5))
    imI = axI.pcolormesh(fgB, ygB, np.ma.masked_invalid(np.maximum(ZB, 0.0)), cmap=cmap, norm=Normalize(vmin=0, vmax=vmx),
                         shading="nearest", rasterized=True)
    axI.set_xlim(*zf); axI.set_ylim(*zy)
    caxI = axA.inset_axes([IX, IY + IH + 0.028, IW, 0.034])
    cbI = fig.colorbar(imI, cax=caxI, orientation="horizontal", extend="max")
    cbI.set_ticks([0, 80, 160])
    caxI.xaxis.set_ticks_position("top"); caxI.xaxis.set_label_position("top")
    caxI.tick_params(length=1.2, width=0.3, pad=0.6, labelsize=3.8, colors="white")
    cbI.outline.set_edgecolor("white"); cbI.outline.set_linewidth(0.3)
    cbI.set_label(r"$\Sigma P$ (nW), linear", fontsize=3.9, color="white", labelpad=0.8)
    axI.axvline(F_EXC, color="#7fdcff", lw=0.4, ls=(0, (3, 2)))
    yk_row = np.array([float(to_dBm(cfB.get(f_, cfB[F_EXC]) * a_ ** eB * fac)) for f_, a_ in zip(fk, ak_f)])
    axI.plot(fk, yk_row, color="white", lw=0.45, ls=(0, (1.2, 1.4)), zorder=4)
    axI.plot(fk[cvwin], yk_row[cvwin], "o", ms=1.2, mfc="white", mec="none", zorder=4)
    axI.plot(fk[~cvwin], yk_row[~cvwin], "o", ms=1.2, mfc="none", mec="white", mew=0.3, zorder=4)
    axI.plot([F_EXC], [yofa(a_th)], "^", ms=2.4, mfc="white", mec="black", mew=0.3, zorder=6)
    axI.plot([F_EXC], [yofa(h1["ak"])], "D", ms=2.5, mfc=C_K1, mec="white", mew=0.4, zorder=6)
    axI.plot([F_EXC], [yofa(h2["ak"])], "D", ms=2.5, mfc=C_K2, mec="white", mew=0.4, zorder=6)
    for ai, c in zip(sel_a, C_BELOW + C_ABOVE):
        axI.plot([F_EXC], [yofa(ai)], "o", ms=2.0, mfc=c, mec="white", mew=0.3, zorder=5)
    axI.set_xticks([4050, 4060, 4070, 4080]); axI.set_yticks([-23, -22, -21])
    axI.tick_params(labelsize=3.9, length=1.4, width=0.3, pad=0.8, colors="white", direction="out")
    for sp in axI.spines.values():
        sp.set_edgecolor("white"); sp.set_linewidth(0.4)
    axI.set_xlabel("$f$ (MHz)", fontsize=4.3, color="white", labelpad=0.8)
    axI.set_ylabel("$P_{\\rm dr}$ (dBm)", fontsize=4.3, color="white", labelpad=0.8)
    axI.text(0.03, 0.97, r"zoom, RBW 510 kHz", transform=axI.transAxes, fontsize=3.7, color="white", va="top", ha="left")
    for yy_, yb_ in ((zy[1], 1.0), (zy[0], 0.0)):
        axA.add_artist(ConnectionPatch(xyA=(zf[0], yy_), coordsA="data", axesA=axA, xyB=(1.0, yb_), coordsB="axes fraction", axesB=axI,
                                       color="#7fdcff", lw=0.35, zorder=3))

    # ---- c: frequency slice
    for sp_ in ("top", "right"):
        axS.spines[sp_].set_visible(False)
    axS.axvspan(zf[0], zf[1], color="#7fdcff", alpha=0.18, lw=0, zorder=0)
    axS.axvline(F_EXC, color="#2a9fc9", lw=0.4, ls=(0, (3, 2)), zorder=1)
    axS.plot(fS[~exact], PS[~exact], ".", ms=2.4, color="0.6", mew=0, zorder=2)
    axS.plot(fS[exact], PS[exact], "o", ms=2.6, mfc="#c0392b", mec="none", zorder=4)
    axS.plot(fg_s, sm_s, "-", color="0.15", lw=0.7, zorder=3)
    for j in jumps:
        axS.plot([j[1], j[2]], [j[3], j[4]], "-", color="#e6550d", lw=1.3, zorder=5, solid_capstyle="round")
        axS.annotate(f"+{j[0]:.0f} nW\nin {j[2] - j[1]:.0f} MHz", xy=(0.5 * (j[1] + j[2]), 0.5 * (j[3] + j[4])), xytext=(-11, 2), textcoords="offset points",
                     fontsize=4.8, color="#e6550d", ha="right", va="center")
    axS.set_xlim(axA.get_xlim()); axS.set_ylim(-0.05 * pk_val, 1.15 * float(np.nanmax(PS)))
    axS.xaxis.set_major_locator(MultipleLocator(100)); axS.xaxis.set_minor_locator(MultipleLocator(25))
    axS.yaxis.set_major_locator(MultipleLocator(100)); axS.yaxis.set_minor_locator(MultipleLocator(50))
    axS.set_xlabel("Drive frequency $f$ (MHz)"); axS.set_ylabel(r"$\Sigma P$ (nW)")
    axS.text(0.012, 0.95, f"slice at $a$ = {A_SLICE_MV:.0f} $\\pm$ {SLICE_BAND_MV:.0f} mV ($P_{{\\rm dr}}\\approx${num['slice_drive_dBm_median']:.1f} dBm), RBW 20 kHz",
             transform=axS.transAxes, fontsize=5.2, color="0.2", va="top", ha="left")
    axS.text(0.995, 0.95, r"$f_{\rm exc}$ and inset window", transform=axS.transAxes, fontsize=4.6, color="#2a9fc9", va="top", ha="right")

    # ---- b: schematic of excitation + measurement
    draw_sequence(axB, WR, 1.95)

    # ---- d / e
    for ax in (axC, axD, axE, axF):
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for h_, col, span in ((h1, C_K1, b1), (h2, C_K2, b2)):
        for ax in (axC, axD):
            ax.axvspan(span[0, 0], span[2, 0], color=col, alpha=0.18, lw=0, zorder=0)
            ax.axvline(h_["ak"], color=col, lw=0.5, zorder=1)
    aa1 = np.linspace(h1["lo"], h1["hi"], 200); aa2 = np.linspace(h2["lo"], h2["hi"], 200)
    axC.axhline(thrB, color="0.7", lw=0.4, ls=":")
    axC.errorbar(a_c, P_c, yerr=np.where(np.isfinite(s_c), s_c, 0), fmt="o", ms=1.9, mfc="#c0392b", mec="none",
                 ecolor="#c0392b", elinewidth=0.35, capsize=0, lw=0, zorder=4)
    axC.plot(a_fine, P_sm, "-", color="0.2", lw=0.55, zorder=3)
    axC.plot(aa1, hinge_model(h1, aa1), "--", color=C_K1, lw=0.8, zorder=5)
    axC.plot(aa2, hinge_model(h2, aa2), "--", color=C_K2, lw=0.8, zorder=5)
    axC.set_xlim(285, 350.5); axC.set_ylim(0, 168)
    axC.set_ylabel(r"$\Sigma P$ (nW)")
    axC.tick_params(labelbottom=False)
    axC.yaxis.set_major_locator(MultipleLocator(40)); axC.yaxis.set_minor_locator(MultipleLocator(20))
    sec = axC.secondary_xaxis("top", functions=(lambda a: yofa(a), lambda y: aofy(y)))
    sec.set_xlabel("Drive power $P_{\\rm dr}$ (dBm)", labelpad=2)
    sec.xaxis.set_major_locator(MultipleLocator(1.0)); sec.xaxis.set_minor_locator(MultipleLocator(0.5))
    sec.tick_params(length=2.0, width=0.5, pad=1.4, labelsize=5.6)
    axC.text(0.015, 0.97, r"$f_{\rm exc}$ = 4065 MHz", transform=axC.transAxes, fontsize=6, va="top", ha="left")
    leg = [Line2D([], [], marker="o", ms=2.2, color="#c0392b", lw=0, label=f"acquired (error bars: {kfac:.1f}$\\times$ s.e.m.)"),
           Line2D([], [], color="0.2", lw=0.55, label=f"Gaussian smoothing, $\\sigma$ = {SIG_CUT_MV:g} mV"),
           Line2D([], [], color=C_K1, ls="--", lw=0.8, label="broken-stick fits")]
    axC.legend(handles=leg, loc="upper left", bbox_to_anchor=(0.0, 0.88), frameon=False, handlelength=1.4,
               borderaxespad=0.0, labelspacing=0.25)
    axC.text(a_c.max() - 0.5, 3, f"3$\\sigma$ = {thrB:.1f} nW", fontsize=4.8, color="0.45", ha="right", va="bottom")

    axD.fill_between(a_c, slope - slope_se, slope + slope_se, color="0.82", lw=0, zorder=1)
    axD.plot(a_c, slope, "-", color="0.15", lw=0.6, zorder=3)
    for x0_, x1_, s_, col in ((K1_WINDOW[0], h1["ak"], h1["s1"], C_K1), (h1["ak"], h2["ak"], h1["s2"], C_K1), (h2["ak"], K2_WINDOW[1], h2["s2"], C_K2)):
        axD.plot([x0_, x1_], [s_, s_], "--", color=col, lw=0.8, zorder=4)
    axD.set_xlim(285, 350.5); axD.set_ylim(0, 7.4)
    axD.set_xlabel("Drive amplitude $a$ (mV)"); axD.set_ylabel(r"d$\Sigma P$/d$a$" + "\n(nW mV$^{-1}$)")
    axD.xaxis.set_major_locator(MultipleLocator(10)); axD.xaxis.set_minor_locator(MultipleLocator(5))
    axD.yaxis.set_major_locator(MultipleLocator(3)); axD.yaxis.set_minor_locator(MultipleLocator(1.5))
    axD.text(0.02, 0.70, f"hinge vs. smooth null: p = {p3:.2f} (linear), {pl3:.2f} (log $\\Sigma P$)", transform=axD.transAxes,
             fontsize=4.9, color="0.3", va="top", ha="left")
    for xx, yy, col, dy_ in ((300.0, h1["s1"], C_K1, 3), (332.5, h1["s2"], C_K1, 3), (348.8, h2["s2"], C_K2, 3)):
        axD.annotate(f"{yy:.2f}", (xx, yy), xytext=(0, dy_), textcoords="offset points", color=col, fontsize=5.2, ha="center")
    for h_, col, txt in ((h1, C_K1, "$a_{k1}$"), (h2, C_K2, "$a_{k2}$")):
        axD.text(h_["ak"] - 1.4, 6.9, f"{txt} = {h_['ak']:.1f} mV", color=col, fontsize=5.2, va="top", ha="right")

    # ---- f / g
    for ai, c in zip(sel_a, C_BELOW + C_ABOVE):
        off, raw_, sm, w, cc, pp = spec510[ai]
        axE.plot(off / 1e3, raw_, lw=0.25, color=c, alpha=0.35)
        axE.plot(off / 1e3, sm, lw=0.8, color=c, label=f"{ai:.0f}")
    axE.set_xlim(-1.5, 1.5); axE.set_ylim(-0.2, 2.5)
    axE.set_xlabel(r"$f_{\rm SA}-f_{\rm exc}/2$ (MHz)"); axE.set_ylabel("Norm. spectrum (MHz$^{-1}$)")
    axE.xaxis.set_major_locator(MultipleLocator(0.5)); axE.xaxis.set_minor_locator(MultipleLocator(0.25))
    axE.yaxis.set_major_locator(MultipleLocator(1)); axE.yaxis.set_minor_locator(MultipleLocator(0.5))
    axE.legend(loc="upper left", frameon=False, handlelength=1.1, borderaxespad=0.0, labelspacing=0.15,
               title="RBW 510 kHz\n$a$ (mV)", title_fontsize=5.0, bbox_to_anchor=(0.0, 1.06))
    ins = axE.inset_axes([0.66, 0.56, 0.33, 0.36])
    for name, col in (("below", "#3182bd"), ("above", "#e6550d")):
        off, raw2, sm2, w2, cc2, n_ = spec20[name]
        ins.plot(off / 1e3, sm2, lw=0.7, color=col)
    ins.set_xlim(-0.3, 0.3); ins.set_ylim(-1, 1.1 * max(spec20[n][2].max() for n in spec20))
    ins.tick_params(labelsize=4.6, length=1.4, pad=0.8, width=0.4)
    for s in ("top", "right"):
        ins.spines[s].set_visible(False)
    ins.xaxis.set_major_locator(MultipleLocator(0.2)); ins.set_yticks([])
    ins.text(0.5, 1.02, "RBW 20 kHz", transform=ins.transAxes, fontsize=4.8, color="0.35", ha="center", va="bottom")
    ins.text(0.98, 0.95, "below $a_{k1}$", transform=ins.transAxes, fontsize=4.6, color="#3182bd", ha="right", va="top")
    ins.text(0.98, 0.78, "above $a_{k1}$", transform=ins.transAxes, fontsize=4.6, color="#e6550d", ha="right", va="top")

    axF.axhline(0.0, color="0.55", lw=0.5)
    for h_, col, span in ((h1, C_K1, b1), (h2, C_K2, b2)):
        axF.axvspan(span[0, 0], span[2, 0], color=col, alpha=0.18, lw=0)
        axF.axvline(h_["ak"], color=col, lw=0.5)
    axF.plot(a_fl, dev_pct, "o-", color="0.15", lw=0.5, ms=1.5)
    axF.set_xlim(285, 350.5); axF.set_ylim(-1.6, 1.6)
    axF.set_xlabel("Drive amplitude $a$ (mV)"); axF.set_ylabel(r"$\Delta P_{\rm dr}/P_{\rm dr}$ (%)")
    axF.xaxis.set_major_locator(MultipleLocator(20)); axF.xaxis.set_minor_locator(MultipleLocator(10))
    axF.yaxis.set_major_locator(MultipleLocator(1)); axF.yaxis.set_minor_locator(MultipleLocator(0.5))
    axF.text(0.02, 0.96, f"deviation from $P_{{\\rm dr}}\\propto a^{{2}}$ (measured: $a^{{{num['control_full_exponent_4065']:.2f}}}$)",
             transform=axF.transAxes, fontsize=4.8, color="0.3", va="top")

    # ---- h (pulse-timing diagram for Figs 3-4) and i (reserved slot)
    draw_timing(axH, WL, 1.42)
    draw_placeholder(axI_, WR, 1.42)

    # ---- panel letters
    for tag, x, y in (("a", 0.03, 6.82 + YOFF), ("b", 4.55, 6.82 + YOFF), ("c", 0.03, 4.43 + YOFF), ("d", 0.03, 2.75 + YOFF),
                      ("e", 0.03, 1.22 + YOFF), ("f", 4.55, 4.36 + YOFF), ("g", 4.55, 2.33 + YOFF),
                      ("h", 0.03, 2.06), ("i", 4.55, 2.06)):
        fig.text(x / W, y / H, tag, fontsize=8, fontweight="bold", va="bottom", ha="left")

    # ---- validation against the raw files
    chk = {}
    rng = np.random.default_rng(0)
    for j in sorted(rng.choice(np.arange(a_c.size), 4, replace=False)):
        P_raw, _ = raw_reduced_point(FILE_ZOOM, "halfaverage", float(a_c[j]), F_EXC, False, None)
        assert np.isclose(P_raw, P_c[j], rtol=1e-4, atol=2e-3), (a_c[j], P_raw, P_c[j])
        chk[f"zoom cut a={a_c[j]:.0f}: pipeline vs independent re-reduction of the raw HDF5 (nW)"] = (float(P_c[j]), float(P_raw))
    for k in sorted(rng.choice(A["a"].size, 3, replace=False)):
        P_raw, _ = raw_reduced_point(FILE_BROAD, "halfaverage", float(A["a"][k]), float(A["f"][k]), False, PROTECT_A_KHZ)
        assert np.isclose(P_raw, PA[k], rtol=1e-4, atol=2e-3), (A["a"][k], A["f"][k], P_raw, PA[k])
        chk[f"broad point a={A['a'][k]:.0f},f={A['f'][k]:.0f}: pipeline vs independent re-reduction (nW)"] = (float(PA[k]), float(P_raw))
    for k in sorted(rng.choice(Bf["a"].size, 2, replace=False)):
        P_raw, _, _ = raw_point_pint(FILE_ZOOM, "full", float(Bf["a"][k]), float(Bf["f"][k]))
        P_c_ = Bf["spec_mean"][k].astype(float).sum() / NW
        assert np.isclose(P_raw, P_c_, rtol=1e-6), (P_raw, P_c_)
        chk[f"full channel a={Bf['a'][k]:.0f},f={Bf['f'][k]:.0f}: cache vs raw HDF5 (nW)"] = (float(P_c_), float(P_raw))
    for ai in (sel_a[0], sel_a[-1]):
        _, sp_raw = raw_reduced_point(FILE_ZOOM, "halfaverage", ai, F_EXC, False, None)
        raw_norm = sp_raw / (sp_raw.sum() * (spec510[ai][0][1] - spec510[ai][0][0]) / 1e3)
        assert np.allclose(raw_norm, spec510[ai][1], rtol=1e-3, atol=1e-4), ai
        chk[f"panel f spectrum a={ai:.0f} mV equals independent re-reduction of the raw HDF5"] = True
    i0 = int(idx4065[0])
    _, _, ax_raw = raw_point_pint(FILE_ZOOM, "halfaverage", float(B["a"][i0]), F_EXC)
    assert np.allclose(ax_raw, B["sa_freq"][i0]), "SA axis differs from raw file"
    assert abs(ax_raw.mean() - F_EXC / 2 * 1e6) < 1e3 and np.all(np.diff(ax_raw) > 0)
    chk["SA axis (halfaverage): centre minus f/2 (kHz)"] = float((ax_raw.mean() - F_EXC / 2 * 1e6) / 1e3)
    chk["SA axis (halfaverage): span (MHz)"] = float((ax_raw[-1] - ax_raw[0]) / 1e6)
    chk["SA axis (full): centre minus f (kHz), max over points"] = float(np.abs(Bf["sa_freq"].mean(1) / 1e3 - Bf["f"] * 1e3).max())
    for (y_, d_, P_, fg_, yg_, nm, m_) in ((yA, A, PA, fgA, ygA, "A", np.ones(A["a"].size, bool)), (yB, B, PB, fgB, ygB, "B", selB)):
        Zt, _ = smooth_map(y_[m_], d_["f"][m_], P_[m_], fg_, yg_, 0.02, 0.002, rel_cut=1.0)
        inside = m_ & (y_ >= yg_[0]) & (y_ <= yg_[-1])
        ii = np.flatnonzero(inside)
        rj = np.round((y_[ii] - yg_[0]) / (yg_[1] - yg_[0])).astype(int); cj = np.round(d_["f"][ii] - fg_[0]).astype(int)
        cell = {}
        for r_, c_, v_ in zip(rj, cj, P_[ii]):
            cell.setdefault((r_, c_), []).append(v_)
        pick = ii[::max(1, ii.size // 200)]
        rp = np.round((y_[pick] - yg_[0]) / (yg_[1] - yg_[0])).astype(int)
        cp = np.round(d_["f"][pick] - fg_[0]).astype(int)
        got = Zt[rp, cp]
        exp_ = np.array([np.mean(cell[(r_, c_)]) for r_, c_ in zip(rp, cp)])
        close = np.isclose(got, exp_, rtol=1e-3)
        assert close.mean() > 0.99, f"map {nm} transposed / mis-indexed ({close.mean():.2f})"
        chk[f"map_{nm}: cell values equal the raw points falling in the cell (rows=P_dr, cols=f), fraction"] = float(close.mean())
    chk["mean_A_raw_vs_smoothed_nW"] = [float(PA.mean()), float(np.nanmean(ZA))]
    chk["synthetic_hinge_test_passed"] = True
    _, v_x, fl_x = cut_at_frequency(B, B["f"].min() - 1.0, PB)
    assert np.all(fl_x == 2) and np.all(np.isnan(v_x)), "extrapolation happened"
    chk["interp_no_extrapolation_ok"] = True
    for ai in sel_a:
        assert ai in set(a4065.tolist()), ai
    chk["panel_f_spectrum_markers_on_acquired_points"] = True
    chk["drive power (dBm) is monotonic in a at f_exc"] = bool(np.all(np.diff(yofa(np.arange(250.0, 351.0))) > 0))
    num["validation"] = chk
    p.step("validation", params={"n_checks": len(chk)}, note="independent re-reads of raw HDF5 groups; all asserts passed",
           result={"n_checks": len(chk), "n_fail": 0})

    # numbers JSON + figure output
    with open(FIG_ROOT / "manifests" / "fig1_overview_numbers.json", "w") as fh:
        json.dump(num, fh, indent=2, default=float)
    p.step("numbers_json", params={"path": str(FIG_ROOT / "manifests" / "fig1_overview_numbers.json")},
           note="ported numbers dump of the reference script")
    p.finish(fig, ["a overview map + zoom inset (single tone)", "b excitation/measurement schematic",
                   "c frequency slice", "d cut at f_exc", "e local slope", "f f/2 spectra",
                   "g drive-chain control", "h two-tone pulse timing (Figs 3-4)", "i reserved: sample characterisation"])


if __name__ == "__main__":
    reextract = "--reextract" in sys.argv
    for fp in (FILE_BROAD, FILE_ZOOM):
        if not fp.exists():
            sys.exit(f"missing {fp}")
    build(reextract=reextract)
