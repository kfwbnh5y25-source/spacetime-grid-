import h5py, numpy as np, json
from scipy.signal import medfilt
from scipy.optimize import curve_fit

OUT = 'reports/mistral/tmp/'
amps = np.arange(200, 1001, 10)   # mV
ls   = np.arange(20, 600, 10)     # AWG samples @64 GSa/s
Pwin  = np.full((amps.size, ls.size), np.nan)
Pkpos = np.full_like(Pwin, np.nan)   # line center MHz
Fwhm  = np.full_like(Pwin, np.nan)   # MHz
Bg    = np.full_like(Pwin, np.nan)   # broadband background (median filtered, W)
Snrint= np.full_like(Pwin, np.nan)   # integrated SNR Pwin/Pnoise
Pnoise_perbin = np.full_like(Pwin, np.nan)
freqs = None
with h5py.File('phase gating/3/3.h5','r') as f:
    for key in f:
        try:
            a_s, l_s = key.split('_'); a, l = int(a_s), int(l_s)
        except ValueError:
            continue
        if l < 20 or l % 10 != 0:  # stray short-probe keys
            continue
        ia = np.where(amps == a)[0]; il = np.where(ls == l)[0]
        if ia.size != 1 or il.size != 1: continue
        tr = f[key]['trace'][0]
        if freqs is None: freqs = tr[0].copy()
        spec = tr[1].astype(np.float64)
        fm = medfilt(spec, 15)
        bg = np.median(fm)
        pk = np.argmax(fm); h = fm[pk] - bg
        # adaptive window: expand at h/4 level
        thr4 = bg + 0.25*h; thr2 = bg + 0.5*h
        lo = pk
        while lo > 0 and fm[lo-1] > thr4: lo -= 1
        hi = pk
        while hi < len(fm)-1 and fm[hi+1] > thr4: hi += 1
        # FWHM
        lw = pk
        while lw > 0 and fm[lw-1] > thr2: lw -= 1
        hw = pk
        while hw < len(fm)-1 and fm[hw+1] > thr2: hw += 1
        seg = spec[lo:hi+1]
        sigma = 1.4826*np.median(np.abs(spec - np.median(spec)))
        Pwin[ia[0], il[0]] = seg.sum()
        Pkpos[ia[0], il[0]] = freqs[pk]/1e6
        Fwhm[ia[0], il[0]] = (hw-lw)*0.075
        Bg[ia[0], il[0]] = bg
        Snrint[ia[0], il[0]] = seg.sum()/(sigma*np.sqrt(len(seg))) if len(seg)>0 else np.nan
        Pnoise_perbin[ia[0], il[0]] = sigma

np.save(OUT+'A_Pwin.npy', Pwin); np.save(OUT+'A_Pkpos.npy', Pkpos)
np.save(OUT+'A_fwhm.npy', Fwhm); np.save(OUT+'A_Bg.npy', Bg)
np.save(OUT+'A_Snrint.npy', Snrint)
np.save(OUT+'A_amps.npy', amps); np.save(OUT+'A_ls.npy', ls); np.save(OUT+'A_freqs.npy', freqs)

withnp = np.isfinite(Pwin)
print('filled:', withnp.sum(), 'of', Pwin.size)
# threshold map and a_th(l)
sigma_map = Pnoise_perbin
width_map = np.where(withnp, 1, 1)  # width not stored; recompute approx? store instead
# recompute width by rerun? store Pnoise already; use Snrint>5 as detection (5-sigma integrated)
det = Snrint > 5
a_th = np.full(ls.size, np.nan)
for il in range(ls.size):
    idx = np.where(det[:, il])[0]
    if idx.size: a_th[il] = amps[idx[0]]
np.save(OUT+'A_ath.npy', a_th)
print('a_th(l):')
for il in range(ls.size):
    l_ns = ls[il]/64
    print(f'l={ls[il]:4d} ({l_ns:6.3f} ns): a_th = {a_th[il]:.0f} mV' if np.isfinite(a_th[il]) else f'l={ls[il]:4d} ({l_ns:6.3f} ns): none')
