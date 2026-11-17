import numpy as np
P0=np.load("P_on_peak.npy"); P4=np.load("P_off_peak.npy")
PHI_ON, PHI_OFF, PHI0_REF = 1.310, 4.550, np.deg2rad(-155.0)
mask = (P0>1e-6)|(P4>1e-6)
i,j = np.where(mask)
m = (P0[i,j]+P4[i,j])/2.0
u1 = P0[i,j]/m; u2 = P4[i,j]/m
n = len(u1)
ths = np.linspace(0,2*np.pi,7201)
c1 = np.cos(PHI_ON-ths); c2 = np.cos(PHI_OFF-ths)
S1, S2 = (u1-1).sum(), (u2-1).sum()
def fit_curve(num):
    g = (c1*num[0]+c2*num[1])/(n*(c1**2+c2**2))
    g = np.where(g>0, g, np.nan)   # A1>0 branch
    ss = n*((u1.mean()-1-g*c1)**2+(u2.mean()-1-g*c2)**2)*0  # placeholder
    # actual: use group means only for phi0 landscape (scatter is group-independent)
    chi = ((u1.mean()-1-g*c1)**2/(u1.var(ddof=1)/n) + (u2.mean()-1-g*c2)**2/(u2.var(ddof=1)/n))
    return g, chi
g, chi = fit_curve((S1,S2))
k = np.nanargmin(chi)
phi0 = ths[k]; print("phi0 = %.1f deg, g=%.4f" % (np.degrees(phi0), g[k]))
# delta-SS CI (chi2+2.3 for 1 dof 68%? use chi2 profile: chi-2.30 for 68% 2-param... phi0 1 param -> 1.0)
ok = chi - chi[k] < 1.0
okd = np.degrees(ths[ok])
print("phi0 68%% CI from mean-level chi2 profile: [%.1f, %.1f] deg (contiguous ranges:)" % (okd.min(), okd.max()))
# segmentation (mod 360)
d = np.diff(okd); brk = np.where(d>1.0)[0]
print("  ranges:", np.round(np.r_[okd[0], np.where(brk)[0], ],1) if len(brk) else "single [%.1f..%.1f]"%(okd.min(),okd.max()))
# bootstrap over cells
rng = np.random.RandomState(7)
bt=[]
m1 = u1.mean()-1; m2 = u2.mean()-1
v1 = u1.var(ddof=1)/n; v2 = u2.var(ddof=1)/n
for _ in range(2000):
    d1 = rng.normal(m1, np.sqrt(v1)); d2 = rng.normal(m2, np.sqrt(v2))
    g_ = (c1*d1+c2*d2)/(c1**2+c2**2)
    chi_ = (d1-g_*c1)**2/v1 + (d2-g_*c2)**2/v2
    kk = np.nanargmin(chi_)
    bt.append(np.degrees(ths[kk]))
bt=np.array(bt)
med = np.median(bt); lo,hi=np.percentile(bt,[16,84])
print("bootstrap phi0: median %.1f deg, 68%% CI [%.1f, %.1f]" % (med, lo, hi))
# cell-level predictive test at phi0=-155
cr = np.cos(PHI_OFF-PHI0_REF) - np.cos(PHI_ON-PHI0_REF)
d = u2-u1
frac = (d>0).mean()
from scipy import stats
print("fixed phi0=-155deg predicts OFF>ON in every cell (A1>0).")
blo, bhi = stats.beta.ppf([0.025, 0.975], frac*n+0.5, (1-frac)*n+0.5)
print("observed OFF>ON fraction: %.3f (binomial 95%% CI: %.3f..%.3f)" % (frac, blo, bhi))
print("mean OFF-ON normalized difference: %+.4f +- %.4f (SEM), t=%.1f" % (d.mean(), d.std(ddof=1)/np.sqrt(n), d.mean()/(d.std(ddof=1)/np.sqrt(n))))
# residual structure vs a2/f2 (strong cells)
res = np.r_[u1-(1+ (u1.mean()-1)), u2-(1+(u2.mean()-1))]
aa = np.r_[np.arange(200,276)[i]*1.0, np.arange(200,276)[i]*1.0]
ff = np.r_[np.arange(4075,4151)[j]*1.0]*1.0+0; ff = np.r_[np.arange(4075,4151)[j], np.arange(4075,4151)[j]].astype(float)
r1 = u1-u1.mean(); r2_ = u2-u2.mean(); rr = np.r_[r1, r2_]
for nm, v in [("a2",aa),("f2",np.r_[np.arange(4075,4151)[j],np.arange(4075,4151)[j]].astype(float))]:
    pr = stats.pearsonr(v, rr)
    print("residual-vs-%s Pearson r = %+.3f (p=%.2g)" % (nm, pr[0], pr[1]))
