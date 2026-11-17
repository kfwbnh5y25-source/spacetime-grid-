import numpy as np
P0=np.load("P_on_peak.npy"); P4=np.load("P_off_peak.npy")
A2=list(range(200,276)); F2=list(range(4075,4151))
PHI_ON, PHI_OFF, PHI0_REF = 1.310, 4.550, np.deg2rad(-155.0)
mask = (P0>5e-7)|(P4>5e-7)
i,j = np.where(mask)
m = (P0[i,j]+P4[i,j])/2.0
u_on = P0[i,j]/m; u_off = P4[i,j]/m
print("cells:", len(m))
print("mean u_on %.4f (sem %.4f), mean u_off %.4f (sem %.4f)" % (u_on.mean(), u_on.std(ddof=1)/np.sqrt(len(u_on)), u_off.mean(), u_off.std(ddof=1)/np.sqrt(len(u_off))))
r = u_off.mean()-1.0; q = u_on.mean()-1.0
# solve r/q = cos(PHI_OFF-ph0)/cos(PHI_ON-ph0) for ph0 in [0,2pi)
ths = np.linspace(0, 2*np.pi, 200001)
ratio = np.cos(PHI_OFF-ths)/np.cos(PHI_ON-ths)
ok = np.where(np.abs(ratio - r/q) < 5e-4)[0]
for t in ths[ok]:
    print("phi0_cand = %.3f rad = %.1f deg; g_on=q/cos(PHI_ON-t)=%.3f" % (t, np.degrees(t), q/np.cos(PHI_ON-t)))
# robust: scan phi0, least squares with shared g
best=None
for t in ths[::2]:
    c_on, c_off = np.cos(PHI_ON-t), np.cos(PHI_OFF-t)
    # design: g shared, fit to u_on-1, u_off-1, all cells
    num = (c_on*(u_on-1) + c_off*(u_off-1)).sum()
    den = (c_on**2*u_on.size + c_off**2*u_off.size)
    g = num/den
    pred_on = 1 + g*c_on; pred_off = 1 + g*c_off
    ss = (((u_on-pred_on)**2).sum() + ((u_off-pred_off)**2).sum())
    tot = (((u_on-u_on.mean())**2).sum() + ((u_off-u_off.mean())**2).sum())
    r2 = 1-ss/tot
    if best is None or r2>best[0]: best=(r2,t,g,ss)
r2,t,g,ss = best
print("BEST: R2=%.4f phi0=%.3f rad (%.1f deg) g=%.4f" % (r2,t,np.degrees(t),g))
print("phi0 ref -155deg = %.3f rad; delta = %.1f deg" % (PHI0_REF, np.degrees(t)-(-155.0)))
# residual structure vs a2, f2
pred_on = 1+g*np.cos(PHI_ON-t); pred_off = 1+g*np.cos(PHI_OFF-t)
res_on = u_on-pred_on; res_off = u_off-pred_off
print("residual rms: on %.4f off %.4f" % (res_on.std(), res_off.std()))
a2g = np.array(A2)[i]; f2g = np.array(F2)[j]
res = np.r_[res_on, res_off]; avv = np.r_[a2g, a2g]; fvv = np.r_[f2g, f2g]
for nm,v in [("a2",avv),("f2",fvv)]:
    bins = np.quantile(v, np.linspace(0,1,7))
    means = [res[(v>=bins[k])&(v<=bins[k+1])].mean() for k in range(6)]
    print("residual mean vs %s sextiles:" % nm, np.round(means,4))
# bootstrap phi0
rng = np.random.RandomState(42)
n=len(u_on); bt=[]
for _ in range(200):
    idx = rng.randint(0,n,n)
    uo, ui = u_off[idx], u_on[idx]
    bb=None
    for t2 in ths[::4]:
        c1,c2 = np.cos(PHI_ON-t2), np.cos(PHI_OFF-t2)
        gg = ((c1*(ui-1)).sum()+(c2*(uo-1)).sum())/(c1*c1*n+c2*c2*n)
        s = (((ui-1-gg*c1)**2).sum()+((uo-1-gg*c2)**2).sum())
        if bb is None or s<bb[0]: bb=(s,t2)
    bt.append(bb[1])
bt=np.array(bt); btdeg=np.degrees(bt)
lo,hi = np.percentile(btdeg,[16,84]); md=np.median(btdeg)
# handle wrap: circular mean
print("bootstrap phi0 deg: median %.1f, CI [%.1f, %.1f] (raw, may wrap)" % (md, lo, hi))
print("bootstrap phi0 mod 360 histogram coarse:", np.histogram(btdeg%360, bins=8, range=(0,360))[0])
