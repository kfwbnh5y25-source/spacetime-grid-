import numpy as np
P0=np.load("P_on_peak.npy"); P4=np.load("P_off_peak.npy")
A2=np.arange(200,276); F2=np.arange(4075,4151)
PHI_ON, PHI_OFF, PHI0_REF = 1.310, 4.550, np.deg2rad(-155.0)

def pooled_fit(mask):
    i,j = np.where(mask)
    m = (P0[i,j]+P4[i,j])/2.0
    u_on = P0[i,j]/m; u_off = P4[i,j]/m
    ths = np.linspace(0, 2*np.pi, 4001)
    c1 = np.cos(PHI_ON-ths)  # (T,)
    c2 = np.cos(PHI_OFF-ths)
    n = len(u_on)
    num = c1*(u_on-1).sum() + c2*(u_off-1).sum()
    den = n*(c1**2 + c2**2)
    g = num/den                                   # (T,)
    pred1 = 1+g[None,:]*c1[None,:]; pred2 = 1+g[None,:]*c2[None,:]
    ss = ((u_on[:,None]-pred1)**2).sum(0) + ((u_off[:,None]-pred2)**2).sum(0)
    tot = ((u_on-u_on.mean())**2).sum() + ((u_off-u_off.mean())**2).sum()
    r2 = 1-ss/tot
    k = int(np.argmax(r2))
    return dict(n=n, r2=r2[k], phi0=ths[k], g=g[k], u_on = u_on, u_off=u_off,
                i=i, j=j, r2_curve=r2, phi0_all=ths, g_all=g, c1=c1, c2=c2)

mask_loose = (P0>5e-7)|(P4>5e-7)
mask_strong = (P0>1e-6)|(P4>1e-6)
for nm,mask in [("loose(5e-7)",mask_loose),("strong(1e-6)",mask_strong)]:
    f = pooled_fit(mask)
    print(f"{nm}: n={f['n']} R2={f['r2']:.4f} phi0={np.degrees(f['phi0']):.1f} deg g={f['g']:.4f}")
    # fixed-phi0 predictive test
    t = PHI0_REF
    c1, c2 = np.cos(PHI_ON-t), np.cos(PHI_OFF-t)
    i,j,f2_ = f['i'], f['j'], None
    uo, ui = f['u_off'], f['u_on']; n=f['n']
    gg = (c1*(ui-1).sum()+c2*(uo-1).sum())/(n*(c1**2+c2**2))
    ss = ((ui-1-gg*c1)**2).sum()+((uo-1-gg*c2)**2).sum()
    tot = ((ui-ui.mean())**2).sum()+((uo-uo.mean())**2).sum()
    print(f"   fixed phi0=-155deg: g={gg:.4f} R2={1-ss/tot:.4f}")
f = pooled_fit(mask_loose)
# residual structure
pred_on = 1+f['g']*f['c1'][f['phi0_all'].searchsorted(f['phi0'])]
res = np.r_[f['u_on']-pred_on, f['u_off']-(2-pred_on)]
a2g = A2[f['i']]; f2g = F2[f['j']]
aa = np.r_[a2g,a2g]; ff=np.r_[f2g,f2g]
for nm,v in [("a2",aa),("f2",ff)]:
    b = np.quantile(v, np.linspace(0,1,7)); mm=[res[(v>=b[k])&(v<=b[k+1])].mean() for k in range(6)]
    print("resid mean vs %s sextiles:" % nm, np.round(mm,4))
print("residual rms %.4f" % res.std())
# -- fast bootstrap phi0 for strong mask
fs = pooled_fit(mask_strong)
rng = np.random.RandomState(42)
uo, ui = fs['u_off'], fs['u_on']; n=fs['n']
T = 1001; ths = np.linspace(0,2*np.pi,T)
C1 = np.cos(PHI_ON-ths)[None,:]; C2 = np.cos(PHI_OFF-ths)[None,:]
U1 = (ui-1)[:,None]; U2 = (uo-1)[:,None]
g_all = (C1*U1 + C2*U2).sum(0)/ (n*(C1[0]**2+C2[0]**2))
SS = ((U1-1-g_all[None,:]*C1)**2).sum(0)+((U2-1-g_all[None,:]*C2)**2).sum(0)
bt=[]
for _ in range(500):
    idx = rng.randint(0,n,n)
    gg = (C1*U1[idx] + C2*U2[idx]).sum(0)/(n*(C1[0]**2+C2[0]**2))
    ss = ((U1[idx]-1-gg[None,:]*C1)**2).sum(0)+((U2[idx]-1-gg[None,:]*C2)**2).sum(0)
    bt.append(ths[ss.argmin()])
bt=np.degrees(np.array(bt))
print("strong-mask bootstrap phi0 deg: median %.1f, 68%% CI [%.1f, %.1f]" % (np.median(bt), *np.percentile(bt,[16,84])))
