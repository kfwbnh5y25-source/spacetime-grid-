import numpy as np
P0=np.load("P_on_peak.npy"); P4=np.load("P_off_peak.npy")
PHI_ON, PHI_OFF, PHI0_REF = 1.310, 4.550, np.deg2rad(-155.0)
mask = (P0>1e-6)|(P4>1e-6)
i,j = np.where(mask)
m = (P0[i,j]+P4[i,j])/2.0
u1 = P0[i,j]/m - 1.0   # ON residual target
u2 = P4[i,j]/m - 1.0
n = len(u1)
print("n", n, "mean u1 %.4f mean u2 %.4f, sd %.3f %.3f" % (u1.mean(), u2.mean(), u1.std(), u2.std()))
ths = np.linspace(0,2*np.pi,3601)
c1 = np.cos(PHI_ON-ths); c2 = np.cos(PHI_OFF-ths)
g = (c1*u1.sum() + c2*u2.sum())/(n*(c1**2+c2**2))
ss = ((u1[:,None]-g[None,:]*c1[None,:])**2).sum(0) + ((u2[:,None]-g[None,:]*c2[None,:])**2).sum(0)
T=(u1**2).sum()+(u2**2).sum()
print("SS min %.2f at %.2f deg (g=%.4f); SS at -155deg (205deg): %.2f" % (ss.min(), np.degrees(ths[ss.argmin()]), g[ss.argmin()], ss[np.argmin(np.abs(np.degrees(ths)-205))]))
# top-5 minima
order = np.argsort(ss)
print("argmin deg top5:", np.round(np.degrees(ths[order[:5]]),2), "SS frac:", np.round(ss[order[:5]]/T,6))
# SS curve profile every 10 deg
print("SS(φ0)/T every 10 deg:", np.round(ss[::10]/T,5))
