import numpy as np
P0=np.load("P_on_peak.npy"); P4=np.load("P_off_peak.npy")
A2=range(200,276); F2=range(4075,4151)
C=10*np.log10(P4/P0)
# where is line strong?
ii,jj = np.unravel_index(np.argmax(P0),P0.shape)
print("max P_on: a2=%d f2=%d P=%.3g ; P_off same cell %.3g" % (200+ii,4075+jj,P0[ii,jj],P4[ii,jj]))
big = np.argwhere(P0>1e-6)
print("cells P_on>1e-6:", len(big), "a2 range", big[:,0].min()+200 if len(big) else '-', big[:,0].max()+200 if len(big) else '-',
      "f2 range", (big[:,1].min()+4075 if len(big) else '-'), (big[:,1].max()+4075 if len(big) else '-'))
# count cells with line > 1e-7 (above ~noise)
print("cells P_on>1e-7:", (P0>1e-7).sum(), " P_off>1e-7:", (P4>1e-7).sum())
# distribution of C in line region
mask = P0>5e-7
print("C in line-region: n=%d mean %.2f min %.2f max %.2f" % (mask.sum(), C[mask].mean(), C[mask].min(), C[mask].max()))
# row 275 and 272 summaries
for r in (272,275,204):
    row = C[r-200]
    print("row a2=%d: C mean %.2f, min %.2f (f2=%d), max %.2f (f2=%d)" % (r, row.mean(), row.min(), F2[row.argmin()] if False else row.argmin()+4075, row.max(), row.argmax()+4075))
# quadrant means
print("C>0 frac:", (C>0.2).mean(), " C<-0.2 frac:", (C<-0.2).mean())
