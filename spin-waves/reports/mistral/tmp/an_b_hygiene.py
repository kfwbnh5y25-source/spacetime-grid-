import numpy as np, pickle
from datetime import datetime
OUT='reports/mistral/tmp/'
blocks=pickle.load(open(OUT+'B_blocks.pkl','rb'))
ts=[]
for b in blocks:
    for t in b['ts']: ts.append(t)
ts=sorted(datetime.fromisoformat(t) for t in ts)
span=(ts[-1]-ts[0]).total_seconds()
print('keys %d, acquisition span %.1f min, median dt %.3f s'%(len(ts),span/60, np.median(np.diff(np.array([t.timestamp() for t in ts])))))
days=sorted({t.date().isoformat() for t in ts})
print('days:',days)
perday={}
for d in days: perday[d]=sum(1 for t in ts if t.date().isoformat()==d)
print(perday)
# corruption threshold sensitivity
print('\ncorruption screen re-run at softer thresholds:')
tot=len(blocks)
for thr in [4,5,6,8,10]:
    nc=0
    for b in blocks:
        n=b['n']; O=b['off']; o=np.argsort(b['ts']); Oq=O[o]
        best=None
        for cut in range(3,n-2):
            s1=Oq[:cut]-Oq[:cut].mean(); s2=Oq[cut:]-Oq[cut:].mean()
            v=(np.sum(s1**2)/max(cut-1,1)+np.sum(s2**2)/max(n-cut-1,1))
            dlev=Oq[cut:].mean()-Oq[:cut].mean()
            se=np.sqrt(v/cut+v/(n-cut))
            t=dlev/se if se>0 else 0
            if best is None or abs(t)>abs(best): best=abs(t)
        if best is not None and best>thr: nc+=1
    print('  thr=%d sigma: %d corrupted of %d (%.1f%%)'%(thr,nc,tot,100*nc/tot))
# which blocks corrupted at 6
print('\ncorrupted blocks (>6 sigma):')
for b in blocks:
    if b['corrupt']: print('  a1=%.0f f1=%d n=%d'%(b['a1'],b['f1'],b['n']))
