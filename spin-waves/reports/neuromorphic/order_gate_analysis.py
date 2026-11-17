"""Analyse order_gate_raw.jsonl: can accumulated modal power classify pump ORDER (AB vs BA)?
Run: ~/.venvs/nvim-python/bin/python order_gate_analysis.py
"""
import json, collections, sys
import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.metrics import roc_auc_score

RAW = sys.argv[1] if len(sys.argv) > 1 else 'order_gate_raw.jsonl'
rng = np.random.default_rng(0)
rows = [json.loads(l) for l in open(RAW)]
by = collections.defaultdict(list)
for r in rows: by[r['cond']].append(r)

FLOOR = 1e-9            # power floor (|a|^2) ~ SA noise floor proxy
MEAS_SIGMA = 0.09       # lognormal sigma on power: ~3.6-4 % CV on amplitude-like quantity -> ~8-9 % on power (conservative)

def feats(p, mode):
    p = np.asarray(p, float)
    p = p * np.exp(rng.normal(0, MEAS_SIGMA, p.shape)) + FLOOR * np.abs(rng.normal(1, 0.3, p.shape))
    if mode == 'all29':     return np.log10(p)
    if mode == 'two_modes': return np.log10(p[[13, 16]])
    if mode == 'm13_only':  return np.log10(p[[13]])       # what a single 5-MHz SA window could see
    if mode == 'energy':    return np.log10(np.array([p.sum()]))   # order-blind control
    raise ValueError

def evaluate(rs, mode, nperm=200):
    pairs = sorted({r['pair'] for r in rs})
    X, y, g = [], [], []
    for r in rs:
        X.append(feats(r['p'], mode)); y.append(1 if r['order'] == 'AB' else 0); g.append(r['pair'])
    X, y, g = np.array(X), np.array(y), np.array(g)
    k = min(6, len(pairs))
    clf = LinearDiscriminantAnalysis(solver='lsqr', shrinkage='auto')
    pr = cross_val_predict(clf, X, y, groups=g, cv=GroupKFold(n_splits=k), method='predict_proba')[:, 1]
    acc = float(((pr > .5) == y).mean()); auc = float(roc_auc_score(y, pr))
    # label-permutation null (permute labels within pair)
    null = []
    for _ in range(nperm):
        y2 = y.copy()
        for pid in pairs:
            idx = np.where(g == pid)[0]
            if rng.random() < .5: y2[idx] = y2[idx][::-1]
        p2 = cross_val_predict(clf, X, y2, groups=g, cv=GroupKFold(n_splits=k), method='predict_proba')[:, 1]
        null.append(((p2 > .5) == y2).mean())
    pval = float((1 + np.sum(np.array(null) >= acc)) / (1 + nperm))
    return acc, auc, pval

out = {}
print(f"{'condition':28s} {'n_pairs':>7s} | " + " | ".join(f'{m:>22s}' for m in ['all29', 'two_modes', 'm13_only', 'energy']))
for cond, rs in by.items():
    n = len(rs) // 2
    # raw modal contrast (no readout noise): log10 p13/p16, AB vs BA
    lr = {o: np.array([np.log10((r['p'][13] + 1e-30) / (r['p'][16] + 1e-30)) for r in rs if r['order'] == o]) for o in ('AB', 'BA')}
    res = {}
    cells = []
    for m in ['all29', 'two_modes', 'm13_only', 'energy']:
        acc, auc, pv = evaluate(rs, m, nperm=100)
        res[m] = dict(acc=acc, auc=auc, p=pv); cells.append(f'acc {acc:.2f} auc {auc:.2f} p={pv:.3f}')
    res['mean_log10_p13_over_p16'] = {o: float(v.mean()) for o, v in lr.items()}
    out[cond] = dict(n_pairs=n, **res)
    print(f'{cond:28s} {n:7d} | ' + ' | '.join(f'{c:>22s}' for c in cells))
json.dump(out, open('order_gate_summary.json', 'w'), indent=1)
