import numpy as np
rng = np.random.RandomState(0)
y_on = rng.randn(1001)*1e-7 + 1e-7
y_ref = rng.randn(1001)*1e-7 + 5e-8
y_off = rng.randn(1001)*1e-7 + 2e-7
def to_dbm(pw):
    pw_pos = np.where(pw > 0, pw, np.nan)
    return 10.0 * np.log10(pw_pos / 1e-3)
t1 = to_dbm(y_on); t2 = to_dbm(y_ref); t3 = to_dbm(y_off)
print("t1 finite:", np.isfinite(t1).sum(), "t2 finite:", np.isfinite(t2).sum(), "t3 finite:", np.isfinite(t3).sum())
with np.errstate(all="ignore"):
    ymin = np.nanmin([t1.min(), t2.min()]) - 3
    ymax = np.nanmax([t1.max(), t3.max()]) + 4
print("ymin", ymin, "ymax", ymax)
ymin2 = np.nanmin([t1.min(), t2.min()]) - 3
print("ymin2", ymin2)
print("numpy", np.__version__)
