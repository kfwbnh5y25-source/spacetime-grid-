src = open("../pipeline/fig3_overlap.py").read()
marker = "ymin = np.nanmin"
dbg = ('print("DBG", ROW, ROW_F2, "y_on n-nan", int(np.isnan(y_on).sum()),'
       ' "y_ref n-nan", int(np.isnan(y_ref).sum()), "td_on allnan", bool(np.isnan(to_dbm(y_on)).all()),'
       ' "td_ref allnan", bool(np.isnan(to_dbm(y_ref)).all()), "td_off allnan", bool(np.isnan(to_dbm(y_off)).all()))\n')
src = src.replace(marker, dbg + marker, 1)
open("dbg3.py","w").write(src)
