import sys, json
sys.path.insert(0, "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/pipeline")
from provenance import Provenance

p = Provenance("smoke_test")
data = p.load_traces("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/A_old_params/short/other_dir/data_broad2.h5",
                     select=lambda k: k.split("_")[1] == "266", pass_agg="mean")
print("loaded", len(data), "keys; sample:", list(data)[:2])
p.record_source("/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating/detune_fine031126.h5", note="fig4 source")
p.step("smoke", params={"x": 1}, result={"n": len(data)})
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(); ax.plot([0, 1], [0, 1])
mp = p.finish(fig, ["panel_a"], formats=("png",))
print("manifest summary:", {k: (len(v) if isinstance(v, list) else str(v)[:60]) for k, v in p.manifest.items()})
