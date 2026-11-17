"""
Generate Publication-Grade Standard Neuromorphic RC Benchmark Summary Figure
from all verified experimental datasets + Digital ESN baselines.
Outputs: spin-waves/reports/neuromorphic/figures/fig_sota_rc_benchmarks.{png,pdf}
"""

import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
NM_DIR = ROOT / "neuromorphic"
FIG_DIR = Path(__file__).resolve().parent / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main():
    params = json.loads((ROOT / "reports" / "idler_math" / "idler_extracted_params.json").read_text())
    me = params["memory_e3fix"]
    v2 = json.loads((NM_DIR / "multisymbol_v2_metrics.json").read_text())
    dig = json.loads((NM_DIR / "digital_baseline_comparison.json").read_text())

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0), dpi=220)

    # (a) Fading-Memory STM Decay Curve (E3FIX)
    ax = axes[0, 0]
    ds = [0, 10, 20, 40, 80, 150, 300]
    acc350 = [100.0, 100.0, 100.0, 43.3, 66.7, 53.3, 40.0]
    acc380 = [100.0, 100.0, 100.0, 63.3, 70.0, 56.7, 46.7]
    ax2 = ax.twinx()
    l1 = ax.plot(ds, acc350, "o-", color="#1f77b4", lw=2.0, ms=6, label="Acc 350 mV (%)")
    l2 = ax.plot(ds, acc380, "s-", color="#2ca02c", lw=2.0, ms=6, label="Acc 380 mV (%)")
    l3 = ax2.semilogy(me["380"]["delays_ns"], me["380"]["fisher_j"], "D--", color="#d62728", lw=1.8, ms=5,
                      label=f"Fisher $J$ 380 mV ($\\tau_J={me['380']['tau_J_ns']:.1f}$ ns)")
    ax.axhline(50.0, color="gray", ls=":", lw=1.2)
    ax.set_xlabel("Inter-Symbol Delay $d$ (ns)")
    ax.set_ylabel("5-Fold CV Accuracy (%)")
    ax2.set_ylabel("Fisher Separability $J$", color="#d62728")
    ax.set_ylim(30, 108)
    ax.set_title("(a) Short-Term Fading Memory (`E3FIX`: $\\tau_J \\approx 4.3\\text{–}9.4\\,\\mathrm{ns}$)")
    lines = l1 + l2 + l3
    ax.legend(lines, [l.get_label() for l in lines], loc="upper right", fontsize=8.0)
    ax.grid(True, alpha=0.3)

    # (b) Standard Physical RC Benchmarks vs Controls
    ax = axes[0, 1]
    names = ["2-Sym\nOrder", "S3\n(6-cls)", "Matched\n2-hot", "Matched\n1-hot", "BIT3\n(8-way)"]
    phys = [100.0, v2["S3"]["acc"]*100, v2["BIT3_2hot3"]["acc"]*100, v2["BIT3_1hot3"]["acc"]*100, v2["BIT3_8way"]["acc"]*100]
    ctrl = [41.9, v2["S3_03x"]["acc"]*100, v2["BIT3_03x_2hot3"]["acc"]*100, v2["BIT3_03x_1hot3"]["acc"]*100, v2["BIT3_03x_8way"]["acc"]*100]
    chances = [50.0, 16.7, 33.3, 33.3, 12.5]
    x = np.arange(len(names))
    ax.bar(x - 0.18, phys, width=0.35, color="#2ca02c", edgecolor="k", label="Supra-Threshold Reservoir")
    ax.bar(x + 0.18, ctrl, width=0.35, color="#7f7f7f", edgecolor="k", label="0.3x Sub-Threshold Ctrl")
    for xi, ch in zip(x, chances):
        ax.plot([xi - 0.38, xi + 0.38], [ch, ch], "r--", lw=1.3, label="Static / Chance Cap" if xi == 0 else None)
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=9)
    ax.set_ylim(0, 112)
    ax.set_ylabel("5-Fold CV Accuracy (%)")
    ax.set_title("(b) Temporal Sequence & Order Benchmarks")
    ax.legend(loc="upper right", fontsize=8.5)
    ax.grid(True, axis="y", alpha=0.3)

    # (c) Coherent Phase-Gated Logic (PHASE_GATE)
    ax = axes[1, 0]
    pg = v2["phase_gate"]
    phis = sorted([int(k) for k in pg.keys()])
    pg_acc = [pg[str(p)]["acc"] * 100 for p in phis]
    pg_j = [pg[str(p)]["fisher_j"] for p in phis]
    ax2c = ax.twinx()
    l1c = ax.plot(phis, pg_acc, "o-", color="#9467bd", lw=2.2, ms=7, label="Order Accuracy (%)")
    l2c = ax2c.plot(phis, pg_j, "D--", color="#ff7f0e", lw=2.0, ms=6, label="Fisher $J$")
    ax.set_xticks(phis); ax.set_xticklabels([f"{p}°" for p in phis])
    ax.set_ylim(75, 105)
    ax.set_xlabel("Relative Carrier Phase $\\Delta\\phi$")
    ax.set_ylabel("Order Accuracy (%)", color="#9467bd")
    ax2c.set_ylabel("Fisher Separability $J$", color="#ff7f0e")
    ax.set_title("(c) Coherent Phase-Gated Interference (`AB` vs `BA`)")
    linesc = l1c + l2c
    ax.legend(linesc, [l.get_label() for l in linesc], loc="upper left", fontsize=8.5)
    ax.grid(True, alpha=0.3)

    # (d) Physical Reservoir vs Digital Baselines (Equal Rank r <= 20)
    ax = axes[1, 1]
    models = ["Order-Blind\nStatic", "Physical\nReservoir", "LTI+NL\nTwin", "Digital ESN\n(29-node)"]
    e2_vals = [dig["E2_ovl"]["Order_Blind_Static"]*100, dig["E2_ovl"]["Physical_Reservoir"]*100,
               dig["E2_ovl"]["LTI_StaticNL_Twin"]*100, dig["E2_ovl"]["Digital_ESN_29node"]*100]
    s3_vals = [dig["S3_6class"]["Order_Blind_Static"]*100, dig["S3_6class"]["Physical_Reservoir_v2"]*100,
               dig["S3_6class"]["LTI_StaticNL_Twin"]*100, dig["S3_6class"]["Digital_ESN_29node"]*100]
    xm = np.arange(len(models))
    ax.bar(xm - 0.18, e2_vals, width=0.35, color="#1f77b4", edgecolor="k", label="2-Symbol Order (`E2_ovl`)")
    ax.bar(xm + 0.18, s3_vals, width=0.35, color="#ff7f0e", edgecolor="k", label="6-Class $S_3$ (`v2`)")
    for xi, (v_e2, v_s3) in enumerate(zip(e2_vals, s3_vals)):
        ax.text(xi - 0.18, v_e2 + 1.5, f"{v_e2:.0f}%", ha="center", fontsize=8)
        ax.text(xi + 0.18, v_s3 + 1.5, f"{v_s3:.1f}%", ha="center", fontsize=8)
    ax.set_xticks(xm); ax.set_xticklabels(models, fontsize=9)
    ax.set_ylim(0, 124)
    ax.set_ylabel("5-Fold CV Accuracy (%)")
    ax.set_title("(d) Physical Reservoir vs. Digital Baselines ($r \\leq 20$)")
    ax.legend(loc="upper left", ncol=2, fontsize=8.0)
    ax.grid(True, axis="y", alpha=0.3)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig_sota_rc_benchmarks.png", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_sota_rc_benchmarks.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Saved fig_sota_rc_benchmarks.{png,pdf}")


if __name__ == "__main__":
    main()
