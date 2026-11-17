"""
meas5_edges.py - EDGE-TUKEY SWEEP (boss request): effect of pulse edge ramp
on phase gating at the sanity operating point, ~7 min.

MOTIVATION: with 32 ns tukey edges the idler line exists (A0 ~ 170 uW) but
phase gating VANISHES (A1/A0 = 0.00, live 1117 run); with 1 ns edges gating
is back (A1/A0 = 0.34). This scan maps A1/A0 vs edge_ns at fixed June
geometry (tone1 4075/327 mV 300 ns, t_gap=0, tone2 4149/588 mV 400 ns,
2 us loop). Edge grid 0.5 -> 10 ns + a 1 ns control repeat at the end
(drift check).

PREDICTION (from the two live points): A1/A0 drops steeply with edge_ns;
roughly 0.35 at 0.5-1 ns, clearly below by 4 ns, near zero by 10 ns,
matching the 32 ns data point (0.00).
OUTCOME A (steep monotone): edge ramp = the phase-smearing knob; record the
falloff, cite in Methods (JuneScripts used 64-sample tukey edge).
OUTCOME B (A1 flat ~0.3 everywhere): edge NOT the knob at small edges -> the
32 ns anomaly needs another cause (revisit before campaign conclusions).

Same analysis + quals as meas0 (line_power, c12). Live per-edge verdict.

Usage: python meas5_edges.py [--dry-run] [--out PATH]
"""
import argparse
import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "")
import meas0_sanity as ms
from hw_core_spectro import SAConfig, SpectroRig

EDGES_NS = [0.5, 1.0, 2.0, 4.0, 7.0, 10.0]
N_PHI = len(ms.PHASES_T)
EXPECT = {0.5: 0.35, 1.0: 0.34, 2.0: 0.2, 4.0: 0.1, 7.0: 0.03, 10.0: 0.02}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--out", default="meas5_edges.h5")
    ap.add_argument("--peak-mv", type=float, default=916.0)
    args = ap.parse_args()

    rig = SpectroRig(args.out, dry_run=args.dry_run,
                     sa_config=SAConfig(center_hz=ms.CENTER_HZ, span_hz=ms.SPAN_HZ,
                                        points=1001, rbw_hz=1e6, avg_count=1,
                                        sweep_time=0.5, ref_level_dbm=-53.8,
                                        scale_div=0.5))
    t0 = time.time()
    plan = [(e, float(ph)) for e in EDGES_NS for ph in ms.PHASES_T]
    plan += [(1.0, float(ph)) for ph in ms.PHASES_T]     # end control at 1 ns
    print(f"meas5: {len(plan)} measures x (BL+SIG), ETA ~ {len(plan)*1.2/60:.0f} min")
    rows_by_edge = {}
    try:
        for i, (edge, ph) in enumerate(plan):
            wf, spec = ms.wf_june(ph, 0.0, edge_ns=edge)
            if wf is None:
                continue
            tag = f"EDGE{edge:.1f}"
            key = f"ED_{tag}_p{ph:.4f}"
            if edge == 1.0 and i >= len(EDGES_NS) * N_PHI:
                key = f"ED_END1.0_p{ph:.4f}"
            rig.register_waveform(key, wf, spec)
            sig, bl, fr = rig.measure(key, wf, meta={"edge_ns": edge, "phi": ph})
            if args.dry_run:
                sig = sig + ms.model_line(fr, ph, 0.0)
            p, n = ms.line_power(sig, bl, fr)
            rows_by_edge.setdefault(edge, []).append((ph, p))
            ms.live(f"{tag} #{i}", ph, p, n, [None])
        print("\n-- edge sweep summary --")
        summary = []
        for edge in EDGES_NS + ([1.0] if len(EDGES_NS) * N_PHI < len(plan) else []):
            rows = rows_by_edge.get(float(edge), [])
            if len(rows) < 8:
                continue
            A1, A0, phi0 = ms.c12(sorted(rows))
            m = A1 / max(A0, 1e-30)
            verdict = "OK" if m >= 0.2 else ("LOW" if m >= 0.08 else "NONE")
            print(f"  edge={edge:5.1f} ns: A0={A0*1e6:7.1f} uW  A1={A1*1e6:7.1f} uW  "
                  f"A1/A0={m:.3f} [{verdict}]  phi0={phi0:+.2f} rad  (expect ~{EXPECT.get(edge, '?')})")
            summary.append((edge, A1, A0, m))
        if summary:
            fig, axs = plt.subplots(1, 2, figsize=(11, 4))
            e = [s[0] for s in summary]
            axs[0].plot(e, [s[3] for s in summary], "o-")
            axs[0].axhline(0.42, color="g", ls="--", label="June anchor 0.42")
            axs[0].axhline(0.0, color="k", ls=":", alpha=0.4)
            axs[0].set_xlabel("edge_ns"); axs[0].set_ylabel("A1/A0")
            axs[0].set_title("phase gating vs tukey edge"); axs[0].legend()
            axs[1].plot(e, [s[1] * 1e6 for s in summary], "o-")
            axs[1].set_xlabel("edge_ns"); axs[1].set_ylabel("A1 (uW)")
            axs[1].set_title("modulation amplitude")
            plt.tight_layout()
            outfig = __import__("os").path.join(
                "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/tmp",
                "meas5_edges.png")
            plt.savefig(outfig, dpi=130)
            print(f"  figure: {outfig}")
    finally:
        print(f"\nmeas5 done in {time.time() - t0:.0f} s, out={args.out}")


if __name__ == "__main__":
    main()
