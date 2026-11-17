"""
runner.py - CENTRAL RUNNER for the overnight phase-gating campaign.

Order (checked against the 8 h budget; ETA per step printed live):
  step 1  meas5_edges.py            ~10 min   edge sweep 0.5-10 ns (incl. 1 ns drift control)
  step 2  meas1 --config 0          ~55 min   long-gap memory, Df=50 MHz, 60 phases
  step 3  meas1 --config 1          ~55 min   long-gap memory, Df=58 MHz, 60 phases
  step 4  meas2_threshold_vs_phase  ~20 min   h2,th(Psi), a2 2 mV steps x phases
  step 5  meas3_depth_repeat        ~25 min   money-OP depth, 60 phases x repeats
  step 6  meas1 --config 2          ~15 min   long-gap, Df=64 MHz, 16 phases
  step 7  meas1 --config 3          ~15 min   long-gap, Df=70 MHz, 16 phases
  step 8  meas4_duty                ~5  min   duty/comb control at probe OP

NOT included: meas0_sanity.py (run manually first; sanity = GO required).

Each step is skipped if its output h5 already exists with >0 datasets
(resume-safe: rerun runner after an interruption). Use --force to redo.

Usage:
  python runner.py                  # full campaign
  python runner.py --list           # show steps, ETAs, output files
  python runner.py --only 2,4      # run only steps 2 and 4
  python runner.py --force-step 1  # redo step 1 even if output exists
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable

STEPS = [
    (1, "edges: tukey edge sweep 0.5-10 ns (incl. 1 ns control)",
     ["meas5_edges.py"], HERE / "meas5_edges.h5", "10 min"),
    (2, "longgap Df=50 MHz, 60 phases, delays -23..+200 ns",
     ["meas1_longgap_memory.py", "--config", "0", "--out", "meas1_longgap_memory_c0.h5"],
     HERE / "meas1_longgap_memory_c0.h5", "55 min"),
    (3, "longgap Df=58 MHz, 60 phases",
     ["meas1_longgap_memory.py", "--config", "1", "--out", "meas1_longgap_memory_c1.h5"],
     HERE / "meas1_longgap_memory_c1.h5", "55 min"),
    (4, "threshold vs phase (a2 2 mV steps)",
     ["meas2_threshold_vs_phase.py"], HERE / "meas2_threshold_vs_phase.h5", "20 min"),
    (5, "depth repeats at money OP (60 phases x 5 trials)",
     ["meas3_depth_repeat.py"], HERE / "meas3_depth_repeat.h5", "25 min"),
    (6, "longgap Df=64 MHz, 16 phases",
     ["meas1_longgap_memory.py", "--config", "2", "--trials", "1",
      "--out", "meas1_longgap_memory_c2.h5"],
     HERE / "meas1_longgap_memory_c2.h5", "15 min"),
    (7, "longgap Df=70 MHz, 16 phases",
     ["meas1_longgap_memory.py", "--config", "3", "--trials", "1",
      "--out", "meas1_longgap_memory_c3.h5"],
     HERE / "meas1_longgap_memory_c3.h5", "15 min"),
    (8, "duty/comb control at probe OP",
     ["meas4_duty.py"], HERE / "meas4_duty.h5", "5 min"),
]


def h5_has_data(p: Path) -> bool:
    if not p.exists():
        return False
    try:
        import h5py
        with h5py.File(p, "r") as f:
            return len(f.keys()) > 0
    except Exception:
        return False


def banner(msg):
    line = "-" * 66
    print(f"\n{line}\n{msg}\n{line}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--only", type=str, default=None, help="comma-separated step ids")
    ap.add_argument("--force-step", type=str, default=None,
                    help="comma-separated step ids to redo even if output exists")
    ap.add_argument("--dry-run-outer", action="store_true", help="print commands only")
    args = ap.parse_args()

    if args.list:
        print(f"{'step':>4}  {'ETA':>7}  {'output':34s}  description")
        for n, desc, cmd, out, eta in STEPS:
            state = "DONE (exists)" if h5_has_data(out) else "todo"
            print(f"{n:>4}  {eta:>7}  {out.name:34s}  {desc}  [{state}]")
        return

    only = set(int(x) for x in args.only.split(",")) if args.only else None
    force = set(int(x) for x in args.force_step.split(",")) if args.force_step else set()
    todo = [(n, desc, cmd, out, eta) for n, desc, cmd, out, eta in STEPS
            if only is None or n in only]

    t_all = time.time()
    log = []
    for n, desc, cmd, out, eta in todo:
        banner(f"STEP {n}: {desc}   (ETA {eta}, output {out.name})")
        if h5_has_data(out) and n not in force:
            print(f"  skip: {out.name} already exists and is non-empty "
                  f"(use --force-step {n} to redo)")
            continue
        full = [PY] + [str(HERE / c) if c.endswith(".py") else c for c in cmd]
        print("  cmd:", " ".join(full[1:]), flush=True)
        if args.dry_run_outer:
            continue
        t0 = time.time()
        try:
            r = subprocess.run(full, cwd=HERE)
            status = f"[{r.returncode}]"
        except KeyboardInterrupt:
            print(f"\n  INTERRUPTED step {n} after {time.time() - t0:.0f} s")
            raise
        line = f"step {n} {status} in {(time.time() - t0) / 60:.1f} min -> {out.name}"
        print(line, flush=True)
        log.append(line)

    banner(f"RUNNER FINISHED in {(time.time() - t_all) / 60:.1f} min")
    for l in log:
        print("  " + l)
    missing = [out.name for _, _, _, out, _ in todo if not h5_has_data(out)]
    if missing:
        print("  MISSING OUTPUTS:", ", ".join(missing))
    else:
        print("  all requested outputs present")


if __name__ == "__main__":
    main()
