"""
provenance.py - trace-to-figure reproducibility layer for Paper A figures.

Every figure script imports this and:
    prov = Provenance("fig3_overlap")
    data = prov.load_traces(h5_path, keys_or_filter, loader=my_loader, pass_agg="mean")
    ... analysis ...
    prov.step("phase_contrast", params={...}, note="...")
    prov.finish(fig, ["panel_a", "panel_b"])   # saves png+pdf+manifest json

The manifest (figures/manifests/<fig>.json) then contains the full path:
    raw h5 (path, sha256, mtime, script_content sha) -> keys used -> every
    processing step + params -> figure file sha256. Rerunning the figure
    script reproduces the figure from the manifest alone.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import os
import platform
import sys
from datetime import datetime, timezone

import h5py
import numpy as np

REPORT_DIRS = [
    "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/phase gating",
    "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral",
]
FIG_ROOT = "/Users/davidzehner/Work/NMM_TOOL/spacetime-grid-/spin-waves/reports/mistral/figures"


def sha256_file(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


class Provenance:
    def __init__(self, fig_name: str):
        self.fig_name = fig_name
        self.manifest = {
            "figure": fig_name,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "env": {"python": sys.version.split()[0],
                    "numpy": np.__version__,
                    "platform": platform.platform(),
                    "modules": {m: getattr(importlib.import_module(m), "__version__", "?")
                                for m in ["h5py", "scipy", "matplotlib"] if _importable(m)}},
            "sources": [],
            "steps": [],
            "outputs": [],
        }
        os.makedirs(FIG_ROOT, exist_ok=True)
        os.makedirs(os.path.join(FIG_ROOT, "manifests"), exist_ok=True)

    # ---- sources -----------------------------------------------------------

    def record_source(self, h5_path: str, note: str = "") -> dict:
        rec = {"path": os.path.abspath(h5_path),
               "sha256": sha256_file(h5_path),
               "size_bytes": os.path.getsize(h5_path),
               "mtime": datetime.fromtimestamp(os.path.getmtime(h5_path)).isoformat(),
               "note": note}
        try:
            with h5py.File(h5_path, "r") as f:
                sc = f.attrs.get("script_content")
                if sc is not None:
                    b = sc if isinstance(sc, bytes) else str(sc).encode()
                    rec["script_content_sha256"] = hashlib.sha256(b).hexdigest()
        except Exception as e:
            rec["script_content_error"] = str(e)
        self.manifest["sources"].append(rec)
        return rec

    # ---- loading with automatic key bookkeeping ----------------------------

    def load_traces(self, h5_path: str, select=None, loader=None, pass_agg: str = "mean",
                    freq_from_attr0: bool = True):
        """Generic loader for the kNET HDF5Buffer layout.
        select: callable(key) -> bool  or list of keys.
        loader: optional callable(key, trace_array) -> anything (default: pass
        the raw (passes, ch, n) array). pass_agg: how to aggregate the
        pass/append axis: 'mean' | 'last' | 'all'.
        Records the selected keys (sorted) into the manifest.
        Returns {key: loaded}."""
        with h5py.File(h5_path, "r") as f:
            all_keys = sorted(f.keys())
            if callable(select):
                keys = [k for k in all_keys if select(k)]
            elif select is not None:
                keys = [k for k in all_keys if k in set(select)]
            else:
                keys = all_keys
            out = {}
            for k in keys:
                tr = f[k]["trace"][:]
                if loader is not None:
                    out[k] = loader(k, tr)
                else:
                    arr = tr[:, 1, :].astype(np.float64) if tr.ndim == 3 and tr.shape[1] >= 2 else tr
                    if pass_agg == "mean":
                        out[k] = arr.mean(axis=0)
                    elif pass_agg == "last":
                        out[k] = arr[-1]
                    else:
                        out[k] = arr
        self.manifest["sources"].append(
            {"path": os.path.abspath(h5_path), "keys_used_n": len(keys),
             "keys_used": keys if len(keys) <= 400 else keys[:200] + ["...truncated..."] + keys[-200:],
             "pass_agg": pass_agg})
        return out

    # ---- processing steps --------------------------------------------------

    def step(self, name: str, params: dict | None = None, note: str = "", result: dict | None = None):
        entry = {"step": name, "time_utc": datetime.now(timezone.utc).isoformat(),
                 "params": _jsonable(params or {}), "note": note,
                 "result": _jsonable(result or {})}
        self.manifest["steps"].append(entry)
        return entry

    def array(self, name: str, arr) -> str:
        """Store intermediate arrays next to the figure for reproducibility.
        Returns relative npy filename recorded in the manifest."""
        fn = f"{self.fig_name}__{name}.npy"
        np.save(os.path.join(FIG_ROOT, fn), np.asarray(arr))
        self.manifest["steps"].append({"step": "save_array", "name": fn,
                                       "shape": list(np.asarray(arr).shape),
                                       "sha256": hashlib.sha256(np.asarray(arr).tobytes()).hexdigest()})
        return fn

    # ---- output ------------------------------------------------------------

    def finish(self, fig, panel_labels: list | None = None, formats=("png", "pdf")):
        import matplotlib.pyplot as plt
        outs = []
        for fmt in formats:
            fn = f"{self.fig_name}.{fmt}"
            fig.savefig(os.path.join(FIG_ROOT, fn), dpi=300, bbox_inches="tight")
            outs.append({"file": fn, "sha256": sha256_file(os.path.join(FIG_ROOT, fn))})
        self.manifest["outputs"] = outs
        self.manifest["panels"] = panel_labels or []
        mpath = os.path.join(FIG_ROOT, "manifests", f"{self.fig_name}.json")
        with open(mpath, "w") as f:
            json.dump(self.manifest, f, indent=2)
        try:
            plt.close(fig)
        except Exception:
            pass
        print(f"[provenance] {self.fig_name}: {len(outs)} outputs, manifest -> {mpath}")
        return mpath


def _importable(m):
    try:
        importlib.import_module(m)
        return True
    except Exception:
        return False


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (str, int, float, bool)) or o is None:
        return o
    return str(o)
