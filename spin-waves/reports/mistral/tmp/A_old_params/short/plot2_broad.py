import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation
from scipy.fft import fft
from tqdm import tqdm

from kNET import load_hdf


# ============================================================
# Configuration
# ============================================================

FILE = "data_broad.h5"
OUTPUT_DIR = Path("interaction_plots_broad")
OUTPUT_DIR.mkdir(exist_ok=True)

TRACE_SCALE = 2e5
REFERENCE_FREQUENCY = 4069

# ============================================================
# Load measurements
# ============================================================

rows = []

with load_hdf(FILE) as hdf:
    all_keys = list(hdf.keys())

    # Retain the same type of keys as in the original script:
    # phase = 0, label = 0, reference = 4069
    keys = [key for key in all_keys if key.endswith("0_4069")]

    print("Matching phase-zero keys:", len(keys))

    for key in tqdm(keys, desc="Loading traces"):
        match = re.fullmatch(
            r"([-+\d.]+)_([-+\d.]+)_([-+\d.]+)_([-+\d.]+)_([-+\d.]+)",
            key,
        )
        if match is None:
            continue

        a1_value, f1_value, phase_value, label_value, ref_value = map(
            float, match.groups()
        )

        if a1_value < 101:
            continue
        if f1_value < 3760:
            continue
        if f1_value == 4069777:
            continue

        try:
            trace = hdf[key]["trace"]
            frequencies = trace[:, 0].mean(axis=0)
            spectrum = trace[:, 1].mean(axis=0) * TRACE_SCALE

            # Reconstruct the other keys using the original convention.
            key_f1 = f"0.0_{f1_value:.1f}_0.000_0_4069"
            key_f2 = f"{a1_value:.1f}_{f1_value:.1f}_0.000_2_4069"

            def load_trace(k):
                return hdf[k]["trace"][:, 1].mean(axis=0) * TRACE_SCALE

            try:
                f1_only = load_trace(key_f1)
            except KeyError:
                f1_only = np.zeros_like(spectrum)

            try:
                f2_only = load_trace(key_f2)
            except KeyError:
                f2_only = np.zeros_like(spectrum)

            rows.append(
                {
                    "key": key,
                    "a1": a1_value,
                    "f1": f1_value,
                    "phase": phase_value,
                    "f12": spectrum,
                    "frequencies": frequencies,
                    "f1_only": f1_only,
                    "f2_only": f2_only,
                }
            )

        except KeyError:
            continue

if not rows:
    raise ValueError("No measurements were loaded. Check the HDF5 key format.")

df = pd.DataFrame(rows)

# ============================================================
# Assemble arrays, preserving the original definitions
# ============================================================

valid_keys = df["key"].tolist()
a1 = df["a1"].to_numpy()
f1 = df["f1"].to_numpy()
p = df["phase"].to_numpy()

res = np.stack(
    [
        np.stack(
            [
                row["f1_only"],
                row["f12"],
                row["f2_only"],
            ]
        )
        for _, row in df.iterrows()
    ]
)

sa_freq = df.iloc[0]["frequencies"]

# Integrated spectra, as in the original script.
P = res.sum(axis=-1)

# Original absolute response.
P_abs = P[:, 2].copy()


P_12 = P[:, 1].copy()

normalization = (P[:, 2]).max()
P_interact = P[:, 1] - P[:, 0] - P[:, 2]
print(P_interact.min(), P_interact.max())
print(f1[np.argsort(P_interact)][:20])
print(a1[np.argsort(P_interact)][:20])
if normalization > 0:
    P_abs = 10 * np.log10(np.maximum(np.abs(P_abs) / normalization, 1e-300))
    P_abs = np.clip(P_abs, -45, 0)
else:
    P_abs[:] = np.nan

P_12 = 10 * np.log10(P_12 /normalization)
dP_rel = np.clip(P[:, 1] / np.clip(P[:, 0] + P[:, 2], 1e-9, np.inf), 1e-9, 1e4)

dP_rel = 10 * np.log10(dP_rel)
dP_rel = np.nan_to_num(dP_rel, nan=0)
dP_rel = np.clip(dP_rel, -20, 20) * ((P[:, 2] + P[:, 1] + P[:, 0]) > 1)

print(f1[np.argsort(dP_rel)][:10])
print(a1[np.argsort(dP_rel)][:10])
print(P[np.argsort(dP_rel)][:10])
print(dP_rel[np.argsort(dP_rel)][:10])

print(dP_rel.min(), dP_rel.max())

dP_rel = np.clip(dP_rel, -40, 40) * (
    (P[:, 2] + P[:, 1] + P[:, 0]) > normalization / 1000
)


# ============================================================
# Original phase-zero maps
# ============================================================

p_unique = np.sort(np.unique(p))
mask0 = np.isclose(p, p_unique[0])

fig, ax = plt.subplots(ncols=3, figsize=(15, 6), sharex=True, sharey=True)

levels2 = np.linspace(
    -max(np.nanmax(np.abs(dP_rel[mask0])), 1e-12),
    max(np.nanmax(np.abs(dP_rel[mask0])), 1e-12),
    200,
)

P_abs = np.nan_to_num(P_abs, nan=np.nanmax(P_abs[mask0]) - 45)
P_12 = np.nan_to_num(P_12, nan=np.nanmax(P_12[mask0]) - 45)
P_12 = np.clip(P_12, -45, 0)

levels3 = np.linspace(
    max(np.nanmax(P_abs) - 45, min(np.nanmin(P_abs), np.nanmin(P_12))),
    max(np.nanmax(P_abs), np.nanmax(P_12)),
    200,
)
levels4 = np.linspace(
    max(np.nanmax(P_12) - 45, np.nanmin(P_12)),
    np.nanmax(P_12),
    200,
)

cf1 = ax[1].tricontourf(
    f1[mask0],
    a1[mask0],
    dP_rel[mask0],
    cmap="seismic",
    levels=levels2,
)
cf2 = ax[0].tricontourf(
    f1[mask0],
    a1[mask0],
    P_abs[mask0],
    cmap="viridis",
    levels=levels3,
)
cf3 = ax[2].tricontourf(
    f1[mask0],
    a1[mask0],
    np.power(10.,P_abs[mask0]/10),
    cmap="viridis",
        levels=100
)

fig.colorbar(cf1, ax=ax[1], label="Relative response (dB)")
fig.colorbar(cf2, ax=ax[0], label="Absolute response (dB)")

for axis in ax:
    axis.set_xlim(f1.min(), f1.max())
    axis.set_ylim(a1.min(), a1.max())
    axis.set_xlabel("Frequency (MHz)")

ax[0].set_ylabel("Input amplitude (mV)")
ax[0].set_title("Absolute response")
ax[1].set_title("Relative response")

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "CSI2.png", dpi=200)
plt.show()
plt.close(fig)


# ============================================================
# NEW: Peak-to-peak phase-dependence map
# ============================================================
# For every (amplitude, frequency) pair:
# ptp = max(P_interact over phase) - min(P_interact over phase).
#
# Only use groups with all expected phases to avoid treating
# incomplete phase sweeps as having zero phase dependence.

phase_df = pd.DataFrame(
    {
        "amplitude": a1,
        "frequency": f1,
        "phase": p,
        "interaction": P_interact,
    }
)

# Round phase values for grouping robustly.
phase_df["phase_group"] = phase_df["phase"].round(3)

expected_phases = np.round(p_unique, 3)
ptp_rows = []

for (amp, freq), group in phase_df.groupby(["amplitude", "frequency"]):
    observed = np.sort(group["phase_group"].unique())

    if not np.all(np.isin(expected_phases, observed)):
        continue

    # One value per phase. Average duplicates if present.
    values = (
        group.groupby("phase_group")["interaction"]
        .mean()
        .reindex(expected_phases)
        .to_numpy()
    )

    if not np.all(np.isfinite(values)):
        continue

    ptp_rows.append(
        {
            "amplitude": amp,
            "frequency": freq,
            "ptp_phase": np.ptp(values),
        }
    )

ptp_df = pd.DataFrame(ptp_rows)

if len(ptp_df) >= 3:
    max_ptp = ptp_df["ptp_phase"].max()
    levels_ptp = np.linspace(-45, 0, 100)

    fig, ax = plt.subplots(figsize=(8, 6))
    ptp_df["ptp_phase_db"] = ptp_df["ptp_phase"]
    #ptp_df["ptp_phase_db"] = np.clip(ptp_df["ptp_phase_db"], -45, 0)
    cf = ax.tricontourf(
        ptp_df["frequency"],
        ptp_df["amplitude"],
        ptp_df["ptp_phase_db"],
        levels=100,
        cmap="viridis",
    )

    fig.colorbar(
        cf,
        ax=ax,
        label="Peak-to-peak phase dependence (interaction units)",
    )
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Input amplitude (mV)")
    ax.set_title("Peak-to-peak phase dependence")

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "ptp_phase_dependence.png", dpi=200)
    plt.close(fig)

    print("Saved ptp map using", len(ptp_df), "complete phase groups.")
else:
    print(
        "Not enough complete phase groups for a ptp contour map:",
        len(ptp_df),
    )


# ============================================================
# Original phase animation, three maps
# ============================================================


max_interaction = np.nanmax(np.abs(P_interact))
if not np.isfinite(max_interaction) or max_interaction == 0:
    max_interaction = 1.0
levels1 = np.linspace(-max_interaction, max_interaction, 201)

fig, ax = plt.subplots(1, 3, figsize=(15, 5))
cb = False


def update(frame):
    global cb

    for axis in ax:
        axis.clear()

    phase_value = p_unique[frame]
    mask = np.isclose(p, phase_value, atol=1e-3)

    if mask.sum() < 3:
        for axis in ax:
            axis.set_title(f"Phase = {phase_value:.2f} rad (insufficient points)")
        return []

    cf = ax[0].tricontourf(
        f1[mask],
        a1[mask],
        P_interact[mask],
        levels=levels1,
        cmap="seismic",
    )

    cf1 = ax[1].tricontourf(
        f1[mask],
        a1[mask],
        dP_rel[mask],
        levels=levels2,
        cmap="seismic",
    )

    cf2 = ax[2].tricontourf(
        f1[mask],
        a1[mask],
        P_12[mask],
        levels=levels3,
        cmap="viridis",
    )
    fig.suptitle(f"Phase = {phase_value:.3f} rad")
    if not cb:
        fig.colorbar(cf, ax=ax[0])
        fig.colorbar(cf1, ax=ax[1])
        fig.colorbar(cf2, ax=ax[2])
        fig.tight_layout()
        cb = True
        for axis in ax:
            axis.set_xlabel("Frequency (MHz)")
            axis.set_ylabel("Input amplitude (mV)")
            axis.set_xlim(f1.min(), f1.max())
            axis.set_ylim(a1.min(), a1.max())
            ax[0].set_title("Interaction")
            ax[1].set_title("Relative response")
            ax[2].set_title("Absolute response")

    return []


ani = FuncAnimation(fig, update, frames=len(p_unique), interval=700, blit=False)
ani.save(OUTPUT_DIR / "ania2f2.gif", dpi=120)
plt.close(fig)


# ============================================================
# Fourier harmonics of the phase dependence
# ============================================================
# Use the phase samples actually present. The FFT coefficients
# are computed only for amplitude-frequency pairs with a complete
# phase sweep and at least three distinct phases.

Nphi = len(p_unique)

if Nphi >= 3:
    reference_mask = np.isclose(p, p_unique[0], atol=1e-3)
    reference_pairs = list(zip(a1[reference_mask], f1[reference_mask]))

    fourier_rows = []

    for amp, freq in reference_pairs:
        mask = np.isclose(a1, amp) & np.isclose(f1, freq)

        group = pd.DataFrame(
            {
                "phase": p[mask],
                "interaction": P_interact[mask],
            }
        )

        signal = (
            group.groupby(group["phase"].round(3))["interaction"]
            .mean()
            .reindex(expected_phases)
            .to_numpy()
        )

        if len(signal) != Nphi or not np.all(np.isfinite(signal)):
            continue

        F = fft(signal) / Nphi

        fourier_rows.append(
            {
                "amplitude": amp,
                "frequency": freq,
                "c0": F[0],
                "c1": F[1],
                "c2": F[2] if Nphi > 2 else np.nan,
            }
        )

    fourier_df = pd.DataFrame(fourier_rows)

    if not fourier_df.empty:
        for col in ["c1", "c2"]:
            scale = np.nanmax(np.abs(fourier_df[col]))
            if scale > 0:
                fourier_df[col] /= scale

        fig, ax = plt.subplots(2, 3, figsize=(12, 7))
        ax = ax.flatten()

        plots = [
            (fourier_df["c1"].apply(np.real), r"$\Re(c_1)$"),
            (fourier_df["c1"].apply(np.imag), r"$\Im(c_1)$"),
            (fourier_df["c2"].apply(np.real), r"$\Re(c_2)$"),
            (fourier_df["c2"].apply(np.imag), r"$\Im(c_2)$"),
            (np.angle(fourier_df["c1"]) / np.pi, r"$\arg(c_1)/\pi$"),
            (np.angle(fourier_df["c2"]) / np.pi*0.99, r"$\arg(c_2)/\pi$"),
        ]
        print(np.angle(fourier_df["c2"]).max() / np.pi)

        for axis, (values, title) in zip(ax, plots):
            cf = axis.tricontourf(
                fourier_df["frequency"],
                fourier_df["amplitude"],
                values,
                levels=np.linspace(-1, 1, 101),
                cmap="Spectral",
            )
            fig.colorbar(cf, ax=axis)
            axis.set_title(title)
            axis.set_xlabel("Frequency (MHz)")
            axis.set_ylabel("Input amplitude (mV)")

        fig.tight_layout()
        fig.savefig(OUTPUT_DIR / "Fourier_4f1_f2.png", dpi=200)
        plt.show()
        plt.close(fig)

        fourier_df.to_csv(OUTPUT_DIR / "phase_harmonics.csv", index=False)

print("All plots saved in:", OUTPUT_DIR.resolve())
