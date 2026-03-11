"""
plot_compton_simulation.py
──────────────────────────────────────────────────────────────────────────────
Auto-plotting script for Geant4 Compton-scattering simulation CSV outputs.

USAGE
-----
    python plot_compton_simulation.py                        # uses ./  as data dir
    python plot_compton_simulation.py /path/to/csv/folder   # explicit folder
    python plot_compton_simulation.py --save                 # save PNGs instead of showing
    python plot_compton_simulation.py /path/to/csv --save   # save to that folder

FILE NAMING CONVENTION
----------------------
The script auto-detects any of these CSV types (base name, optional suffix like --1):
    angular_distribution[*].csv
    detector_spectrum[*].csv
    phantom_dose_distribution[*].csv
    photon_spectrum[*].csv
    volume_edep[*].csv
    volume_entries[*].csv

Files sharing the same base name but different suffixes are treated as separate
run configurations and overlaid / compared automatically.

PLOTS GENERATED
---------------
Per file type:
  angular_distribution   → θ histogram | φ histogram | θ vs φ scatter | θ by volume
  detector_spectrum      → Entry energy spectrum per volume | XY hit map per volume
  phantom_dose_distribution → XY dose map | Z dose profile | edep histogram by particle
  photon_spectrum        → Energy spectrum per volume | XY photon positions
  volume_edep            → edep per volume | edep per particle type | cumulative edep
  volume_entries         → Entry energy per volume | particle fractions | process breakdown
──────────────────────────────────────────────────────────────────────────────
"""

import sys
import os
import glob
import re
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator

# ─── Aesthetics ──────────────────────────────────────────────────────────────
matplotlib.rcParams.update({
    "figure.dpi": 130,
    "font.family": "DejaVu Sans",
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

PALETTE = [
    "#2E86AB", "#E84855", "#3BB273", "#F4A261", "#9B5DE5",
    "#F72585", "#4CC9F0", "#FB8500", "#8AC926", "#FF595E",
]

VOLUME_ORDER = ["Cherenkov", "BremsTarget", "Collimator", "Phantom",
                "Scintillator1", "DWC1", "DWC2", "Scintillator2"]

PARTICLE_MARKERS = {"e-": "o", "e+": "s", "gamma": "^", "proton": "D", "neutron": "v"}


# ─── Helpers ─────────────────────────────────────────────────────────────────

def _color(i): return PALETTE[i % len(PALETTE)]

def _log_bins(series, n=50):
    lo = max(series[series > 0].min(), 1e-4)
    hi = series.max()
    if lo >= hi:
        return np.linspace(lo, hi + 1, n)
    return np.logspace(np.log10(lo), np.log10(hi), n)

def _auto_bins(series, n=50):
    return np.linspace(series.min(), series.max(), n)

def _label(run_tag):
    """Convert a run tag like '' or '--1' into a readable label."""
    if run_tag in ("", "0"):
        return "Run A"
    return f"Run {run_tag.lstrip('-')}"

def _vline_stats(ax, series, color, lw=1.2):
    mu = series.mean()
    ax.axvline(mu, color=color, linestyle="--", linewidth=lw,
               label=f"μ = {mu:.2f}")

def savefig_or_show(fig, save_dir, filename, save):
    try:
        fig.tight_layout()
    except Exception:
        pass
    if save:
        out = Path(save_dir) / filename
        fig.savefig(out, bbox_inches="tight")
        print(f"  Saved → {out}")
        plt.close(fig)
    else:
        plt.show()


# ─── Loaders ─────────────────────────────────────────────────────────────────

BASE_NAMES = [
    "angular_distribution",
    "detector_spectrum",
    "phantom_dose_distribution",
    "photon_spectrum",
    "volume_edep",
    "volume_entries",
]

def discover_files(data_dir):
    """Return {base_name: {tag: DataFrame}} for every recognised CSV in data_dir."""
    data_dir = Path(data_dir)
    result = defaultdict(dict)
    for csv_path in sorted(data_dir.glob("*.csv")):
        name = csv_path.stem  # e.g.  "angular_distribution--1"
        for base in BASE_NAMES:
            if name == base or name.startswith(base + "-") or name.startswith(base + "_"):
                tag = name[len(base):]  # "", "--1", "_2", …
                try:
                    df = pd.read_csv(csv_path, on_bad_lines="skip", engine="python")
                    df.columns = df.columns.str.strip()
                    # drop rows where numeric columns got corrupted strings
                    for col in df.columns:
                        if col not in ("volume", "particle", "process", "fromVolume"):
                            df[col] = pd.to_numeric(df[col], errors="coerce")
                    df = df.dropna()
                    result[base][tag] = df
                    print(f"  Loaded {csv_path.name}  →  {base!r}  tag={tag!r}  "
                          f"({len(df)} rows)")
                except Exception as exc:
                    print(f"  WARN: could not load {csv_path.name}: {exc}")
                break
    return result


# ─── Plot functions ───────────────────────────────────────────────────────────

# ── 1. Angular distribution ──────────────────────────────────────────────────

def plot_angular_distribution(runs: dict, save_dir, save):
    """θ histogram | φ histogram | θ-φ scatter | θ by volume."""

    fig = plt.figure(figsize=(14, 10))
    fig.suptitle("Angular Distribution of Compton Photons", fontsize=13, fontweight="bold")
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)

    ax_theta = fig.add_subplot(gs[0, 0])
    ax_phi   = fig.add_subplot(gs[0, 1])
    ax_scat  = fig.add_subplot(gs[1, 0])
    ax_vol   = fig.add_subplot(gs[1, 1])

    for i, (tag, df) in enumerate(runs.items()):
        col = _color(i)
        lbl = _label(tag)
        # filter to photons only if 'particle' column exists
        dg = df[df["particle"] == "gamma"] if "particle" in df.columns else df
        if dg.empty:
            dg = df  # fall back to all

        theta = dg["theta_deg"]
        phi   = dg["phi_deg"]

        # θ histogram
        bins_t = _auto_bins(theta, 40)
        ax_theta.hist(theta, bins=bins_t, color=col, alpha=0.7, label=lbl, density=False)
        _vline_stats(ax_theta, theta, col)

        # φ histogram
        bins_p = _auto_bins(phi, 40)
        ax_phi.hist(phi, bins=bins_p, color=col, alpha=0.7, label=lbl)

        # θ-φ scatter
        ax_scat.scatter(phi, theta, s=10, alpha=0.5, color=col, label=lbl)

    ax_theta.set(xlabel="Scattering angle θ (°)", ylabel="Counts",
                 title="Compton Scattering Angle (θ)")
    ax_theta.legend()

    ax_phi.set(xlabel="Azimuthal angle φ (°)", ylabel="Counts",
               title="Azimuthal Angle Distribution (φ)")
    ax_phi.legend()

    ax_scat.set(xlabel="φ (°)", ylabel="θ (°)",
                title="θ vs φ Scatter")
    ax_scat.legend()

    # θ per volume (all runs combined)
    all_df = pd.concat(list(runs.values()), ignore_index=True)
    if "volume" in all_df.columns:
        volumes = all_df["volume"].unique()
        for j, vol in enumerate(volumes):
            sub = all_df[all_df["volume"] == vol]
            if "particle" in sub.columns:
                sub = sub[sub["particle"] == "gamma"]
            if sub.empty:
                continue
            ax_vol.hist(sub["theta_deg"], bins=_auto_bins(sub["theta_deg"], 30),
                        alpha=0.6, label=vol, color=_color(j))
        ax_vol.set(xlabel="θ (°)", ylabel="Counts",
                   title="θ Distribution by Volume")
        ax_vol.legend()
    else:
        ax_vol.axis("off")

    savefig_or_show(fig, save_dir, "01_angular_distribution.png", save)


# ── 2. Detector spectrum ──────────────────────────────────────────────────────

def plot_detector_spectrum(runs: dict, save_dir, save):
    """Entry energy spectrum per volume | XY hit maps."""

    all_df  = pd.concat(list(runs.values()), ignore_index=True)
    volumes = [v for v in VOLUME_ORDER if v in all_df["volume"].unique()]
    if not volumes:
        volumes = list(all_df["volume"].unique())

    n_vol = len(volumes)
    ncols = min(3, n_vol)
    nrows = int(np.ceil(n_vol / ncols))

    # ── Energy spectra per volume ──
    fig1, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    fig1.suptitle("Entry Energy Spectrum per Detector Volume", fontsize=13, fontweight="bold")
    axes = np.array(axes).flatten()

    for idx, vol in enumerate(volumes):
        ax = axes[idx]
        for i, (tag, df) in enumerate(runs.items()):
            sub = df[df["volume"] == vol]
            if sub.empty:
                continue
            en = sub["energy_MeV"]
            bins = _log_bins(en)
            ax.hist(en, bins=bins, alpha=0.65, label=_label(tag), color=_color(i))
        ax.set_xscale("log")
        ax.set(xlabel="Energy (MeV)", ylabel="Counts", title=vol)
        ax.legend()

    for ax in axes[n_vol:]:
        ax.axis("off")

    savefig_or_show(fig1, save_dir, "02a_detector_spectrum_energy.png", save)

    # ── XY hit maps per volume ──
    fig2, axes2 = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows))
    fig2.suptitle("XY Hit Positions per Detector Volume", fontsize=13, fontweight="bold")
    axes2 = np.array(axes2).flatten()

    for idx, vol in enumerate(volumes):
        ax = axes2[idx]
        for i, (tag, df) in enumerate(runs.items()):
            sub = df[df["volume"] == vol]
            if sub.empty or "x_mm" not in sub.columns:
                continue
            ax.scatter(sub["x_mm"], sub["y_mm"], s=8, alpha=0.5,
                       label=_label(tag), color=_color(i))
        ax.set(xlabel="x (mm)", ylabel="y (mm)", title=f"{vol} — Hit Map",
               aspect="equal")
        ax.legend()

    for ax in axes2[n_vol:]:
        ax.axis("off")

    savefig_or_show(fig2, save_dir, "02b_detector_spectrum_xy.png", save)


# ── 3. Phantom dose distribution ─────────────────────────────────────────────

def plot_phantom_dose(runs: dict, save_dir, save):
    """XY dose map | Z dose profile | edep histogram by particle."""

    fig = plt.figure(figsize=(15, 10))
    fig.suptitle("Phantom Dose Distribution", fontsize=13, fontweight="bold")
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.38)

    ax_xy   = fig.add_subplot(gs[0, 0])
    ax_xz   = fig.add_subplot(gs[0, 1])
    ax_z    = fig.add_subplot(gs[0, 2])
    ax_edep = fig.add_subplot(gs[1, 0])
    ax_part = fig.add_subplot(gs[1, 1])
    ax_3d   = fig.add_subplot(gs[1, 2])  # optional 2-D heatmap

    for i, (tag, df) in enumerate(runs.items()):
        col = _color(i)
        lbl = _label(tag)

        # XY scatter coloured by edep
        sc = ax_xy.scatter(df["x_mm"], df["y_mm"], c=df["edep_MeV"],
                           cmap="plasma", s=12, alpha=0.6, label=lbl)

        # XZ (side view)
        ax_xz.scatter(df["z_mm"], df["x_mm"], c=df["edep_MeV"],
                      cmap="viridis", s=12, alpha=0.5)

        # Z profile
        ax_z.hist(df["z_mm"], bins=_auto_bins(df["z_mm"], 40),
                  alpha=0.65, color=col, label=lbl)

        # edep histogram
        ax_edep.hist(df["edep_MeV"], bins=_auto_bins(df["edep_MeV"], 40),
                     alpha=0.65, color=col, label=lbl)

    # edep by particle
    all_df = pd.concat(list(runs.values()), ignore_index=True)
    if "particle" in all_df.columns:
        for j, pt in enumerate(all_df["particle"].unique()):
            sub = all_df[all_df["particle"] == pt]
            ax_part.hist(sub["edep_MeV"], bins=_auto_bins(sub["edep_MeV"], 40),
                         alpha=0.65, label=pt, color=_color(j))
    ax_part.set(xlabel="Energy deposition (MeV)", ylabel="Counts",
                title="edep by Particle Type")
    ax_part.legend()

    # 2D edep heatmap in XY
    hx = np.linspace(all_df["x_mm"].min(), all_df["x_mm"].max(), 40)
    hy = np.linspace(all_df["y_mm"].min(), all_df["y_mm"].max(), 40)
    H, xedge, yedge = np.histogram2d(all_df["x_mm"], all_df["y_mm"],
                                      bins=[hx, hy], weights=all_df["edep_MeV"])
    im = ax_3d.pcolormesh(xedge, yedge, H.T, cmap="hot")
    fig.colorbar(im, ax=ax_3d, label="ΣEdep (MeV)")

    ax_xy.set(xlabel="x (mm)", ylabel="y (mm)", title="XY Dose Map (colour = edep)")
    ax_xy.legend()
    ax_xz.set(xlabel="z (mm)", ylabel="x (mm)", title="ZX Side View")
    ax_z.set(xlabel="z (mm)", ylabel="Counts", title="Dose along Beam Axis (z)")
    ax_z.legend()
    ax_edep.set(xlabel="Energy deposition (MeV)", ylabel="Counts",
                title="edep Histogram")
    ax_edep.legend()
    ax_3d.set(xlabel="x (mm)", ylabel="y (mm)", title="2D edep Heatmap (XY)")

    savefig_or_show(fig, save_dir, "03_phantom_dose.png", save)


# ── 4. Photon spectrum ────────────────────────────────────────────────────────

def plot_photon_spectrum(runs: dict, save_dir, save):
    """Photon energy spectrum per volume | XY positions coloured by energy."""

    all_df  = pd.concat(list(runs.values()), ignore_index=True)
    volumes = [v for v in VOLUME_ORDER if v in all_df["volume"].unique()]
    if not volumes:
        volumes = list(all_df["volume"].unique())

    fig, axes = plt.subplots(2, len(volumes), figsize=(5 * len(volumes), 9))
    if len(volumes) == 1:
        axes = axes.reshape(2, 1)
    fig.suptitle("Photon Spectrum per Volume", fontsize=13, fontweight="bold")

    for col_idx, vol in enumerate(volumes):
        ax_e  = axes[0, col_idx]
        ax_xy = axes[1, col_idx]

        for i, (tag, df) in enumerate(runs.items()):
            sub = df[df["volume"] == vol]
            if sub.empty:
                continue
            bins = _log_bins(sub["energy_MeV"])
            ax_e.hist(sub["energy_MeV"], bins=bins, alpha=0.65,
                      label=_label(tag), color=_color(i))

            if "x_mm" in sub.columns:
                ax_xy.scatter(sub["x_mm"], sub["y_mm"], c=sub["energy_MeV"],
                              cmap="coolwarm", norm=LogNorm(), s=12, alpha=0.6)

        ax_e.set_xscale("log")
        ax_e.set(xlabel="Energy (MeV)", ylabel="Counts", title=f"{vol}\nEnergy")
        ax_e.legend()
        ax_xy.set(xlabel="x (mm)", ylabel="y (mm)", title=f"{vol}\nXY (colour=E)",
                  aspect="equal")

    savefig_or_show(fig, save_dir, "04_photon_spectrum.png", save)


# ── 5. Volume energy deposition ───────────────────────────────────────────────

def plot_volume_edep(runs: dict, save_dir, save):
    """edep per volume (stacked bars) | edep per particle | cumulative edep."""

    fig = plt.figure(figsize=(15, 10))
    fig.suptitle("Energy Deposition by Volume & Particle Type", fontsize=13, fontweight="bold")
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.38)

    ax_vol  = fig.add_subplot(gs[0, :2])
    ax_part = fig.add_subplot(gs[0, 2])
    ax_hist = fig.add_subplot(gs[1, 0])
    ax_run  = fig.add_subplot(gs[1, 1])
    ax_cum  = fig.add_subplot(gs[1, 2])

    all_df  = pd.concat(list(runs.values()), ignore_index=True)
    volumes = [v for v in VOLUME_ORDER if v in all_df["volume"].unique()]
    if not volumes:
        volumes = list(all_df["volume"].unique())
    particles = list(all_df["particle"].unique()) if "particle" in all_df.columns else []

    # Stacked bar: total edep per volume × particle
    x = np.arange(len(volumes))
    bottom = np.zeros(len(volumes))
    for j, pt in enumerate(particles):
        heights = []
        for vol in volumes:
            sub = all_df[(all_df["volume"] == vol) & (all_df["particle"] == pt)]
            heights.append(sub["edep_MeV"].sum())
        ax_vol.bar(x, heights, bottom=bottom, label=pt, color=_color(j), alpha=0.85)
        bottom += np.array(heights)

    ax_vol.set_xticks(x)
    ax_vol.set_xticklabels(volumes, rotation=25, ha="right")
    ax_vol.set(ylabel="Total edep (MeV)", title="Total Energy Deposition per Volume")
    ax_vol.legend()

    # Per-particle histogram
    for j, pt in enumerate(particles):
        sub = all_df[all_df["particle"] == pt]
        ax_part.hist(sub["edep_MeV"], bins=_auto_bins(sub["edep_MeV"], 40),
                     alpha=0.65, label=pt, color=_color(j))
    ax_part.set(xlabel="edep (MeV)", ylabel="Counts", title="edep by Particle")
    ax_part.legend()

    # Log-scale edep histogram per volume
    for j, vol in enumerate(volumes):
        sub = all_df[all_df["volume"] == vol]
        if sub.empty or sub["edep_MeV"].max() <= 0:
            continue
        ax_hist.hist(sub["edep_MeV"], bins=_log_bins(sub["edep_MeV"], 40),
                     alpha=0.5, label=vol, color=_color(j))
    ax_hist.set_xscale("log")
    ax_hist.set(xlabel="edep (MeV)", ylabel="Counts", title="edep Distribution (log x)")
    ax_hist.legend()

    # Run comparison: mean edep per volume
    for i, (tag, df) in enumerate(runs.items()):
        means = [df[df["volume"] == v]["edep_MeV"].mean() for v in volumes]
        ax_run.bar(x + i * 0.35, means, width=0.35, label=_label(tag), color=_color(i), alpha=0.8)
    ax_run.set_xticks(x + 0.175)
    ax_run.set_xticklabels(volumes, rotation=25, ha="right")
    ax_run.set(ylabel="Mean edep (MeV)", title="Mean edep per Volume (Run Comparison)")
    ax_run.legend()

    # Cumulative edep (sorted) per run
    for i, (tag, df) in enumerate(runs.items()):
        sorted_edep = np.sort(df["edep_MeV"].values)
        cum = np.cumsum(sorted_edep) / sorted_edep.sum()
        ax_cum.plot(sorted_edep, cum, color=_color(i), label=_label(tag))
    ax_cum.set_xscale("log")
    ax_cum.set(xlabel="edep (MeV)", ylabel="Cumulative fraction",
               title="Cumulative edep Distribution")
    ax_cum.legend()

    savefig_or_show(fig, save_dir, "05_volume_edep.png", save)


# ── 6. Volume entries ─────────────────────────────────────────────────────────

def plot_volume_entries(runs: dict, save_dir, save):
    """Entry energy per volume | particle fractions | process breakdown."""

    fig = plt.figure(figsize=(15, 10))
    fig.suptitle("Particle Entries into Detector Volumes", fontsize=13, fontweight="bold")
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.4)

    ax_cnt  = fig.add_subplot(gs[0, 0])
    ax_en   = fig.add_subplot(gs[0, 1])
    ax_part = fig.add_subplot(gs[0, 2])
    ax_proc = fig.add_subplot(gs[1, 0])
    ax_xy   = fig.add_subplot(gs[1, 1])
    ax_run  = fig.add_subplot(gs[1, 2])

    all_df  = pd.concat(list(runs.values()), ignore_index=True)
    volumes = [v for v in VOLUME_ORDER if v in all_df["volume"].unique()]
    if not volumes:
        volumes = list(all_df["volume"].unique())
    x = np.arange(len(volumes))

    # Entry counts per volume per run
    for i, (tag, df) in enumerate(runs.items()):
        counts = [len(df[df["volume"] == v]) for v in volumes]
        ax_cnt.bar(x + i * 0.35, counts, width=0.35, label=_label(tag),
                   color=_color(i), alpha=0.8)
    ax_cnt.set_xticks(x + 0.175)
    ax_cnt.set_xticklabels(volumes, rotation=25, ha="right")
    ax_cnt.set(ylabel="Entry count", title="Particle Entries per Volume")
    ax_cnt.legend()

    # Entry energy distribution per volume (log)
    for j, vol in enumerate(volumes):
        sub = all_df[all_df["volume"] == vol]
        if sub.empty or "energy_MeV" not in sub.columns:
            continue
        en = sub["energy_MeV"]
        if en.max() <= 0:
            continue
        ax_en.hist(en, bins=_log_bins(en, 40), alpha=0.5, label=vol, color=_color(j))
    ax_en.set_xscale("log")
    ax_en.set(xlabel="Energy (MeV)", ylabel="Counts", title="Entry Energy by Volume")
    ax_en.legend()

    # Particle type fractions per volume (stacked bar)
    particles = list(all_df["particle"].unique()) if "particle" in all_df.columns else []
    bottom = np.zeros(len(volumes))
    for j, pt in enumerate(particles):
        heights = [len(all_df[(all_df["volume"] == v) & (all_df["particle"] == pt)])
                   for v in volumes]
        ax_part.bar(x, heights, bottom=bottom, label=pt, color=_color(j), alpha=0.85)
        bottom += np.array(heights)
    ax_part.set_xticks(x)
    ax_part.set_xticklabels(volumes, rotation=25, ha="right")
    ax_part.set(ylabel="Counts", title="Particle Type per Volume")
    ax_part.legend()

    # Process breakdown (if column exists)
    if "process" in all_df.columns:
        procs = all_df["process"].value_counts()
        ax_proc.barh(procs.index, procs.values, color=PALETTE[:len(procs)], alpha=0.8)
        ax_proc.set(xlabel="Count", title="Physics Process Breakdown")
    else:
        ax_proc.axis("off")

    # XY entry positions (all volumes combined)
    if "x_mm" in all_df.columns:
        for j, vol in enumerate(volumes):
            sub = all_df[all_df["volume"] == vol]
            ax_xy.scatter(sub["x_mm"], sub["y_mm"], s=8, alpha=0.5,
                          label=vol, color=_color(j))
        ax_xy.set(xlabel="x (mm)", ylabel="y (mm)", title="XY Entry Positions (all volumes)")
        ax_xy.legend()
    else:
        ax_xy.axis("off")

    # fromVolume transition matrix (if column exists)
    if "fromVolume" in all_df.columns:
        fv = all_df["fromVolume"].value_counts().head(10)
        ax_run.barh(fv.index, fv.values, color=PALETTE[:len(fv)], alpha=0.8)
        ax_run.set(xlabel="Count", title="Origin Volume of Entering Particles")
    else:
        # fallback: mean entry energy per volume per run comparison
        for i, (tag, df) in enumerate(runs.items()):
            if "energy_MeV" not in df.columns:
                continue
            means = [df[df["volume"] == v]["energy_MeV"].mean() for v in volumes]
            ax_run.bar(x + i * 0.35, means, width=0.35, label=_label(tag),
                       color=_color(i), alpha=0.8)
        ax_run.set_xticks(x + 0.175)
        ax_run.set_xticklabels(volumes, rotation=25, ha="right")
        ax_run.set(ylabel="Mean energy (MeV)", title="Mean Entry Energy (Run Comparison)")
        ax_run.legend()

    savefig_or_show(fig, save_dir, "06_volume_entries.png", save)


# ── 7. Cross-run summary dashboard ───────────────────────────────────────────

def plot_summary_dashboard(file_data: dict, save_dir, save):
    """One-page overview comparing the two runs across all key observables."""

    fig = plt.figure(figsize=(16, 12))
    fig.suptitle(
        "Geant4 Compton Simulation — Summary Dashboard",
        fontsize=14, fontweight="bold"
    )
    gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.55, wspace=0.4)

    # Row 0: Angular θ | Photon spectrum (Phantom) | Phantom edep Z
    ax_theta = fig.add_subplot(gs[0, 0])
    ax_phot  = fig.add_subplot(gs[0, 1])
    ax_zprof = fig.add_subplot(gs[0, 2])
    ax_proc  = fig.add_subplot(gs[0, 3])

    # Row 1: Volume total edep | Particle fractions | Entry energy Phantom
    ax_edvol = fig.add_subplot(gs[1, :2])
    ax_pfrac = fig.add_subplot(gs[1, 2])
    ax_enph  = fig.add_subplot(gs[1, 3])

    # Row 2: XY phantom hits | Klein-Nishina overlay
    ax_xyph  = fig.add_subplot(gs[2, 0])
    ax_kn    = fig.add_subplot(gs[2, 1:3])
    ax_cum   = fig.add_subplot(gs[2, 3])

    # ── Populate ──

    # θ histogram
    if "angular_distribution" in file_data:
        for i, (tag, df) in enumerate(file_data["angular_distribution"].items()):
            dg = df[df["particle"] == "gamma"] if "particle" in df.columns else df
            ax_theta.hist(dg["theta_deg"], bins=40, alpha=0.65,
                          label=_label(tag), color=_color(i))
        ax_theta.set(xlabel="θ (°)", ylabel="Counts", title="Scattering Angle θ")
        ax_theta.legend()

    # Photon spectrum (Phantom volume)
    if "photon_spectrum" in file_data:
        for i, (tag, df) in enumerate(file_data["photon_spectrum"].items()):
            sub = df[df["volume"] == "Phantom"] if "volume" in df.columns else df
            if sub.empty:
                sub = df
            if not sub.empty:
                ax_phot.hist(sub["energy_MeV"], bins=_log_bins(sub["energy_MeV"], 40),
                             alpha=0.65, label=_label(tag), color=_color(i))
        ax_phot.set_xscale("log")
        ax_phot.set(xlabel="E (MeV)", ylabel="Counts", title="Photon Energy (Phantom)")
        ax_phot.legend()

    # Phantom Z profile
    if "phantom_dose_distribution" in file_data:
        for i, (tag, df) in enumerate(file_data["phantom_dose_distribution"].items()):
            ax_zprof.hist(df["z_mm"], bins=40, alpha=0.65,
                          label=_label(tag), color=_color(i))
        ax_zprof.set(xlabel="z (mm)", ylabel="Counts", title="Dose Along z")
        ax_zprof.legend()

    # Process breakdown
    if "volume_entries" in file_data:
        all_ent = pd.concat(list(file_data["volume_entries"].values()), ignore_index=True)
        if "process" in all_ent.columns:
            procs = all_ent["process"].value_counts().head(8)
            ax_proc.barh(procs.index, procs.values, color=PALETTE[:len(procs)], alpha=0.8)
            ax_proc.set(xlabel="Count", title="Physics Processes")
        else:
            ax_proc.axis("off")

    # Volume total edep (stacked bar)
    if "volume_edep" in file_data:
        all_edep = pd.concat(list(file_data["volume_edep"].values()), ignore_index=True)
        volumes = [v for v in VOLUME_ORDER if v in all_edep["volume"].unique()]
        if not volumes:
            volumes = list(all_edep["volume"].unique())
        x = np.arange(len(volumes))
        particles = list(all_edep["particle"].unique()) if "particle" in all_edep.columns else []
        bottom = np.zeros(len(volumes))
        for j, pt in enumerate(particles):
            heights = [all_edep[(all_edep["volume"] == v) &
                                (all_edep["particle"] == pt)]["edep_MeV"].sum()
                       for v in volumes]
            ax_edvol.bar(x, heights, bottom=bottom, label=pt, color=_color(j), alpha=0.85)
            bottom += np.array(heights)
        ax_edvol.set_xticks(x)
        ax_edvol.set_xticklabels(volumes, rotation=20, ha="right")
        ax_edvol.set(ylabel="Total edep (MeV)", title="Energy Deposition per Volume")
        ax_edvol.legend()

    # Particle fractions (pie)
    if "volume_entries" in file_data:
        all_ent = pd.concat(list(file_data["volume_entries"].values()), ignore_index=True)
        if "particle" in all_ent.columns:
            counts = all_ent["particle"].value_counts()
            ax_pfrac.pie(counts.values, labels=counts.index,
                         colors=PALETTE[:len(counts)],
                         autopct="%1.0f%%", startangle=90)
            ax_pfrac.set_title("Particle Type Fractions")

    # Entry energy – Phantom
    if "volume_entries" in file_data:
        all_ent = pd.concat(list(file_data["volume_entries"].values()), ignore_index=True)
        sub = all_ent[all_ent["volume"] == "Phantom"] if "volume" in all_ent.columns else all_ent
        if not sub.empty and "energy_MeV" in sub.columns:
            ax_enph.hist(sub["energy_MeV"], bins=_log_bins(sub["energy_MeV"], 40),
                         color=_color(0), alpha=0.8)
            ax_enph.set_xscale("log")
            ax_enph.set(xlabel="E (MeV)", ylabel="Counts", title="Phantom Entry Energy")

    # XY phantom hit map
    if "phantom_dose_distribution" in file_data:
        for i, (tag, df) in enumerate(file_data["phantom_dose_distribution"].items()):
            ax_xyph.scatter(df["x_mm"], df["y_mm"], c=df["edep_MeV"],
                            cmap="plasma", s=14, alpha=0.6, label=_label(tag))
        ax_xyph.set(xlabel="x (mm)", ylabel="y (mm)", title="Phantom XY (colour=edep)")
        ax_xyph.legend()

    # Klein-Nishina curve overlay on angular distribution
    theta_kn = np.linspace(0, np.pi, 500)
    # differential cross section (unnormalised, r0=1) at a few energies
    for ei, (E_MeV, lbl_kn) in enumerate([(1, "1 MeV"), (10, "10 MeV"), (100, "100 MeV")]):
        k = E_MeV / 0.511  # E / m_e c²
        eps = 1 / (1 + k * (1 - np.cos(theta_kn)))
        dsig = (eps ** 2) * (eps + 1/eps - np.sin(theta_kn) ** 2)
        dsig /= dsig.max()
        ax_kn.plot(np.degrees(theta_kn), dsig, color=_color(ei + 2), label=lbl_kn)

    if "angular_distribution" in file_data:
        # overlay normalised histogram
        all_ang = pd.concat(list(file_data["angular_distribution"].values()), ignore_index=True)
        dg = all_ang[all_ang["particle"] == "gamma"] if "particle" in all_ang.columns else all_ang
        if not dg.empty:
            ax_kn.hist(dg["theta_deg"], bins=40, density=True,
                       alpha=0.3, color="gray", label="Sim data")

    ax_kn.set(xlabel="θ (°)", ylabel="dσ/dΩ (normalised)", xlim=(0, 180),
              title="Klein-Nishina Cross-section vs Simulated θ")
    ax_kn.legend()

    # Cumulative edep
    if "volume_edep" in file_data:
        for i, (tag, df) in enumerate(file_data["volume_edep"].items()):
            sorted_e = np.sort(df["edep_MeV"].values)
            cum = np.cumsum(sorted_e) / sorted_e.sum()
            ax_cum.plot(sorted_e, cum, color=_color(i), label=_label(tag))
        ax_cum.set_xscale("log")
        ax_cum.set(xlabel="edep (MeV)", ylabel="Cumulative fraction",
                   title="Cumulative edep")
        ax_cum.legend()

    savefig_or_show(fig, save_dir, "00_summary_dashboard.png", save)


# ─── Main ─────────────────────────────────────────────────────────────────────

PLOT_REGISTRY = {
    "angular_distribution":      plot_angular_distribution,
    "detector_spectrum":         plot_detector_spectrum,
    "phantom_dose_distribution": plot_phantom_dose,
    "photon_spectrum":           plot_photon_spectrum,
    "volume_edep":               plot_volume_edep,
    "volume_entries":            plot_volume_entries,
}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = [a for a in sys.argv[1:] if a.startswith("--")]
    save = "--save" in flags

    data_dir = Path(args[0]) if args else Path(".")
    save_dir = Path(args[1]) if len(args) > 1 else Path(".")  # output dir separate from data

    print(f"\n{'='*60}")
    print(f"  Compton Geant4 Auto-Plotter")
    print(f"  Data dir : {data_dir.resolve()}")
    print(f"  Save PNGs: {save}")
    print(f"{'='*60}\n")

    file_data = discover_files(data_dir)

    if not file_data:
        print("  No recognised CSV files found. Exiting.")
        sys.exit(1)

    print(f"\nGenerating individual plots …\n")
    for base_name, plot_fn in PLOT_REGISTRY.items():
        if base_name in file_data:
            print(f"  → {base_name}")
            plot_fn(file_data[base_name], save_dir, save)

    print(f"\nGenerating summary dashboard …\n")
    plot_summary_dashboard(file_data, save_dir, save)

    print(f"\nDone. {len(file_data)} file type(s) plotted.\n")


if __name__ == "__main__":
    main()