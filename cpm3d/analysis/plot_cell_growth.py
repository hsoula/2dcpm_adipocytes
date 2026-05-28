#!/usr/bin/env python3
"""
plot_cell_growth.py — cell growth analysis for cpm3d
-----------------------------------------------------
Reads  cells.csv  produced by `cargo run --bin analyze`
and saves two figures:

  cell_growth.png          — actual & target volume trajectories per cell
                             with population mean ± 1 SD overlay

  growth_speed_vs_r.png    — instantaneous growth speed dV/dt (central
                             finite differences over adjacent snapshots)
                             plotted as a function of spherical radius r,
                             one point per cell per time step, with a
                             binned mean ± 1 SD curve

If  growth_trace.csv  is present (produced by --track-growth in simulate_life),
a third figure is produced:

  within_mcs_growth.png    — per-cell volume step-function at sub-MCS
                             resolution; x-axis is mcs + fraction-through-MCS

Usage
-----
  python analysis/plot_cell_growth.py <dir>
  python analysis/plot_cell_growth.py <dir> --out figures/
  python analysis/plot_cell_growth.py <dir> --r-bins 30 --min-points 4
  python analysis/plot_cell_growth.py <dir> --mcs-range 10 50
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# ── Geometry ──────────────────────────────────────────────────────────────────

def vol_to_radius(v):
    """Spherical approximation: r = (3V / 4π)^(1/3)"""
    v = np.asarray(v, dtype=float)
    return np.cbrt(3.0 * v / (4.0 * np.pi))


# ── Speed computation ─────────────────────────────────────────────────────────

def compute_speeds(df: pd.DataFrame, min_points: int = 3) -> pd.DataFrame:
    """
    For each continuous cell track (sigma × birth_mcs), compute the
    instantaneous growth speed at every interior snapshot using central
    finite differences:

        speed[i] = (V[i+1] - V[i-1]) / (mcs[i+1] - mcs[i-1])

    Returns a DataFrame with columns:
        sigma, mcs, volume, radius, speed
    """
    records = []
    for (sigma, birth_mcs), track in df.groupby(["sigma", "birth_mcs"], sort=False):
        track = track.sort_values("mcs")
        mcs = track["mcs"].to_numpy(dtype=float)
        vol = track["volume"].to_numpy(dtype=float)
        if len(mcs) < min_points:
            continue
        dt    = mcs[2:] - mcs[:-2]          # denominator: 2*Δt
        dv    = vol[2:] - vol[:-2]
        speed = np.where(dt > 0, dv / dt, np.nan)
        r     = vol_to_radius(vol[1:-1])
        for i in range(len(speed)):
            records.append({
                "sigma":  sigma,
                "mcs":    mcs[1:-1][i],
                "volume": vol[1:-1][i],
                "radius": r[i],
                "speed":  speed[i],
            })
    return pd.DataFrame(records) if records else pd.DataFrame(
        columns=["sigma", "mcs", "volume", "radius", "speed"]
    )


# ── Figure 1: volume trajectories ────────────────────────────────────────────

def plot_trajectories(df: pd.DataFrame, out_path: Path):
    """Actual volume + target volume per cell, with mean ± SD."""
    tracks = df.groupby(["sigma", "birth_mcs"])
    n_tracks = len(tracks)

    fig, (ax_v, ax_tv) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    cmap = cm.get_cmap("tab20", max(n_tracks, 1))

    for i, ((sigma, bm), track) in enumerate(sorted(tracks)):
        track = track.sort_values("mcs")
        c = cmap(i % 20)
        ax_v.plot(track["mcs"], track["volume"],
                  color=c, lw=0.9, alpha=0.5)
        ax_tv.plot(track["mcs"], track["target_volume"],
                   color=c, lw=0.9, alpha=0.5, ls="--")

    # Population mean ± SD at each MCS (non-dying cells only)
    g_v  = df.groupby("mcs")["volume"]
    g_tv = df.groupby("mcs")["target_volume"]
    t    = g_v.mean().index.to_numpy()

    for ax, g, label in [(ax_v, g_v, "actual volume"),
                         (ax_tv, g_tv, "target volume")]:
        mean = g.mean().to_numpy()
        std  = g.std().fillna(0).to_numpy()
        ax.plot(t, mean, color="black", lw=2.2, zorder=5, label="mean")
        ax.fill_between(t, mean - std, mean + std,
                        color="black", alpha=0.12, label="±1 SD")
        ax.set_ylabel(f"{label.capitalize()} (voxels)", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9, loc="upper left")

    ax_v.set_title(f"Volume trajectories  ({n_tracks} cell tracks)", fontsize=12)
    ax_tv.set_xlabel("MCS", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"→ {out_path}")


# ── Figure 2: growth speed vs radius ─────────────────────────────────────────

def plot_speed_vs_radius(speed_df: pd.DataFrame, out_path: Path, n_bins: int = 20):
    """
    Scatter of (radius, dV/dt) for every cell × time point,
    overlaid with a binned mean ± 1 SD curve.
    """
    sd = speed_df.dropna(subset=["speed"])
    if sd.empty:
        print("No speed data — skipping growth_speed_vs_r.png")
        return

    r     = sd["radius"].to_numpy()
    speed = sd["speed"].to_numpy()

    # Colour by sigma so each cell gets its own hue
    sigmas      = sd["sigma"].to_numpy()
    unique_sigs = np.unique(sigmas)
    cmap        = cm.get_cmap("tab20", len(unique_sigs))
    sig_idx     = {s: i for i, s in enumerate(unique_sigs)}
    colours     = [cmap(sig_idx[s] % 20) for s in sigmas]

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(r, speed, c=colours, alpha=0.25, s=8, linewidths=0,
               label="per-cell point")

    # Binned mean ± SD
    edges   = np.linspace(r.min(), r.max(), n_bins + 1)
    centres, means, stds, counts = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (r >= lo) & (r < hi)
        sub  = speed[mask]
        if len(sub) < 3:
            continue
        centres.append((lo + hi) / 2)
        means.append(sub.mean())
        stds.append(sub.std())
        counts.append(len(sub))

    if centres:
        c_arr = np.array(centres)
        m_arr = np.array(means)
        s_arr = np.array(stds)
        ax.plot(c_arr, m_arr, color="crimson", lw=2.2, zorder=5,
                label="binned mean")
        ax.fill_between(c_arr, m_arr - s_arr, m_arr + s_arr,
                        color="crimson", alpha=0.20, label="±1 SD")

    ax.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
    ax.set_xlabel(r"Radius  $r = \left(\frac{3V}{4\pi}\right)^{1/3}$  (voxels)",
                  fontsize=12)
    ax.set_ylabel(r"Growth speed  $\mathrm{d}V/\mathrm{d}t$  (voxels / MCS)",
                  fontsize=12)
    ax.set_title(f"Growth speed vs cell radius  "
                 f"({len(sd)} points, {len(unique_sigs)} cells)",
                 fontsize=12)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"→ {out_path}")


# ── Growth-trace loading and figure ──────────────────────────────────────────

def load_growth_trace(path: Path, mcs_range=None) -> pd.DataFrame:
    """
    Load growth_trace.csv (kind,mcs,sigma,area_at_event,birth_mcs,lifetime_mcs).

    Assigns a continuous time axis:
        time = mcs + event_local_index / events_in_this_mcs
    so each MCS maps to [mcs, mcs+1) and events are evenly spread within it.
    """
    df = pd.read_csv(path)
    if mcs_range is not None:
        lo, hi = mcs_range
        df = df[(df["mcs"] >= lo) & (df["mcs"] <= hi)]
    if df.empty:
        return df

    # fractional position within each MCS (based on row order in the CSV)
    df = df.reset_index(drop=True)
    local_idx  = df.groupby("mcs").cumcount()
    mcs_counts = df.groupby("mcs")["mcs"].transform("count")
    df["time"] = df["mcs"] + local_idx / mcs_counts
    return df


def plot_within_mcs_growth(gt: pd.DataFrame, out_path: Path):
    """
    Step-function volume trace for every cell at sub-MCS resolution.

    Each accepted voxel flip (grow/shrink) is one point; the step connects
    consecutive events for the same cell, giving true within-MCS dynamics.
    """
    if gt.empty:
        print("growth_trace.csv is empty — skipping within_mcs_growth.png")
        return

    unique_sigs = sorted(gt["sigma"].unique())
    cmap = cm.get_cmap("tab20", max(len(unique_sigs), 1))
    sig_color = {s: cmap(i % 20) for i, s in enumerate(unique_sigs)}

    fig, ax = plt.subplots(figsize=(12, 5))

    for sigma, track in gt.groupby("sigma"):
        track = track.sort_values("time")
        c = sig_color[sigma]
        # drawstyle='steps-post' draws a horizontal segment until the next event
        ax.plot(track["time"], track["area_at_event"],
                color=c, lw=0.7, alpha=0.7, drawstyle="steps-post",
                label=f"σ={sigma}")

    # mark MCS boundaries
    mcs_vals = np.arange(int(gt["mcs"].min()), int(gt["mcs"].max()) + 2)
    for m in mcs_vals:
        ax.axvline(m, color="gray", lw=0.4, ls=":", alpha=0.5)

    ax.set_xlabel("Time  (MCS + fraction through MCS)", fontsize=11)
    ax.set_ylabel("Actual volume (voxels)", fontsize=11)
    ax.set_title(
        f"Within-MCS volume fluctuations  "
        f"({len(unique_sigs)} cells, "
        f"{int(gt['mcs'].min())}–{int(gt['mcs'].max())} MCS)",
        fontsize=12,
    )
    ax.grid(True, alpha=0.25)

    # legend only if few cells
    if len(unique_sigs) <= 12:
        ax.legend(fontsize=7, ncol=2, loc="upper left")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"→ {out_path}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Plot cell growth trajectories and growth speed vs radius")
    ap.add_argument("data_dir",
                    help="Directory containing cells.csv (output of 'cargo run --bin analyze')")
    ap.add_argument("--out", default=None,
                    help="Output directory for figures (default: same as data_dir)")
    ap.add_argument("--min-points", type=int, default=3,
                    help="Min snapshots per track required for speed estimates (default: 3)")
    ap.add_argument("--r-bins", type=int, default=20,
                    help="Number of radius bins for the mean speed curve (default: 20)")
    ap.add_argument("--mcs-range", type=int, nargs=2, default=None,
                    metavar=("LO", "HI"),
                    help="MCS window for within-MCS figure, e.g. --mcs-range 10 50")
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    out_dir  = Path(args.out) if args.out else data_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = data_dir / "cells.csv"
    if not csv_path.exists():
        sys.exit(f"Error: {csv_path} not found — run 'cargo run --bin analyze -- --dir {data_dir}' first")

    df = pd.read_csv(csv_path)
    df = df[(df["dying"] == 0) & (df["volume"] > 0)].copy()
    if df.empty:
        sys.exit("No usable cell rows in cells.csv")

    print(f"Loaded {len(df)} cell-snapshot rows  "
          f"({df['sigma'].nunique()} unique sigmas, "
          f"{df['mcs'].nunique()} snapshots)")

    # ── Figure 1 ─────────────────────────────────────────────────────────────
    plot_trajectories(df, out_dir / "cell_growth.png")

    # ── Figure 2 ─────────────────────────────────────────────────────────────
    speed_df = compute_speeds(df, min_points=args.min_points)
    plot_speed_vs_radius(speed_df, out_dir / "growth_speed_vs_r.png",
                         n_bins=args.r_bins)

    # ── Figure 3: within-MCS growth (optional) ───────────────────────────────
    gt_path = data_dir / "growth_trace.csv"
    if gt_path.exists():
        print(f"Found {gt_path} — loading sub-MCS growth trace …")
        gt = load_growth_trace(gt_path, mcs_range=args.mcs_range)
        plot_within_mcs_growth(gt, out_dir / "within_mcs_growth.png")
    else:
        print(f"No growth_trace.csv in {data_dir} — skipping within-MCS figure "
              f"(re-run simulate_life with --track-growth)")


if __name__ == "__main__":
    main()
