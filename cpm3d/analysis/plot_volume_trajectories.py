#!/usr/bin/env python3
"""
plot_volume_trajectories.py
---------------------------
Plot per-cell volume trajectories and population mean from growth_trace.csv,
plus a growth-rate-vs-radius figure.

The CSV is produced by simulate_life with --track-growth.
Column meaning after the user's sub-MCS timing change:

    kind            "grow" | "shrink"
    mcs             integer MCS index
    sigma           cell id
    area_at_event   cell volume immediately after the accepted flip
    birth_mcs       attempt index t within this MCS  (0 … mcs_size-1)
    lifetime_mcs    same as mcs (redundant, kept for compatibility)

Continuous time axis:
    time = birth_mcs / (w * h * d) + mcs

w*h*d (= mcs_size) is read from the first state JSON in the same directory.
You can override it with --mcs-size.

Usage
-----
  python analysis/plot_volume_trajectories.py <data_dir>
  python analysis/plot_volume_trajectories.py <data_dir> --out figures/
  python analysis/plot_volume_trajectories.py <data_dir> --mcs-range 0 50
  python analysis/plot_volume_trajectories.py <data_dir> --mcs-size 8000
  python analysis/plot_volume_trajectories.py <data_dir> --window 0.5 --r-bins 25
"""

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# ── Helpers ───────────────────────────────────────────────────────────────────

def read_mcs_size(data_dir: Path) -> int:
    """Try to read w*h*d from the earliest state_mcs*.json in data_dir."""
    jsons = sorted(data_dir.glob("state_mcs*.json"))
    if not jsons:
        return None
    with open(jsons[0]) as f:
        s = json.load(f)
    p = s["params"]
    if p.get("mcs_per_step") is not None:
        return int(p["mcs_per_step"])
    return int(p["grid_w"]) * int(p["grid_h"]) * int(p["grid_d"])


def load_trace(path: Path, mcs_size: int, mcs_range=None) -> pd.DataFrame:
    """
    Load growth_trace.csv and compute a continuous time axis:
        time = birth_mcs / mcs_size + mcs
    """
    df = pd.read_csv(path)
    if mcs_range is not None:
        lo, hi = mcs_range
        df = df[(df["mcs"] >= lo) & (df["mcs"] <= hi)]
    if df.empty:
        return df
    df["time"] = df["birth_mcs"] / mcs_size + df["mcs"]
    return df


# ── Main figure ───────────────────────────────────────────────────────────────

def plot_trajectories(df: pd.DataFrame, mcs_size: int, out_path: Path):
    """
    One subplot: individual cell volume traces (thin, transparent) +
    population mean (thick black) ± 1 SD (shaded).

    The mean is computed by resampling every cell to a shared time grid
    and averaging (forward-fill between events).
    """
    # Use both grow and shrink events: area_at_event is the volume *after* the flip
    sigmas = sorted(df["sigma"].unique())
    n      = len(sigmas)

    cmap      = matplotlib.colormaps.get_cmap("tab20").resampled(max(n, 1))
    sig_color = {s: cmap(i % 20) for i, s in enumerate(sigmas)}

    # Build per-cell time-sorted series ──────────────────────────────────────
    tracks = {}
    for sigma, grp in df.groupby("sigma"):
        grp = grp.sort_values("time")
        tracks[sigma] = (grp["time"].to_numpy(), grp["area_at_event"].to_numpy())

    # Common fine time grid for mean ─────────────────────────────────────────
    t_min = df["time"].min()
    t_max = df["time"].max()
    # ~4 points per MCS is enough for a smooth mean
    n_grid   = max(int((t_max - t_min) * 4), 200)
    t_grid   = np.linspace(t_min, t_max, n_grid)

    # Resample each track onto t_grid with forward-fill
    mat = np.full((n, n_grid), np.nan)
    for row, sigma in enumerate(sigmas):
        t, v = tracks[sigma]
        if len(t) == 0:
            continue
        # np.searchsorted gives insert position; we want the last event ≤ grid pt
        idx = np.searchsorted(t, t_grid, side="right") - 1
        valid = idx >= 0
        mat[row, valid] = v[idx[valid]]

    mean_v = np.nanmean(mat, axis=0)
    std_v  = np.nanstd(mat, axis=0)
    n_obs  = np.sum(~np.isnan(mat), axis=0)

    # ── Figure ───────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(13, 6))

    # Individual traces
    for sigma in sigmas:
        t, v = tracks[sigma]
        ax.plot(t, v,
                color=sig_color[sigma], lw=0.6, alpha=0.35,
                drawstyle="steps-post")

    # Population mean ± SD
    valid = n_obs > 0
    ax.plot(t_grid[valid], mean_v[valid],
            color="black", lw=2.2, zorder=6, label=f"mean (n={n})")
    ax.fill_between(t_grid[valid],
                    mean_v[valid] - std_v[valid],
                    mean_v[valid] + std_v[valid],
                    color="black", alpha=0.15, zorder=5, label="±1 SD")

    # MCS boundary lines
    mcs_lo = int(df["mcs"].min())
    mcs_hi = int(df["mcs"].max()) + 1
    for m in range(mcs_lo, mcs_hi + 1):
        ax.axvline(m, color="steelblue", lw=0.35, ls=":", alpha=0.4)

    # Legend: only label individual cells if few enough
    if n <= 16:
        handles = [
            plt.Line2D([0], [0], color=sig_color[s], lw=1.5, label=f"σ={s}")
            for s in sigmas
        ]
        handles += ax.get_legend_handles_labels()[0]  # mean + SD
        ax.legend(handles=handles, fontsize=7, ncol=2, loc="upper left")
    else:
        ax.legend(fontsize=9, loc="upper left")

    ax.set_xlabel(
        r"Time  $\left(t_{\mathrm{attempt}} \,/\, WHD + \mathrm{MCS}\right)$",
        fontsize=12,
    )
    ax.set_ylabel("Volume (voxels)", fontsize=12)
    ax.set_title(
        f"Cell volume trajectories  "
        f"({n} cells,  {mcs_lo}–{mcs_hi-1} MCS,  "
        f"mcs_size={mcs_size})",
        fontsize=12,
    )
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"→ {out_path}")


# ── Growth rate vs radius ─────────────────────────────────────────────────────

def vol_to_radius(v: np.ndarray) -> np.ndarray:
    """Spherical approximation: r = (3V / 4π)^(1/3)."""
    return np.cbrt(3.0 * np.asarray(v, dtype=float) / (4.0 * np.pi))


def compute_growth_rates(
    df: pd.DataFrame,
    window: float = 0.5,
    min_events: int = 4,
) -> pd.DataFrame:
    """
    Sliding-window linear regression of volume vs time for each cell.

    For every window of width `window` MCS units (step = window/2) that
    contains at least `min_events` events, fit a line to (time, volume) and
    record:
        sigma        cell id
        t_center     midpoint of the window (MCS units)
        mean_volume  average volume in the window
        radius       spherical radius from mean_volume
        growth_rate  slope dV/dt  (voxels / MCS)

    Parameters
    ----------
    window : float
        Width of the sliding window in MCS units.
    min_events : int
        Minimum number of events required inside a window.
    """
    step = window / 2.0
    records = []

    for sigma, grp in df.groupby("sigma"):
        grp   = grp.sort_values("time")
        times = grp["time"].to_numpy(dtype=float)
        vols  = grp["area_at_event"].to_numpy(dtype=float)

        t_lo = times[0]
        t_hi_max = times[-1] - window
        if t_hi_max < t_lo:
            continue  # track shorter than one window

        t_start = t_lo
        while t_start <= t_hi_max:
            t_end = t_start + window
            mask  = (times >= t_start) & (times < t_end)
            n_ev  = mask.sum()
            if n_ev >= min_events:
                t_win = times[mask]
                v_win = vols[mask]
                # Weighted linear regression: more weight to later events
                # (they represent the "current" state better in noisy CPM)
                slope, _ = np.polyfit(t_win, v_win, 1)
                mean_v   = v_win.mean()
                records.append({
                    "sigma":       sigma,
                    "t_center":    (t_start + t_end) / 2.0,
                    "mean_volume": mean_v,
                    "radius":      float(vol_to_radius(mean_v)),
                    "growth_rate": slope,
                })
            t_start += step

    return pd.DataFrame(records) if records else pd.DataFrame(
        columns=["sigma", "t_center", "mean_volume", "radius", "growth_rate"]
    )


def plot_growth_rate_vs_radius(
    gr: pd.DataFrame,
    out_path: Path,
    window: float,
    n_bins: int = 20,
):
    """
    Scatter of (radius, dV/dt) coloured by cell, with binned mean ± 1 SD.

    Each point is one sliding-window estimate from `compute_growth_rates`.
    The binned mean smooths out fluctuations and reveals the population trend.
    """
    if gr.empty:
        print("No growth-rate estimates — skipping growth_rate_vs_radius.png")
        return

    sigmas      = sorted(gr["sigma"].unique())
    n_cells     = len(sigmas)
    cmap        = matplotlib.colormaps.get_cmap("tab20").resampled(max(n_cells, 1))
    sig_color   = {s: cmap(i % 20) for i, s in enumerate(sigmas)}

    r    = gr["radius"].to_numpy()
    rate = gr["growth_rate"].to_numpy()

    fig, ax = plt.subplots(figsize=(9, 6))

    # ── Per-cell scatter ──────────────────────────────────────────────────────
    for sigma in sigmas:
        sub = gr[gr["sigma"] == sigma]
        ax.scatter(
            sub["radius"], sub["growth_rate"],
            color=sig_color[sigma], s=10, alpha=0.25,
            linewidths=0, label=f"σ={sigma}",
            zorder=2,
        )

    # ── Binned mean ± SD ─────────────────────────────────────────────────────
    r_min, r_max = r.min(), r.max()
    edges   = np.linspace(r_min, r_max, n_bins + 1)
    centres, means, stds, counts = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (r >= lo) & (r < hi)
        sub  = rate[mask]
        if len(sub) < 3:
            continue
        centres.append((lo + hi) / 2.0)
        means.append(sub.mean())
        stds.append(sub.std())
        counts.append(len(sub))

    if centres:
        c_arr = np.array(centres)
        m_arr = np.array(means)
        s_arr = np.array(stds)
        ax.plot(c_arr, m_arr, color="crimson", lw=2.5, zorder=6,
                label="binned mean")
        ax.fill_between(c_arr, m_arr - s_arr, m_arr + s_arr,
                        color="crimson", alpha=0.20, zorder=5, label="±1 SD")

    ax.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)

    ax.set_xlabel(
        r"Radius  $r = \left(\frac{3\,\langle V \rangle}{4\pi}\right)^{1/3}$"
        "  (voxels)",
        fontsize=12,
    )
    ax.set_ylabel(r"Growth rate  $\mathrm{d}V/\mathrm{d}t$  (voxels / MCS)",
                  fontsize=12)
    ax.set_title(
        f"Growth rate vs cell radius  "
        f"(window={window} MCS,  {len(gr)} estimates,  {n_cells} cells)",
        fontsize=12,
    )

    if n_cells <= 12:
        ax.legend(fontsize=7, ncol=2, loc="upper left")
    else:
        # Only show the binned-mean legend entry
        handles = [
            plt.Line2D([0], [0], color="crimson", lw=2, label="binned mean"),
            plt.Patch(color="crimson", alpha=0.25, label="±1 SD"),
        ]
        ax.legend(handles=handles, fontsize=9, loc="upper left")

    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"→ {out_path}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir",
                    help="Directory containing growth_trace.csv (and state JSON files)")
    ap.add_argument("--out", default=None,
                    help="Output directory for figures (default: data_dir)")
    ap.add_argument("--mcs-range", type=int, nargs=2, default=None,
                    metavar=("LO", "HI"),
                    help="Restrict to this MCS window, e.g. --mcs-range 0 100")
    ap.add_argument("--mcs-size", type=int, default=None,
                    help="Override w*h*d (attempts per MCS).  "
                         "Auto-detected from state JSON if omitted.")
    ap.add_argument("--window", type=float, default=0.5,
                    help="Sliding-window width for growth-rate estimation, "
                         "in MCS units (default: 0.5).  "
                         "Expressed as a fraction of the time axis t/WHD+MCS.")
    ap.add_argument("--r-bins", type=int, default=20,
                    help="Number of radius bins for the binned mean curve "
                         "(default: 20).")
    ap.add_argument("--min-events", type=int, default=4,
                    help="Minimum events per window for a growth-rate estimate "
                         "(default: 4).")
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    out_dir  = Path(args.out) if args.out else data_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Locate growth_trace.csv ───────────────────────────────────────────────
    gt_path = data_dir / "growth_trace.csv"
    if not gt_path.exists():
        sys.exit(
            f"Error: {gt_path} not found.\n"
            f"Re-run simulate_life with --track-growth to produce it."
        )

    # ── Grid size ─────────────────────────────────────────────────────────────
    mcs_size = args.mcs_size
    if mcs_size is None:
        mcs_size = read_mcs_size(data_dir)
    if mcs_size is None:
        sys.exit(
            "Cannot determine mcs_size: no state_mcs*.json found in data_dir.\n"
            "Pass --mcs-size W*H*D explicitly."
        )
    print(f"mcs_size (w*h*d) = {mcs_size}")

    # ── Load ─────────────────────────────────────────────────────────────────
    df = load_trace(gt_path, mcs_size, mcs_range=args.mcs_range)
    if df.empty:
        sys.exit("growth_trace.csv is empty (or MCS range matches no rows).")

    print(
        f"Loaded {len(df):,} events  "
        f"({df['sigma'].nunique()} cells, "
        f"MCS {int(df['mcs'].min())}–{int(df['mcs'].max())})"
    )

    # ── Figure 1: volume trajectories ────────────────────────────────────────
    plot_trajectories(df, mcs_size, out_dir / "volume_trajectories.png")

    # ── Figure 2: growth rate vs radius ──────────────────────────────────────
    print(f"Computing growth rates  (window={args.window} MCS, "
          f"min_events={args.min_events}) …")
    gr = compute_growth_rates(df, window=args.window, min_events=args.min_events)
    if gr.empty:
        print("  No windows had enough events — try a larger --window or "
              "a longer simulation with --track-growth.")
    else:
        print(f"  {len(gr)} window estimates across {gr['sigma'].nunique()} cells")
        plot_growth_rate_vs_radius(
            gr, out_dir / "growth_rate_vs_radius.png",
            window=args.window, n_bins=args.r_bins,
        )


if __name__ == "__main__":
    main()
