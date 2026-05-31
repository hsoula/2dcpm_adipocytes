#!/usr/bin/env python3
"""
plot_volume_trajectories.py
---------------------------
Plot per-cell volume trajectories and population mean from growth_trace.csv.

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

    # ── Plot ─────────────────────────────────────────────────────────────────
    plot_trajectories(df, mcs_size, out_dir / "volume_trajectories.png")


if __name__ == "__main__":
    main()
