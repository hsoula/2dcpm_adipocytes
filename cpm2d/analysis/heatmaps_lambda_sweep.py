"""
heatmaps_lambda_sweep.py
========================
Visualise the lambda parameter sweep produced by `analyze_sweep`.

For every unique value of lambda_iso the script produces one figure with
six heatmaps (lambda_area × lambda_perimeter).  The colour encodes a
different shape metric in each panel, making it easy to identify the
parameter regime that gives realistic adipocyte morphology.

Usage
-----
    python analysis/heatmaps_lambda_sweep.py                       # default path
    python analysis/heatmaps_lambda_sweep.py --csv results/sweep.csv --out figs/

Metrics shown
-------------
  area_ratio        mean(area) / target_area          → 1.0 = on target
  iso_ratio         mean(p² / 4πA)                    → 1.0 = circle
  solidity          mean(area / convex_hull_area)      → 1.0 = convex
  profile_convexity (flat_edge + outer_corner) / total_boundary
  filament_frac     (filament_tip + filament_neck) / total_boundary  (want ~0)
  concave_frac      (inner_corner + concave_bay)  / total_boundary
"""

import argparse
import os
import warnings

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns

# ── CLI ──────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--csv", default="../results/sweep.csv",
                    help="Path to the analyze_sweep CSV output")
parser.add_argument("--out", default="../figs/heatmaps",
                    help="Output directory for saved figures")
parser.add_argument("--fmt", default="pdf", choices=["pdf", "png", "svg"],
                    help="Figure file format")
parser.add_argument("--last-mcs-only", action="store_true",
                    help="Keep only the final MCS snapshot per run (default: all)")
parser.add_argument("--annot", action="store_true",
                    help="Print numeric values inside each heatmap cell")
args = parser.parse_args()

os.makedirs(args.out, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
print(f"Loading {args.csv} …")
df = pd.read_csv(args.csv)
print(f"  {len(df):,} rows  |  columns: {list(df.columns)}")

# Drop medium cell (cell_id == 0)
df = df[df["cell_id"] > 0].copy()

# Optionally restrict to last MCS per (lambda combo, run)
if args.last_mcs_only:
    idx = df.groupby(["lambda_area", "lambda_perimeter", "lambda_iso", "run"])["mcs"].transform("max")
    df = df[df["mcs"] == idx].copy()
    print(f"  After last-MCS filter: {len(df):,} rows")

# ── Derived metrics ───────────────────────────────────────────────────────────
PI = np.pi

df["iso_ratio"]      = (df["perimeter"] ** 2) / (4 * PI * df["area"].clip(lower=1))
df["area_ratio"]     = df["area"] / df["target_area"].clip(lower=1)
df["filament_frac"]  = (df["filament_tip"] + df["filament_neck"]) / df["total_boundary"].clip(lower=1)
df["concave_frac"]   = (df["inner_corner"] + df["concave_bay"])   / df["total_boundary"].clip(lower=1)

# ── Metrics to plot ───────────────────────────────────────────────────────────
METRICS = [
    # (column,           label,                    cmap,        vmin, vmax,  fmt)
    ("area_ratio",       "Area ratio\n(mean/target)",    "RdYlGn",    0.5,  1.5,  ".2f"),
    ("iso_ratio",        "Isoperimetric ratio\n(p²/4πA)","RdYlGn_r",  1.0,  3.0,  ".2f"),
    ("solidity",         "Solidity\n(area/hull area)",   "YlGnBu",    0.4,  1.0,  ".2f"),
    ("profile_convexity","Profile convexity\n(flat+outer/total)", "YlGnBu", 0.0, 1.0, ".2f"),
    ("filament_frac",    "Filament fraction\n(tip+neck/total)",   "YlOrRd",  0.0, 0.4, ".2f"),
    ("concave_frac",     "Concave fraction\n(inner+bay/total)",   "YlOrRd",  0.0, 0.6, ".2f"),
]

GROUP_COLS = ["lambda_area", "lambda_perimeter", "lambda_iso"]
agg = df.groupby(GROUP_COLS)[[m[0] for m in METRICS]].mean().reset_index()

iso_vals = sorted(agg["lambda_iso"].unique())
la_vals  = sorted(agg["lambda_area"].unique())
lp_vals  = sorted(agg["lambda_perimeter"].unique())

print(f"\n  lambda_area   : {la_vals}")
print(f"  lambda_perim  : {lp_vals}")
print(f"  lambda_iso    : {iso_vals}")
print(f"\nGenerating {len(iso_vals)} figure(s) …")

# ── Plot ──────────────────────────────────────────────────────────────────────
NCOLS = 3
NROWS = int(np.ceil(len(METRICS) / NCOLS))

for iso in iso_vals:
    sub = agg[agg["lambda_iso"] == iso]

    fig, axes = plt.subplots(NROWS, NCOLS,
                             figsize=(5 * NCOLS, 4.5 * NROWS),
                             constrained_layout=True)
    axes = np.array(axes).reshape(NROWS, NCOLS)

    fig.suptitle(f"Lambda sweep — $\\lambda_{{iso}}$ = {iso:.2f}",
                 fontsize=15, fontweight="bold")

    for ax_idx, (col, label, cmap, vmin, vmax, fmt) in enumerate(METRICS):
        row, col_pos = divmod(ax_idx, NCOLS)
        ax = axes[row, col_pos]

        pivot = sub.pivot(index="lambda_perimeter",
                          columns="lambda_area",
                          values=col)
        # Ensure axes are in ascending order
        pivot = pivot.sort_index(ascending=False)   # perimeter on y, high at top
        pivot.columns = pivot.columns.astype(float)
        pivot.index   = pivot.index.astype(float)

        sns.heatmap(
            pivot,
            ax=ax,
            cmap=cmap,
            vmin=vmin, vmax=vmax,
            annot=args.annot,
            fmt=fmt if args.annot else "",
            linewidths=0.4,
            linecolor="white",
            cbar_kws={"shrink": 0.8, "format": mticker.FormatStrFormatter(f"%{fmt}")},
        )

        ax.set_title(label, fontsize=11)
        ax.set_xlabel("$\\lambda_{area}$", fontsize=10)
        ax.set_ylabel("$\\lambda_{perim}$", fontsize=10)

        # Format tick labels to 2 d.p.
        ax.set_xticklabels([f"{float(t.get_text()):.2f}" for t in ax.get_xticklabels()],
                           rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels([f"{float(t.get_text()):.2f}" for t in ax.get_yticklabels()],
                           rotation=0, fontsize=8)

    # Hide any unused axes
    for ax_idx in range(len(METRICS), NROWS * NCOLS):
        row, col_pos = divmod(ax_idx, NCOLS)
        axes[row, col_pos].set_visible(False)

    out_path = os.path.join(args.out, f"heatmap_iso{iso:.2f}.{args.fmt}")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved → {out_path}")

# ── Summary scatter: good-regime candidates ───────────────────────────────────
# A parameter set is "candidate" if:
#   area_ratio ∈ [0.8, 1.2]   (area close to target)
#   iso_ratio  ∈ [1.0, 1.5]   (roughly circular, not spiky)
#   solidity   ≥ 0.7           (compact shape)
#   filament_frac ≤ 0.1        (few pathological pixels)

good = agg[
    agg["area_ratio"].between(0.8, 1.2) &
    agg["iso_ratio"].between(1.0, 1.5)  &
    (agg["solidity"] >= 0.7)            &
    (agg["filament_frac"] <= 0.1)
].copy()

if good.empty:
    print("\nNo parameter combinations passed all quality filters.")
else:
    print(f"\n{len(good)} candidate parameter combination(s):")
    print(good[["lambda_area", "lambda_perimeter", "lambda_iso",
                "area_ratio", "iso_ratio", "solidity", "filament_frac"]]
          .sort_values(["lambda_iso", "lambda_area", "lambda_perimeter"])
          .to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # Save candidates to CSV
    cand_path = os.path.join(args.out, "good_params.csv")
    good.to_csv(cand_path, index=False)
    print(f"\nCandidates saved → {cand_path}")

print("\nDone.")
