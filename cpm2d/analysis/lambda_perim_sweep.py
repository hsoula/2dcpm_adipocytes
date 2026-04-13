"""
lambda_perim_sweep.py
=====================
Focused analysis following the observation that:
  - lambda_iso  has almost no effect → collapsed (averaged over)
  - lambda_area has a binary effect  → split into "off" (=0) vs "on" (>0)
  - lambda_perimeter is the dominant axis

Produces:
  1. Line plots of every shape metric vs lambda_perimeter,
     with two traces: lambda_area = 0 and lambda_area > 0.
     Shaded band = ±1 std across runs / cells.

  2. Confirmation panels showing that lambda_iso and lambda_area
     (within the non-zero range) genuinely don't matter:
       - metric vs lambda_iso  (one line per lambda_perimeter)
       - metric vs lambda_area (one line per lambda_perimeter)

Usage
-----
    python analysis/lambda_perim_sweep.py
    python analysis/lambda_perim_sweep.py --csv results/sweep.csv --out figs/
"""

import argparse
import os

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

# ── CLI ──────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--csv", default="../results/sweep.csv")
parser.add_argument("--out", default="../figs/perim_sweep")
parser.add_argument("--fmt", default="pdf", choices=["pdf", "png", "svg"])
parser.add_argument("--last-mcs-only", action="store_true",
                    help="Restrict to final MCS snapshot per run")
args = parser.parse_args()

os.makedirs(args.out, exist_ok=True)

# ── Load ──────────────────────────────────────────────────────────────────────
print(f"Loading {args.csv} …")
df = pd.read_csv(args.csv)
df = df[df["cell_id"] > 0].copy()

if args.last_mcs_only:
    idx = df.groupby(["lambda_area", "lambda_perimeter", "lambda_iso", "run"])["mcs"].transform("max")
    df  = df[df["mcs"] == idx].copy()

PI = np.pi
df["iso_ratio"]     = df["perimeter"]**2 / (4 * PI * df["area"].clip(lower=1))
df["area_ratio"]    = df["area"] / df["target_area"].clip(lower=1)
df["filament_frac"] = (df["filament_tip"]  + df["filament_neck"]) / df["total_boundary"].clip(lower=1)
df["concave_frac"]  = (df["inner_corner"]  + df["concave_bay"])   / df["total_boundary"].clip(lower=1)

# Binary lambda_area label
df["la_group"] = df["lambda_area"].apply(lambda x: "off (λ_a = 0)" if x == 0.0 else "on  (λ_a > 0)")

METRICS = [
    ("area_ratio",        "Area ratio (mean / target)",         (0.5, 1.5)),
    ("iso_ratio",         "Isoperimetric ratio (p² / 4πA)",     (1.0, 4.0)),
    ("solidity",          "Solidity (area / hull area)",         (0.3, 1.05)),
    ("profile_convexity", "Profile convexity",                   (0.0, 1.0)),
    ("filament_frac",     "Filament fraction (tip + neck)",      (0.0, 0.5)),
    ("concave_frac",      "Concave fraction (inner + bay)",      (0.0, 0.7)),
]

COLORS = {"off (λ_a = 0)": "#e74c3c", "on  (λ_a > 0)": "#2980b9"}
LP_VALS = sorted(df["lambda_perimeter"].unique())

# ══════════════════════════════════════════════════════════════════════════════
# Figure 1 — metrics vs lambda_perimeter, la_group as two traces
# ══════════════════════════════════════════════════════════════════════════════
NCOLS = 2
NROWS = int(np.ceil(len(METRICS) / NCOLS))
fig1, axes1 = plt.subplots(NROWS, NCOLS, figsize=(7 * NCOLS, 4 * NROWS),
                            constrained_layout=True)
axes1 = np.array(axes1).reshape(NROWS, NCOLS)
fig1.suptitle("Metrics vs $\\lambda_{perim}$  (collapsed over $\\lambda_{iso}$, $\\lambda_{area}$ binary)",
              fontsize=13, fontweight="bold")

grp_lp = df.groupby(["lambda_perimeter", "la_group"])

for idx, (col, label, (ylo, yhi)) in enumerate(METRICS):
    row, c = divmod(idx, NCOLS)
    ax = axes1[row, c]

    stats = grp_lp[col].agg(["mean", "std"]).reset_index()

    for group, colour in COLORS.items():
        sub = stats[stats["la_group"] == group].sort_values("lambda_perimeter")
        if sub.empty:
            continue
        ax.plot(sub["lambda_perimeter"], sub["mean"],
                color=colour, marker="o", linewidth=2, label=group)
        ax.fill_between(sub["lambda_perimeter"],
                        sub["mean"] - sub["std"],
                        sub["mean"] + sub["std"],
                        color=colour, alpha=0.15)

    ax.set_xlabel("$\\lambda_{perim}$", fontsize=11)
    ax.set_ylabel(label, fontsize=10)
    ax.set_ylim(ylo, yhi)
    ax.axhline(1.0 if "ratio" in col or "convexity" in col else ylo,
               color="grey", linewidth=0.8, linestyle="--", alpha=0.6)
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.legend(fontsize=8)

for idx in range(len(METRICS), NROWS * NCOLS):
    row, c = divmod(idx, NCOLS)
    axes1[row, c].set_visible(False)

path1 = os.path.join(args.out, f"metrics_vs_lambda_perim.{args.fmt}")
fig1.savefig(path1, dpi=150)
plt.close(fig1)
print(f"Saved → {path1}")

# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 — confirm lambda_iso has no effect
#   one subplot per metric; x = lambda_iso; one line per lambda_perimeter
# ══════════════════════════════════════════════════════════════════════════════
# Restrict to lambda_area > 0 so we see the "active" regime
df_on = df[df["lambda_area"] > 0]

lp_uniq = sorted(df_on["lambda_perimeter"].unique())
cmap_lp  = plt.cm.get_cmap("viridis", len(lp_uniq))
lp_color = {lp: cmap_lp(i) for i, lp in enumerate(lp_uniq)}

fig2, axes2 = plt.subplots(NROWS, NCOLS, figsize=(7 * NCOLS, 4 * NROWS),
                            constrained_layout=True)
axes2 = np.array(axes2).reshape(NROWS, NCOLS)
fig2.suptitle("Metrics vs $\\lambda_{iso}$  — confirming negligible effect\n"
              "(λ_area > 0, one line per λ_perim)",
              fontsize=13, fontweight="bold")

grp_li = df_on.groupby(["lambda_iso", "lambda_perimeter"])

for idx, (col, label, (ylo, yhi)) in enumerate(METRICS):
    row, c = divmod(idx, NCOLS)
    ax = axes2[row, c]

    stats = grp_li[col].mean().reset_index()

    for lp in lp_uniq:
        sub = stats[stats["lambda_perimeter"] == lp].sort_values("lambda_iso")
        ax.plot(sub["lambda_iso"], sub[col],
                color=lp_color[lp], marker="o", linewidth=1.8, markersize=4)

    ax.set_xlabel("$\\lambda_{iso}$", fontsize=11)
    ax.set_ylabel(label, fontsize=10)
    ax.set_ylim(ylo, yhi)
    ax.grid(axis="y", linestyle=":", alpha=0.4)

# Shared legend for lambda_perimeter values
handles = [mlines.Line2D([], [], color=lp_color[lp], marker="o",
                         label=f"λ_p = {lp:.2f}") for lp in lp_uniq]
fig2.legend(handles=handles, title="$\\lambda_{perim}$",
            loc="lower right", ncol=min(4, len(lp_uniq)), fontsize=8)

for idx in range(len(METRICS), NROWS * NCOLS):
    row, c = divmod(idx, NCOLS)
    axes2[row, c].set_visible(False)

path2 = os.path.join(args.out, f"metrics_vs_lambda_iso.{args.fmt}")
fig2.savefig(path2, dpi=150)
plt.close(fig2)
print(f"Saved → {path2}")

# ══════════════════════════════════════════════════════════════════════════════
# Figure 3 — confirm lambda_area has only binary effect
#   x = lambda_area; one line per lambda_perimeter
# ══════════════════════════════════════════════════════════════════════════════
fig3, axes3 = plt.subplots(NROWS, NCOLS, figsize=(7 * NCOLS, 4 * NROWS),
                            constrained_layout=True)
axes3 = np.array(axes3).reshape(NROWS, NCOLS)
fig3.suptitle("Metrics vs $\\lambda_{area}$  — confirming binary / scale-invariant effect\n"
              "(one line per λ_perim, averaged over λ_iso)",
              fontsize=13, fontweight="bold")

grp_la = df.groupby(["lambda_area", "lambda_perimeter"])

for idx, (col, label, (ylo, yhi)) in enumerate(METRICS):
    row, c = divmod(idx, NCOLS)
    ax = axes3[row, c]

    stats = grp_la[col].mean().reset_index()

    for lp in lp_uniq:
        sub = stats[stats["lambda_perimeter"] == lp].sort_values("lambda_area")
        ax.plot(sub["lambda_area"], sub[col],
                color=lp_color[lp], marker="o", linewidth=1.8, markersize=4)

    # Mark the 0→non-zero transition
    ax.axvline(0, color="grey", linewidth=1.0, linestyle="--", alpha=0.7)
    ax.set_xlabel("$\\lambda_{area}$", fontsize=11)
    ax.set_ylabel(label, fontsize=10)
    ax.set_ylim(ylo, yhi)
    ax.grid(axis="y", linestyle=":", alpha=0.4)

fig3.legend(handles=handles, title="$\\lambda_{perim}$",
            loc="lower right", ncol=min(4, len(lp_uniq)), fontsize=8)

for idx in range(len(METRICS), NROWS * NCOLS):
    row, c = divmod(idx, NCOLS)
    axes3[row, c].set_visible(False)

path3 = os.path.join(args.out, f"metrics_vs_lambda_area.{args.fmt}")
fig3.savefig(path3, dpi=150)
plt.close(fig3)
print(f"Saved → {path3}")

# ══════════════════════════════════════════════════════════════════════════════
# Summary table — best lambda_perimeter (lambda_area > 0, any lambda_iso)
# ══════════════════════════════════════════════════════════════════════════════
summary = (df_on
           .groupby("lambda_perimeter")[
               ["area_ratio", "iso_ratio", "solidity",
                "profile_convexity", "filament_frac", "concave_frac"]
           ]
           .agg(["mean", "std"])
           .round(3))
summary.columns = ["_".join(c) for c in summary.columns]

print("\n── Summary by lambda_perimeter (λ_area > 0, λ_iso collapsed) ──")
print(summary.to_string())

summary.to_csv(os.path.join(args.out, "summary_by_lambda_perim.csv"))
print(f"\nSaved → {os.path.join(args.out, 'summary_by_lambda_perim.csv')}")
print("\nDone.")
