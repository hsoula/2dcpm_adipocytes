"""
Plot demography statistics from simulate_live output.

Expected inputs (produced by extract_demography):
  results_dir/population.csv   — mcs, n_living, n_dying, n_dead, ...
  results_dir/events.csv       — kind, mcs, sigma, area_at_event, birth_mcs, lifetime_mcs
  results_dir/cell_areas.csv   — mcs, sigma, area, target_area, status, birth_mcs

Outputs:
  <out>                        — population summary figure (3 panels)
  <out stem>_areas.<ext>       — area distribution + cell trajectories (3 panels)

Usage:
  python analysis/plot_demography.py --results results/live/ --out figures/demography.pdf
"""

import argparse
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from pathlib import Path

# ── CLI ───────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument("--results", default="results/live/",
                    help="Directory containing population.csv and events.csv")
parser.add_argument("--out", default="figures/demography.pdf",
                    help="Output figure path (.pdf or .png)")
args = parser.parse_args()

results_dir = Path(args.results)
out_path    = Path(args.out)
out_path.parent.mkdir(parents=True, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────

pop = pd.read_csv(results_dir / "population.csv")
events = pd.read_csv(results_dir / "events.csv")

births = events[events["kind"] == "birth"].copy()
deaths = events[events["kind"] == "dying"].copy()

# ── Figure ────────────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(14, 10))
gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)

# ── Panel 1: population over time ────────────────────────────────────────────
ax1 = fig.add_subplot(gs[0, :])  # full top row

ax1.stackplot(pop["mcs"],
              pop["n_dying"],
              pop["n_living"],
              labels=["dying", "living"],
              colors=["#f28b30", "#4e9fd4"],
              alpha=0.8)

# Overlay birth/death event rug
if len(births):
    ax1.plot(births["mcs"], np.zeros(len(births)) - 0.3,
             "|", color="#2ca02c", ms=6, alpha=0.6, label="birth events")
if len(deaths):
    ax1.plot(deaths["mcs"], np.zeros(len(deaths)) - 0.6,
             "|", color="#d62728", ms=6, alpha=0.6, label="death events")

# Mean area on secondary axis
ax1b = ax1.twinx()
ax1b.plot(pop["mcs"], pop["mean_area_living"],
          color="#9467bd", lw=1.5, ls="--", label="mean area (living)")
ax1b.plot(pop["mcs"], pop["mean_target_area_living"],
          color="#c5b0d5", lw=1, ls=":", label="mean target area (living)")
ax1b.set_ylabel("Area (pixels)", color="#9467bd")
ax1b.tick_params(axis="y", labelcolor="#9467bd")

ax1.set_xlabel("MCS")
ax1.set_ylabel("Cell count")
ax1.set_title("Population dynamics")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax1b.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left", fontsize=8)

# ── Panel 2: cell lifetime distribution ──────────────────────────────────────
ax2 = fig.add_subplot(gs[1, 0])

if len(deaths) > 0:
    bins = min(40, max(10, len(deaths) // 5))
    ax2.hist(deaths["lifetime_mcs"], bins=bins,
             color="#d62728", alpha=0.75, edgecolor="white")
    median_lt = deaths["lifetime_mcs"].median()
    ax2.axvline(median_lt, color="black", ls="--", lw=1.2,
                label=f"median = {median_lt:.0f} MCS")
    ax2.legend(fontsize=8)
else:
    ax2.text(0.5, 0.5, "no deaths recorded", ha="center", va="center",
             transform=ax2.transAxes)

ax2.set_xlabel("Lifetime (MCS)")
ax2.set_ylabel("Count")
ax2.set_title("Cell lifetime distribution")

# ── Panel 3: area at death vs birth MCS ──────────────────────────────────────
ax3 = fig.add_subplot(gs[1, 1])

if len(deaths) > 0:
    sc = ax3.scatter(deaths["birth_mcs"], deaths["area_at_event"],
                     c=deaths["lifetime_mcs"], cmap="viridis",
                     s=30, alpha=0.7, edgecolors="none")
    cbar = fig.colorbar(sc, ax=ax3)
    cbar.set_label("Lifetime (MCS)", fontsize=8)

    # Growth rate estimate: (area_at_death - 4) / lifetime
    deaths_nonzero = deaths[deaths["lifetime_mcs"] > 0].copy()
    if len(deaths_nonzero):
        deaths_nonzero["growth_achieved"] = (
            (deaths_nonzero["area_at_event"] - 4) / deaths_nonzero["lifetime_mcs"]
        )
        mean_gr = deaths_nonzero["growth_achieved"].mean()
        ax3.set_title(f"Area at death  (mean growth ≈ {mean_gr:.3f} px/MCS)")
    else:
        ax3.set_title("Area at death")
else:
    ax3.text(0.5, 0.5, "no deaths recorded", ha="center", va="center",
             transform=ax3.transAxes)
    ax3.set_title("Area at death")

ax3.set_xlabel("Birth MCS")
ax3.set_ylabel("Area at death (pixels)")

# ── Save figure 1 ─────────────────────────────────────────────────────────────
fig.suptitle("CPM demography summary", fontsize=13, y=1.01)
fig.savefig(out_path, bbox_inches="tight", dpi=150)
print(f"Saved → {out_path}")

# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 — area distribution over time + per-cell area trajectories
# ══════════════════════════════════════════════════════════════════════════════

areas_path = results_dir / "cell_areas.csv"
if not areas_path.exists():
    print(f"cell_areas.csv not found at {areas_path} — skipping area figures")
    plt.show()
    exit(0)

ca = pd.read_csv(areas_path)

out_areas = out_path.with_stem(out_path.stem + "_areas") \
            if hasattr(out_path, "with_stem") \
            else out_path.parent / (out_path.stem + "_areas" + out_path.suffix)

fig2 = plt.figure(figsize=(16, 12))
gs2  = gridspec.GridSpec(2, 2, figure=fig2, hspace=0.45, wspace=0.35)

# ── Panel A: per-cell area trajectories ──────────────────────────────────────
axA = fig2.add_subplot(gs2[0, :])   # full top row

# Determine which cells died (appear in events) vs survived
dead_sigmas = set(deaths["sigma"].unique()) if len(deaths) else set()

# Thin line per cell, colour by fate
MAX_TRACES = 200   # cap to avoid overplotting
sigmas = ca["sigma"].unique()
if len(sigmas) > MAX_TRACES:
    rng = np.random.default_rng(42)
    sigmas = rng.choice(sigmas, size=MAX_TRACES, replace=False)

for sigma in sigmas:
    traj = ca[ca["sigma"] == sigma].sort_values("mcs")
    color = "#d62728" if sigma in dead_sigmas else "#4e9fd4"
    axA.plot(traj["mcs"], traj["area"], color=color, lw=0.6, alpha=0.35)

# Mean ± std across all cells at each MCS
agg = ca.groupby("mcs")["area"].agg(["mean", "std"]).reset_index()
axA.plot(agg["mcs"], agg["mean"], color="black", lw=2, label="mean area")
axA.fill_between(agg["mcs"],
                 agg["mean"] - agg["std"],
                 agg["mean"] + agg["std"],
                 color="black", alpha=0.15, label="±1 std")

# Legend proxies for cell colours
from matplotlib.lines import Line2D
proxy_alive = Line2D([0], [0], color="#4e9fd4", lw=1.5, label="survived")
proxy_dead  = Line2D([0], [0], color="#d62728", lw=1.5, label="died")
handles_A, labels_A = axA.get_legend_handles_labels()
axA.legend(handles=[proxy_alive, proxy_dead] + handles_A,
           labels=["survived", "died"] + labels_A,
           fontsize=8, loc="upper left")

axA.set_xlabel("MCS")
axA.set_ylabel("Area (pixels)")
axA.set_title("Per-cell area trajectories  (red = died, blue = survived)")

# ── Panel B: area distribution at evenly-spaced snapshots ────────────────────
axB = fig2.add_subplot(gs2[1, 0])

mcs_vals   = np.sort(ca["mcs"].unique())
n_snapshots = min(8, len(mcs_vals))
snapshot_mcs = mcs_vals[np.linspace(0, len(mcs_vals) - 1, n_snapshots, dtype=int)]

cmap_t = matplotlib.colormaps["plasma"]#plt.get_cmap("plasma", n_snapshots)

# for i, t in enumerate(snapshot_mcs):
#     areas_t = ca[ca["mcs"] == t]["area"]
#     if len(areas_t) < 2:
#         continue
#     areas_t.plot.kde(ax=axB, color=cmap_t(i), lw=1.8,
#                      label=f"MCS {t}")
#     axB.axvline(areas_t.median(), color=cmap_t(i), lw=0.7, ls="--", alpha=0.6)

# axB.set_xlabel("Area (pixels)")
# axB.set_ylabel("Density")
# axB.set_title("Area distribution at selected time points")
# axB.legend(fontsize=7, title="snapshot", title_fontsize=7)
# axB.set_xlim(left=0)

# ── Panel C: final area by fate (living vs dying vs dead) ────────────────────
axC = fig2.add_subplot(gs2[1, 1])

last_per_cell = ca.sort_values("mcs").groupby("sigma").last().reset_index()

# Merge fate from events
last_per_cell["fate"] = last_per_cell["sigma"].apply(
    lambda s: "died" if s in dead_sigmas else "survived"
)

for fate, colour in [("survived", "#4e9fd4"), ("died", "#d62728")]:
    sub = last_per_cell[last_per_cell["fate"] == fate]["area"]
    if len(sub) == 0:
        continue
    bins = np.linspace(0, last_per_cell["area"].max() * 1.05, 30)
    axC.hist(sub, bins=bins, color=colour, alpha=0.65,
             edgecolor="white", label=f"{fate}  (n={len(sub)})")

axC.set_xlabel("Final recorded area (pixels)")
axC.set_ylabel("Count")
axC.set_title("Final area distribution by fate")
axC.legend(fontsize=8)

# ── Save figure 2 ─────────────────────────────────────────────────────────────
fig2.suptitle("CPM cell area dynamics", fontsize=13, y=1.01)
fig2.savefig(out_areas, bbox_inches="tight", dpi=150)
print(f"Saved → {out_areas}")

plt.show()
