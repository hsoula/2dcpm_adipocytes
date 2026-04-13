"""
Plot demography statistics from simulate_live output.

Expected inputs (produced by extract_demography):
  results_dir/population.csv   — mcs, n_living, n_dying, n_dead, ...
  results_dir/events.csv       — kind, mcs, sigma, area_at_event, birth_mcs, lifetime_mcs

Usage:
  python analysis/plot_demography.py --results results/live/ --out figures/demography.pdf
"""

import argparse
import pandas as pd
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
deaths = events[events["kind"] == "death"].copy()

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

# ── Save ──────────────────────────────────────────────────────────────────────
fig.suptitle("CPM demography summary", fontsize=13, y=1.01)
fig.savefig(out_path, bbox_inches="tight", dpi=150)
print(f"Saved → {out_path}")
plt.show()
