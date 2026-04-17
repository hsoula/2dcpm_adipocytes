#!/usr/bin/env python3
"""
plot_cell_growth.py — per-cell volume growth tracking for cpm3d
----------------------------------------------------------------
Reads the JSON state files produced by `simulate` (state_mcs*.json)
and plots per-cell volume trajectories over time.

Usage:
    python analysis/plot_cell_growth.py data/sim3d/ [--out figure.png]
"""

import argparse
import json
import glob
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


# ── I/O ──────────────────────────────────────────────────────────────────────

def load_states(data_dir: str) -> list[dict]:
    """Load all state_mcs*.json files from data_dir, sorted by MCS."""
    pattern = os.path.join(data_dir.rstrip("/"), "state_mcs*.json")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No state_mcs*.json files found in {data_dir!r}")
    records = []
    for f in files:
        with open(f) as fh:
            records.append(json.load(fh))
    return sorted(records, key=lambda r: r["mcs"])


# ── Trajectory extraction ─────────────────────────────────────────────────────

def extract_trajectories(records: list[dict]):
    """
    Returns
    -------
    mcs_per_cell  : dict[int, list[int]]   cell_id -> MCS values
    vol_per_cell  : dict[int, list[int]]   cell_id -> actual volumes
    tvol_per_cell : dict[int, list[int]]   cell_id -> target volumes
    """
    mcs_per_cell  = defaultdict(list)
    vol_per_cell  = defaultdict(list)
    tvol_per_cell = defaultdict(list)

    for rec in records:
        mcs = rec["mcs"]
        for cs in rec.get("cells", []):
            cid = cs["id"]
            if cid == 0 or not cs.get("alive", True):
                continue
            mcs_per_cell[cid].append(mcs)
            vol_per_cell[cid].append(cs["volume"])
            tvol_per_cell[cid].append(cs["target_volume"])

    return mcs_per_cell, vol_per_cell, tvol_per_cell


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_growth(
    mcs_per_cell: dict,
    vol_per_cell: dict,
    tvol_per_cell: dict,
    out_path: str,
):
    n_cells = len(vol_per_cell)
    if n_cells == 0:
        print("No live cells found.")
        return

    fig, (ax_vol, ax_tvol) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    cmap = cm.get_cmap("tab20", max(n_cells, 1))

    for i, cid in enumerate(sorted(vol_per_cell.keys())):
        color = cmap(i / max(n_cells - 1, 1))
        ax_vol.plot(mcs_per_cell[cid], vol_per_cell[cid],
                    color=color, lw=1.2, label=f"cell {cid}")
        ax_tvol.plot(mcs_per_cell[cid], tvol_per_cell[cid],
                     color=color, lw=1.0, ls="--")

    # Ensemble mean ± std
    all_mcs = sorted({m for v in mcs_per_cell.values() for m in v})
    vol_at   = defaultdict(list)
    tvol_at  = defaultdict(list)
    for cid in vol_per_cell:
        for t, v, tv in zip(mcs_per_cell[cid], vol_per_cell[cid], tvol_per_cell[cid]):
            vol_at[t].append(v)
            tvol_at[t].append(tv)

    ts      = [t for t in all_mcs if vol_at[t]]
    mean_v  = np.array([np.mean(vol_at[t])  for t in ts])
    std_v   = np.array([np.std(vol_at[t])   for t in ts])
    mean_tv = np.array([np.mean(tvol_at[t]) for t in ts])
    std_tv  = np.array([np.std(tvol_at[t])  for t in ts])

    ax_vol.plot(ts, mean_v, color="black", lw=2.2, label="mean", zorder=5)
    ax_vol.fill_between(ts, mean_v - std_v, mean_v + std_v,
                        color="black", alpha=0.12, label="±1 std")

    ax_tvol.plot(ts, mean_tv, color="black", lw=2.2, ls="--", zorder=5)
    ax_tvol.fill_between(ts, mean_tv - std_tv, mean_tv + std_tv,
                         color="black", alpha=0.12)

    ax_vol.set_ylabel("Volume (voxels)", fontsize=11)
    ax_vol.set_title("Per-cell volume trajectories", fontsize=12)
    #ax_vol.legend(fontsize=7, ncol=min(4, max(1, n_cells // 5 + 1)),
    #              loc="upper left", framealpha=0.7)
    ax_vol.grid(True, alpha=0.3)

    ax_tvol.set_ylabel("Target volume (voxels)", fontsize=11)
    ax_tvol.set_xlabel("MCS", fontsize=11)
    ax_tvol.set_title("Per-cell target volume (dashed) trajectories", fontsize=12)
    ax_tvol.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Saved: {out_path}")
    plt.close(fig)


# ── Summary table ─────────────────────────────────────────────────────────────

def print_summary(mcs_per_cell, vol_per_cell, tvol_per_cell):
    print(f"\n{'─'*65}")
    print(f"{'cell':>6}  {'n_pts':>5}  {'tv_init':>8}  "
          f"{'v_init':>8}  {'v_final':>8}  {'Δv':>8}")
    print(f"{'─'*65}")
    for cid in sorted(vol_per_cell.keys()):
        vs   = vol_per_cell[cid]
        tvs  = tvol_per_cell[cid]
        if not vs:
            continue
        dv = vs[-1] - vs[0] if len(vs) > 1 else 0
        print(f"{cid:>6}  {len(vs):>5}  {tvs[0]:>8}  "
              f"{vs[0]:>8}  {vs[-1]:>8}  {dv:>+8}")
    print(f"{'─'*65}\n")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plot per-cell volume growth from cpm3d state JSON files")
    parser.add_argument("data_dir",
                        help="Directory containing state_mcs*.json files")
    parser.add_argument("--out", default=None,
                        help="Output figure path (default: <data_dir>/cell_growth.png)")
    args = parser.parse_args()

    records = load_states(args.data_dir)
    print(f"Loaded {len(records)} snapshots, "
          f"MCS {records[0]['mcs']}–{records[-1]['mcs']}")

    mcs_pc, vol_pc, tvol_pc = extract_trajectories(records)

    print_summary(mcs_pc, vol_pc, tvol_pc)

    out_path = args.out or os.path.join(args.data_dir.rstrip("/"), "cell_growth.png")
    plot_growth(mcs_pc, vol_pc, tvol_pc, out_path)


if __name__ == "__main__":
    main()
