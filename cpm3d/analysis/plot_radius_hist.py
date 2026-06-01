#!/usr/bin/env python3
"""
Usage: python analysis/plot_radius_hist.py <state_mcs*.json> [--bins N] [--out fig.png]
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("json")
ap.add_argument("--bins", type=int, default=20)
ap.add_argument("--out", default=None)
args = ap.parse_args()

with open(args.json) as f:
    state = json.load(f)

cells = [c for c in state["cells"] if c["id"] > 0 and c["alive"] and c["volume"] > 0]
if not cells:
    sys.exit("No living cells found.")

volumes = np.array([c["volume"] for c in cells])
radii   = (3.0 * volumes / (4.0 * np.pi)) ** (1.0 / 3.0)

out = args.out or Path(args.json).with_suffix(".radius_hist.png")

fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(radii, bins=args.bins, color="steelblue", edgecolor="white", linewidth=0.5)
ax.set_xlabel("Radius  (3V/4π)^{1/3}  (voxels)", fontsize=12)
ax.set_ylabel("Cell count", fontsize=12)
ax.set_title(f"Cell radius distribution  (n={len(radii)}, MCS={state['mcs']})", fontsize=12)
ax.axvline(radii.mean(), color="crimson", lw=1.5, ls="--", label=f"mean = {radii.mean():.2f}")
ax.legend(fontsize=10)
fig.tight_layout()
fig.savefig(out, dpi=150)
print(f"→ {out}")
