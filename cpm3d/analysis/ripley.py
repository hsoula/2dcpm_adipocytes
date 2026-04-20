import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform

# ---------------------------------------------------------------
def ripley_K_3d(points, r_vals, box_mins, box_maxs, box_volume,
                edge_correction='border'):
    """
    Ripley's K function for 3D points.

    Parameters
    ----------
    points      : (N,3) array
    r_vals      : 1-D array of distances
    box_mins    : (3,) lower corner of bounding box
    box_maxs    : (3,) upper corner of bounding box
    box_volume  : scalar, volume of bounding box
    edge_correction : 'border' or 'none'

    Returns
    -------
    K : 1-D array, K(r) values
    """
    N    = len(points)
    dist = squareform(pdist(points))
    K    = []

    for r_val in r_vals:
        if edge_correction == 'border':
            inside  = (np.all(points - box_mins >= r_val, axis=1) &
                       np.all(box_maxs - points >= r_val, axis=1))
            indices = np.where(inside)[0]
            if len(indices) == 0:
                K.append(np.nan)
                continue
            count = sum(
                ((dist[i] > 0) & (dist[i] <= r_val)).sum()
                for i in indices
            )
            K_val = (box_volume / (len(indices) * N)) * count
        else:
            pairs = np.sum(dist[np.triu_indices(N, k=1)] <= r_val)
            K_val = (box_volume / (N * (N - 1))) * 2 * pairs
        K.append(K_val)

    return np.array(K)

rf = pd.read_csv('data/sim3d_life/volume_timeseries.csv')

last = max(rf.mcs)
rfe = rf[rf.mcs==last]

centroids  = rfe[['cent_x', 'cent_y', 'cent_z']].astype(float).values
N          = len(centroids)
box_mins   = centroids.min(axis=0)
box_maxs   = centroids.max(axis=0)
box_volume = float(np.prod(box_maxs - box_mins))
step = 0.1 
top = 10 
r_vals = np.arange(step, top + step, step)
K_obs = ripley_K_3d(centroids, r_vals, box_mins, box_maxs, box_volume)
print(K_obs)
L_obs = (np.where(K_obs > 0, K_obs, np.nan) / ((4/3) * np.pi)) ** (1/3) - r_vals
plt.plot(r_vals, L_obs)
plt.savefig('results/ripley_life.png')
plt.show()