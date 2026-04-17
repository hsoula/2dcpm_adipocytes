import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
rf = pd.read_csv('data/sim3d/volume_timeseries.csv')
last = max(rf.mcs)
rfe = rf[rf.mcs==last]
rfs = rf[rf.mcs==0]

plt.scatter(rfs.target_volume, rfe.volume)
plt.show()

bins = np.arange(0,6,0.1)
he, _ = np.histogram(rfe.volume**(1/3), bins=bins)
hs, _ = np.histogram(rfs.target_volume**(1/3), bins=bins)
b = 0.5 * (bins[1:]+bins[:-1])
plt.plot(b, he)
plt.plot(b, hs)
plt.legend(['final','start'])
plt.savefig('results/distrib.png')
plt.show()
