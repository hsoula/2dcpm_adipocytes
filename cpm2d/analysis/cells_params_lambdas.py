import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

data_dir = "../results/lambdas"

df = pd.read_csv(f'{data_dir}/out.csv')
df['iso'] = df['perimeter']**2/df['area']

dtfmean = df.groupby(['lambda_area','lambda_perimeter','lambda_iso','target_area','n_cells']).mean()


perims = np.unique(df['lambda_perimeter'])

dtfmean = df.groupby(['lambda_perimeter']).mean()

plt.scatter(dtfmean.index,dtfmean['perimeter'])
#plt.plot(dtfmean.index, 4*np.pi + 0 * dtfmean.index)
plt.savefig('perim.pdf')
plt.show()


##@ 
# 00000000
# 00011100
# 00111110
# 01111110
# 00111100
# 00011000
# 00000000

# 2+4+6+5+3 = 20 
# 22
