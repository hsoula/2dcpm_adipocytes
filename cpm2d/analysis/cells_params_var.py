import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

data_dir = "../results"

df = pd.read_csv(f'{data_dir}/sweep_var.csv')
df.drop(columns='lambda_area', inplace=True)
df = df[df['lambda_iso']==0]
df['iso'] = df['perimeter']**2/df['area']

# plt.scatter(dtfmean.index,dtfmean['perimeter'])
# #plt.plot(dtfmean.index, 4*np.pi + 0 * dtfmean.index)
# plt.savefig('perim.pdf')
# plt.show()


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
