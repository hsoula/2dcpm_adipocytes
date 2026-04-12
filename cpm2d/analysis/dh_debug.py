import pandas as pd
data_dir = "../results"

df = pd.read_csv(f"{data_dir}/debug_energies.csv")

# 1. Are accepted moves dominated by area term?
print(df[df.accepted & df.area_dominated][['dh_contact','dh_area_new','dh_total']].describe())

# 2. How many runaway negative dH moves were accepted?
print(df[df.accepted & df.negative_runaway].shape[0], "runaway accepts out of", df.accepted.sum())

# 3. Distribution of dh_total for accepted moves
df[df.accepted]['dh_total'].hist(bins=50)

# 4. Which cells are the aggressors (gaining pixels most)?
df[df.accepted].groupby('s_new')['dh_area_new'].mean()