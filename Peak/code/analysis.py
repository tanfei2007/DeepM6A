import pandas as pd
import numpy as np

pos = pd.read_csv('peak_positive.csv', header=0, sep=',')
neg = pd.read_csv('peak_negative.csv', header=0, sep=',')
all =  pd.concat([pos, neg])


peak_max = all.groupby(['peak'])['pred'].max()

df = pd.DataFrame()
df['peak'] = list(peak_max.index)
df['prob'] = peak_max.values

df.to_csv('peak_prob.csv', index=None)
