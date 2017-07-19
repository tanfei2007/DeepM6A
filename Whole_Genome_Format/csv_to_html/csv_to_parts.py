#*****************************************
# split big .csv file into small .csv file
# Xiurui Hou
# xh256@njit.edu
#*****************************************
import pandas as pd
import math
import sys

species = sys.argv[1]
df = pd.read_csv('../'+species+'_merge.csv')
line_num = df.shape[0]
part_num = int(math.ceil(line_num / 1000000.0))
for i in range(0, part_num-1):
    start = i * 1000000
    end = start + 1000000 - 1
    df.loc[start: end].to_csv('./csv_parts/'+species+'_part'+str(i)+'.csv', index = False)
i += 1
df.loc[end+1:].to_csv('./csv_parts/'+species+'_part'+str(i)+'.csv', index = False)


