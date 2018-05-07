import numpy as np
import math as m
import pandas as pd

with open('UV-Vis 4_26_2018 2_24_50 PM.tsv') as f:
    lines = f.readlines()

data2 = pd.read_csv('UV-Vis 4_26_2018 2_24_50 PM.tsv', sep = '\t', names
        = list(range(0,2)))

wavelenghts = data2[0][10:1330]



#repeat every 1333
data = [x.strip() for x in lines]

n =m.ceil(len(data)/1333)

col = []

for i in range(n):
    col.append(data2.loc[1331*(i),0])

abval = pd.Series()

for i in range(n):
    abval += data2.loc[10 +1331*i:1329 + 1331*i,1]


#index of first data point
start = 10

#index of last data point
stop = 1330

#ormalize to interband at 380 nm 
interband = 390

data_name_list = []

for i in range(n):
    data_name_list.append(data[i*1333])

print(data_name_list)

data_list = [i for i in data[10:1330]]



