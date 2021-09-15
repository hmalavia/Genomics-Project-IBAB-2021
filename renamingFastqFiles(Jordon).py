import pandas as pd
import os
data = pd.read_csv('JordonNameDict.csv')

run = data.org
filename = data.new

d = dict(zip(run, filename))

cwd = os.getcwd()
base = cwd+'/'
for i in os.listdir(cwd):
    src = base+i
    if i in d.keys():
        dst = d[i]
        print('{}->{}'.format(src, dst))