import pandas as pd
import os

namedict=pd.read_csv('GkountelaNamedict.csv')

runName=namedict.Run
newName=namedict.filename

d=dict(zip(runName,newName))


cwd=os.getcwd()

base=cwd+'/'

l=os.listdir()
print(base)

for i in l:
        src=base+i
        if i in d.keys():
                dst=d[i]
                os.rename(src,dst)

