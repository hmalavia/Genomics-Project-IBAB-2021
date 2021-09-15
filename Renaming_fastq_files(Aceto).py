import pandas as pd
import os

df=pd.read_csv('SraRunTable.csv')

#df.head(10)


dfFilt=df[['Run','patient','also_known_as','single_cells_or_cluster']]

#dfFilt

dfFilt['filename']=dfFilt['patient']+dfFilt['single_cells_or_cluster']

c=dfFilt.groupby('filename').cumcount()+1

dfFilt['c']=c

dfFilt['newfilename']=dfFilt['filename'] + dfFilt['c'].astype('str')

ogfilename=dfFilt['Run']+'.fastq'
newfilename=dfFilt['newfilename']+'.fastq'

d=dict(zip(ogfilename,newfilename))
d

for a,b,c in os.walk('/home/group_viji02/Harsh_Arnav_PGDB/DataSets/Aceto2014/'):
    os.chdir(a)
    cwd=a
    l=os.listdir(cwd)
    for i in l:
        base=a
        src=base+'/'+i
        if i in d.keys():
            dst=d[i]
            os.rename(src,dst)

    os.chdir('../')
