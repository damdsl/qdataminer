import pandas as pd
import json

df=pd.read_csv("/home/damiendsl/data_analysis_toolbox/datasets/GSE71250/macs_SRR2124925/df.csv",  sep=',')
df=df.rename(columns={'X.10.log10.pvalue.': 'X_10_log10_pvalue'})
df=df.rename(columns={'FDR...': 'FDR'})
df=df.drop('V4', 1)
df=df.drop('V5', 1)
df=df.drop('seqnames', 1)
df=df.drop('Unnamed: 0', 1)
df = pd.merge(df, PeakDNAseq)
df=df.rename(columns={'x': 'DNA_sequence'})
df = df.set_index('SYMBOL')

with open('tmp.json', 'w') as f:
    f.write(df.to_json(orient='records'))
    
    
from pymongo import MongoClient

client = MongoClient('localhost',27017)
db=client.chip
collection=db.peaks
config=json.loads(open('tmp.json').read())

for x in config:
    collection.insert(x)
