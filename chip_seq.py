import pandas as pd
import json

PeakDNAseq = read_csv(Peak_DNA_seq.csv, sep=';')
df=pd.read_csv("df.csv",  sep=',')
df = pd.merge(df, PeakDNAseq)

df=df.rename(columns={'X.10.log10.pvalue.': 'X_10_log10_pvalue'})
df=df.rename(columns={'FDR...': 'FDR'})
df=df.drop('V4', 1)
df=df.drop('V5', 1)
df=df.drop('seqnames', 1)
df=df.drop('Unnamed: 0', 1)

df=df.rename(columns={'x': 'DNA_sequence'})
df = df.set_index('SYMBOL')

df['annotation'] = df.annotation.str.split(' ')
df['annotation'] = tuple(x[0] for x in df['annotation'])

with open('tmp.json', 'w') as f:
    f.write(df.to_json(orient='records'))
    
    
from pymongo import MongoClient

client = MongoClient('localhost',27017)
db=client.chip
collection=db.peaks
config=json.loads(open('tmp.json').read())

for x in config:
    collection.insert(x)
