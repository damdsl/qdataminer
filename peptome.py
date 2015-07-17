
import pandas as pd
import toolz
import json
from pymongo import MongoClient
from bson import ObjectId
import re

df = pd.read_csv('/home/elasticsearch/A84D-39CD/example_protein_lists_MaxQuant/peptides_SGPN.txt', sep=',')
df["Proteins"] = df.Proteins.str.split(";")
df = df.set_index('Sequence')

null_values = df.Proteins.isnull()
for idx in null_values[null_values].index:
    df['Proteins'].set_value(idx, [])

cols = [col for col in df.columns if col.startswith('Intensity ')]
cols = [col.split() for col in cols if len(col.split()) == 2]
samples = toolz.pluck(1, cols)
samples = [sample for sample in samples if sample.startswith('Q')]

intensity_cols = ['Intensity ' + sample for sample in samples]
intensities = df[intensity_cols]
intensities.columns = samples

cols_selection = df.loc[:,['N-term cleavage window','C-term cleavage window', 'Length', 'Mass', 'Proteins', 'Acetyl (K) site IDs', 'Oxidation (M) site IDs', 'Phospho (STY) site IDs']]


uniprot_re = re.compile('\S*\|([a-zA-Z0-9-]+)\|')

def extract_prot_id(header):
    match = re.match(uniprot_re, header)
    if not match:
        #match = re.match('CON__([a-zA-Z0-9-]+)', header)
        #if not match:
        return None
    return match.groups()[0]

def convert_headers(headers):
    ids = [extract_prot_id(header) for header in headers]
    return [{'header': header,
             'uniprot': id} for header, id in zip(headers, ids)]
    
client = MongoClient('localhost', 27017)
db = client.peptome
analysis = db.peptide_analysis
 
for col in intensities.columns:
    df3 = intensities[col]
    result = pd.concat([df3, cols_selection], axis=1, join_axes=[cols_selection.index])
    result["sample"] = col
    result.columns = ["intensity", "N-term_cleavage_window", "C-term_cleavage_window", "length", "mass", "proteinIds", "acetyl(K)_siteIDs", "oxidation(M)_siteIDs", "phospho(STY)_siteIDs", "sample"] 
    result['proteinGroup'] = result.proteinIds.apply(convert_headers)
    result  = result.drop("proteinIds", 1)
    result = result.reset_index()

    with open('tmp_%s.json' % col, 'w') as f:
        f.write(result.to_json(orient='records'))
    config = json.loads(open(f.name).read())
    for x in config:
        analysis.insert(x)
