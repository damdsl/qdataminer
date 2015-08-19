import pandas as pd
import toolz
import json
from bson import ObjectId
import re

df = pd.read_csv('<yourfilepath>.txt', sep='\t')
df = df.set_index('Protein IDs')
df.index.name = 'proteinIds'

#Extraction of quantitative information per QBiC sample identifier
cols = [col for col in df.columns if col.startswith('Intensity ')]
cols = [col.split() for col in cols if len(col.split()) == 2]
samples = toolz.pluck(1, cols)
samples = [sample for sample in samples if sample.startswith('Q')]
intensity_cols = ['Intensity ' + sample for sample in samples]
intensities = df[intensity_cols]
intensities.columns = samples

#Preparation
flat_int = intensities.unstack()
flat_int.name = "intensity"
flat_int = flat_int.reset_index()
flat_int.columns = ["sample", "proteinIds", "intensity"]
flat_int["proteinIds"] = flat_int.proteinIds.str.split(";")

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
    
#Data are saved in a json file, one file per sample containing its quantified proteome
flat_int['proteinGroup'] = flat_int.proteinIds.apply(convert_headers)
for group, df in flat_int.groupby("sample"):
    data = {}
    data["sample"] = group
    data["proteins"] = json.loads(df[["proteinGroup", "intensity"]].to_json(orient='records'))
    with open('output_%s.json' % group, 'w') as f:
        json.dump(data, f)
        analysis.insert(data)

for group in df.index:
    data={}
    data['jaspar_number']=group
    data['annotations']= matrix_annotation2.get_group(group).to_json(orient='records')
    file= df.loc[group, 'matrixID']
    with open('/root/jaspar_database/jaspar_database/sites/%s.1.sites' % file, 'w') as handle:
        motif = motifs.read(handle, "sites")
        sequence=[]
        for instance in motif.instances:
            sequence.append(str(instance))
        data['motifs']=sequence
    with open('output_%s.json' % group, 'w') as f:
        json.dump(data, f)
