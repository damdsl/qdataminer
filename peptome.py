
import pandas as pd
import toolz
import json
from pymongo import MongoClient
from bson import ObjectId
import re

df = pd.read_csv({MaxQuant_document_path.csv}, sep=',')
df = df.set_index('Sequence')

#Cleaning
#Protein names are splitted and saved in an array
df["Proteins"] = df.Proteins.str.split(";")

#replace NaN in Protein name field by an empty array
null_values = df.Proteins.isnull()
for idx in null_values[null_values].index:
    df['Proteins'].set_value(idx, [])


#Extraction of intensity per QBiC samples identifiers (Barcode) and save it in a dataframe
cols = [col for col in df.columns if col.startswith('Intensity ')]
cols = [col.split() for col in cols if len(col.split()) == 2]
samples = toolz.pluck(1, cols)
samples = [sample for sample in samples if sample.startswith('Q')]
intensity_cols = ['Intensity ' + sample for sample in samples]
intensities = df[intensity_cols]
intensities.columns = samples

#Extraction of the fields of interest 
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
    
#Concatenation of intensity dataframe and cols_selection dataframe in order to get one dataframe with all quantitative informations per sample
#protein list is embedded in a group protein sub-document contatining specific annotations (uniprot code and fasta header) 
for col in intensities.columns:
    df3 = intensities[col]
    result = pd.concat([df3, cols_selection], axis=1, join_axes=[cols_selection.index])
    result["sample"] = col
    result.columns = ["intensity", "N-term_cleavage_window", "C-term_cleavage_window", "length", "mass", "proteinIds", "acetyl(K)_siteIDs", "oxidation(M)_siteIDs", "phospho(STY)_siteIDs", "sample"] 
    result['proteinGroup'] = result.proteinIds.apply(convert_headers)
    result  = result.drop("proteinIds", 1)
    result = result.reset_index()
    #results are saved in a json file / one file per sample containing the list of sequenced peptides and related quantitative informations
    with open('tmp_%s.json' % col, 'w') as f:
        f.write(result.to_json(orient='records'))


