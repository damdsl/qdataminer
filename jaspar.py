for group in df.index:
    data={}
    meta = matrix_annotation2.get_group(group)
    meta=meta.set_index('index')
    data=meta.to_json(orient='records'))

for group in df.index:
    data={}
    data= df.loc[group, ['collection', 'matrixID', 'nameTF', 'species', 'acc']].to_json()
    meta = matrix_annotation2.get_group(group)
    meta=meta.drop('number', 1)
    meta=meta.set_index('index')
    metainformation=meta.to_json()
    print(group)
    files=df.loc[group,'matrixID']
    print(files)
    with open("/root/jaspar_database/jaspar_database/sites/%s.1.sites" % files, 'w') as handle:
        motif = motifs.read(handle, "sites")
        motifs=[]
        for instance in motif.instances:
            motifs.append(str(instance))
    
