for group in df.index:
    data={}
    meta = matrix_annotation2.get_group(group)
    meta=meta.set_index('index')
    data=meta.to_json(orient='records'))
