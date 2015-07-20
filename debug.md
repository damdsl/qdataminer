In case data are not transfered from mongodb to elasticsearch by transporter

error message: 

#Use elasticsearch to 

#create a tmp database with search results:

client = MongoClient('localhost', 27017)
db = client.search
data = db.data

for col in intensities.columns:
    for i in result["Sequence"]:
        sequence = i
        #Find in elasticsearch what documents are not transferred
        res = es.search(index="peptome", body={
            "query": {
                "bool": {
                    "must": [
                        {
                            "query_string": {
                                "default_field": "myindex.myfield",
                                "query": sequence
                            }
                        }
                        ,
                        {
                            "query_string": {
                                "default_field": "mytype.myfield2",
                                "query": col
                            }
                        }
                    ],
                    "must_not": [ ],
                    "should": [ ]
                }
            },
            "from": 0,
            "size": 10,
            "sort": [ ],
            "facets": { }
        })
        #Save information with unique identifier permitting to retrieve in mongodb none-transferred documents
        new_doc = {}
        new_doc["query_sequence"] = sequence
        new_doc["query_sample"] = col
        new_doc["ES_query_return"] = res
        data.insert(new_doc)
