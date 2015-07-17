# qdataminer

This repo contains python scripts to process and extract valuable quantitative information from MaxQuant (MQ) proteomics reports

peptome script takes MQ quantitative result files as an input, extracts values of interest and returns one json file per sample containing a list of json documents with some pre-selected metadata informations and quantified intensity results. With such format data can be loaded through iteration into a document oriented database such MongoDB.

proteome script creates a json file per analyzed samples containing its entire proteome (expressed protein list with measured intensity)
