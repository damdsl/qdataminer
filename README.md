# qdataminer

qdataminer is a search system for the mining of quantitative omics data

Data is extracted and processed through python scripts and then stored in a Mongodb database. Mining of data can be done through elasticsearch API synchronized to Mongodb via transporter application (Compose) 

This repo contains python scripts to process and extract valuable quantitative information (from MaxQuant (MQ) proteomics reports and chip_seq data so far)

- peptome script takes MQ quantitative result files as an input, extracts values of interest and returns one json file per sample containing a list of json documents with some pre-selected metadata informations and quantified intensity results. With such format data can be loaded through iteration into a document oriented database such MongoDB.

- proteome script creates a json file per analyzed samples containing its entire proteome (expressed protein list with measured intensity)

- chip_seq.R and chip_seq.py take BED files from MACS peak calling BED file and returns one json document with values of intersest
