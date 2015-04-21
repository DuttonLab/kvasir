#!/usr/bin/env python

from pymongo import MongoClient

client = MongoClient()
db = client['refactor_test']
all_species = db.collection_names(False)

for species in all_species:
    current_species_collection = db[species]
    for gene in current_species_collection.find():
        print gene['species']
        print gene['locus_tag']
