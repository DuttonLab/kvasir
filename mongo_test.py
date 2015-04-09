#!/usr/bin/env python

from pymongo import MongoClient

db = MongoClient().demeter_test_db

all_species = db.collection_names(False)

print all_species

for species in all_species:
    current_species_collection = db[species]
    print current_species_collection
    for gene in current_species_collection.find():
        print gene

    #print current_species