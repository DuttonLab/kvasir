##!/usr/bin/env python

from pymongo import MongoClient
from subprocess import Popen
import os

mongod = Popen(
        ["mongod", "--dbpath", '/Users/KBLaptop/computation/db/'],
    )
client = MongoClient()
db = client['location_import_test']
collection = db['test_collection']

all_species = db.collection_names(False)
    
for species in all_species:
    current_species_collection = db[species]
    for gene in current_species_collection.find():
        print gene['location']

        
mongod.terminate()