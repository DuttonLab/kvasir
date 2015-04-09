#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

from __future__ import print_function
from pymongo import MongoClient
import subprocess

db = MongoClient().demeter_test_db

all_species = db.collection_names(False)

print(all_species)

for species in all_species:
    current_species_collection = db[species]
    print(current_species_collection)
    for gene in current_species_collection.find():
        print('>', gene['species'], ' | ', gene['locus_tag'], ' | ', gene['annotation'], '\n', gene['translation'], sep='')


#print current_species