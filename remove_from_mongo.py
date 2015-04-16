#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
Checks for and eliminates duplicate gene entries in database. 
'''

def remove(mongo_db_name):
    from pymongo import MongoClient

    client = MongoClient()
    db = client[mongo_db_name]
    all_species = db.collection_names(False)

    for species in all_species:
        current_species_collection = db[species]
        for gene in current_species_collection.find():
            print gene['species']
        

# For testing
remove('refactor_test')
