#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
This script is designed to read a mongo database created with DataImport,
create a faa file, build a BLAST database with that file, then delete the 
faa file. Must have Mongod running, in terminal: `mongod --dbpath path/to/db`
'''

def db_debug(mongo_db_name):
    from pymongo import MongoClient
    import os

    client = MongoClient()
    
    hits_db = client[mongo_db_name]
    all_species = hits_db.collection_names(False)

    for species in all_species:
        current_species = hits_db[species]
        print current_species
        for entry in current_species.find():
            print entry['locus_tag']
    

db_debug('pipe_test')