#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

from __future__ import print_function
from pymongo import MongoClient
import subprocess

def create_mongo_database(database_name, path_to_database):
    mongod = subprocess.Popen(
        "mongod --dbpath {0}".format(path_to_database),
        shell=True
    )
    client = MongoClient()
    db = client[database_name]
    collection = db['test_collection']
    collection.insert_one({'something new':'here\'s something else'})
    print('err... what\'s up?')
    mongod.terminate()


create_mongo_database('test5_database', '~/computation/db')
    
#db = MongoClient().demeter_test_db
#
#all_species = db.collection_names(False)
#
#print(all_species)
#
#for species in all_species:
#    current_species_collection = db[species]
#    print(current_species_collection)
#    for gene in current_species_collection.find():
#        print('>', gene['species'], ' | ', gene['locus_tag'], ' | ', gene['annotation'], '\n', gene['translation'], sep='')


#print current_species