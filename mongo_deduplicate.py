#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
Checks for and eliminates duplicate gene entries in database. 
'''

def dedupe(mongo_db_name):
    from pymongo import MongoClient

    client = MongoClient()
    db = client[mongo_db_name]
    all_species = db.collection_names(False)

    #for species in all_species:
    current_species_collection = db['brachy_IMG']

    for gene in current_species_collection.find():
        current_tag, current_id = gene['locus_tag'], gene['_id']
        current_tags = current_species_collection.find({'locus_tag':current_tag})
        if current_tags.count() > 1:
            for dupe in current_tags.distinct('_id')[1:]:
                print 'removing:'
                print dupe
                current_species_collection.remove({'_id':dupe})

# For testing
dedupe('refactor_test')
