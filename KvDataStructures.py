#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''Class objects for use in other scripts, and functions for extracting information from MongoDB'''

import pymongo
from bson.objectid import ObjectId
import re

class gene_location(object):
    '''Takes BioPython location object eg. `[1030:1460](-)` and extracts integer values'''
    def __init__(self, location):
        self.location = location

        location_parse = re.search(r'\[(\d+)\:(\d+)\](\S+)', location)                
        self.start = int(location_parse.group(1))
        self.end = int(location_parse.group(2))
        self.direction = location_parse.group(3)

class mongo_iter(object):
    """Iterator that steps through species in mongoDB. Call with:
    `for current_species_collection in mongo_iter(mongo_db_name):`"""
    def __init__(self):
        self.index = -1
        self.collections = get_collections()

    def __iter__(self):
        return self

    def next(self):
        if self.index == len(self.collections) - 1:
            raise StopIteration
        else:
            self.index += 1
            return get_species_collection(self.collections[self.index])
                
client = pymongo.MongoClient()
db = None

def mongo_init(mongo_db_name):
    global db
    db = client[mongo_db_name]

def get_collections():
    return db.collection_names(False)

def get_species_collection(species):
    return db[species]

def get_mongo_record(species_collection, mongo_id):
    return species_collection.find_one({'_id':ObjectId(mongo_id)})
     
        