#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''Class objects for use in other scripts, and functions for extracting information from MongoDB'''

import pymongo

class gene(object):
    def __init__(self, identifier, contig, location):
        self.identifier = identifier
        self.contig = contig
        self.location = gene_location(location)

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
    `for current_species collection, species in mongo_iter(mongo_db_name):`"""
    def __init__(self, db):
        self.db = db
        self.index = -1
        self.collections = get_collections(db)

    def __iter__(self):
        return self

    def next(self):
        if self.index == len(self.collections) - 1:
            raise StopIteration
        else:
            self.index += 1
            return get_species_collection(self.db, self.collections[self.index])
                
            
def get_collections(mongo_db_name):
    client = pymongo.MongoClient()
    db = client[mongo_db_name]
    return db.collection_names(False)

def get_species_collection(mongo_db_name, species):
    client = pymongo.MongoClient()
    db = client[mongo_db_name]
    return db[species]

for item in enumerate(mongo_iter('full_pipe_test')):
    print item[1].name


            
        