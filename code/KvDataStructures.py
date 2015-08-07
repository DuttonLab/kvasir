#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''Class objects for use in other scripts, and functions for extracting information from MongoDB'''

import pymongo
from bson.objectid import ObjectId
from operator import itemgetter
import re


class mongo_iter(object):
    """Iterator that steps through species in mongoDB. Call with:
    `for current_species_collection in mongo_iter(mongo_db_name):`"""
    def __init__(self):
        self.index = -1
        self.collections = get_species_collections()

    def __iter__(self):
        return self

    def next(self):
        if self.index == len(self.collections) - 1:
            raise StopIteration
        else:
            self.index += 1
            return get_collection(self.collections[self.index])
                
client = pymongo.MongoClient()
db = None

def get_gene_location(location):
    '''
    Takes BioPython location object eg. `[1030:1460](-)` and extracts integer values.
    Returns (start, stop, direction)
    '''
    location_parse = re.search(r'\[<?(\d+)\:>?(\d+)\]\((\+|-)\)', str(location))
    if location_parse.group(3) == '+':
        direction = 1
    elif location_parse.group(3) == '-':
        direction = -1
    return (int(location_parse.group(1)), int(location_parse.group(2)), direction)

def mongo_init(mongo_db_name):
    global db
    db = client[mongo_db_name]
    return db.name

def get_collections():
    return db.collection_names(False)

def get_collection(collection):
    return db[collection]

def remove_collection(collection):
    print get_collections()
    db.drop_collection(collection)

def reset_database(database):
    mongo_init(database)
    print "Now you see it:"
    print db.collection_names(False)
    for collection in db.collection_names(False):
        db.drop_collection(collection)
    print "Now you don't!"
    print db.collection_names(False)

def get_species_collections():
    collections = db.collection_names(False)
    if '16S' in collections:
        collections.remove('16S')
    if 'hits' in collections:
        collections.remove('hits')
    return collections

def get_mongo_record(species, mongo_id):
    species_collection = get_collection(species)
    return species_collection.find_one({'_id':ObjectId(mongo_id)})

def get_genus(species_name):
    return re.search(r'^(\w+)', species_name).group(1)

def view_record():
    for item in get_species_collections():
        print item
        counter = 1
        for collection in get_collection(item).find():
            counter += 1
            if counter > 3:
                break

def fasta_id_parse(fasta_id):
    """"species_name|_id" -> (species_name, _id)"""
    parsed = re.search(r'(\w+)\|(\w+)', fasta_id)
    return (parsed.group(1), parsed.group(2))

def index_contigs(species_collection):
    return species_collection.find().sort([("location.contig", pymongo.ASCENDING), ("location.start", pymongo.ASCENDING)])

def concat_contigs(species_collection):
    current_contig = 0
    last_contig_end = 0
    last_gene_end = 0
    concatenated = {}
    for gene in index_contigs(species_collection):
        if gene['location']['contig'] == current_contig:
            gene['location']['start'] += last_contig_end
            gene['location']['end'] += last_contig_end
            last_gene_end = gene['location']['end']
            
            concatenated[gene['_id']] = gene
        else:
            current_contig = gene['location']['contig']
            last_contig_end = last_gene_end

            gene['location']['start'] += last_contig_end
            gene['location']['end'] += last_contig_end
            last_gene_end = gene['location']['end']
        
            concatenated[gene['_id']] = gene

    return concatenated


if __name__ == '__main__':
    db = str(raw_input("Please enter database name: "))
    mongo_init(db)
    print get_species_collections()    
    