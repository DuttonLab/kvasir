
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

import pymongo
import os
from itertools import groupby
from operator import itemgetter
from bson.objectid import ObjectId
import re
import pandas as pd

def output_tsv(mongo_db_name):
    client = pymongo.MongoClient()
    db = client[mongo_db_name]

    all_species = db.collection_names(False)

    for species in all_species:
        with open('kvasir/{0}_hits.tsv'.format(species), 'w+') as output_handle:
            output_handle.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n'.format(
                'parent_locus',
                'parent_annotation',
                'parent_seq',
                'parent_contig',
                'parent_loc',
                'hit_species',
                'hit_tag',
                'hit_annotation',
                'hit_seq',
                'hit_contig',
                'hit_loc',
                )
            )

            current_species_collection = db[species]
            for record in current_species_collection.find():
                if record['hits']:
                    for hit in record['hits']:
                        hit_db_record = db[hit['hit_species']].find_one({'_id':ObjectId(hit['hit_id'])})
                        output_handle.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n'.format(
                            record['locus_tag'],
                            record['annotation'],
                            record['dna_seq'],
                            record['contig'],
                            record['location'],
                            hit_db_record['species'],
                            hit_db_record['locus_tag'],
                            hit_db_record['annotation'],
                            hit_db_record['dna_seq'],
                            hit_db_record['contig'],
                            hit_db_record['location'],
                            )
                        )
                            
def output_islands(db, species):
    current_species_collection = db[species]
    islands = []
    species_hits = []
    
    for record in current_species_collection.find():
        if record['hits']:
            species_hits.append(
                gene(
                (record['species'], str(record['_id'])),
                record['contig'],
                record['location']
                )
            )
    for entry_1 in species_hits:
        entry_recorded = -1
        for entry_2 in species_hits:
            if entry_1 == entry_2:
                pass
            elif entry_1.contig != entry_2.contig:
                pass
            elif abs(entry_1.location_start - entry_2.location_end) <= 5000:
                entry_recorded = 1
                islands.append([entry_1.identifier, entry_2.identifier])
        if entry_recorded == -1:
            islands.append([entry_1.identifier])

    return collapse_lists(islands)

class gene(object):
    def __init__(self, identifier, contig, location):
        self.identifier = identifier
        self.contig = contig
        self.location = location
        
        location_parse = re.search(r'\[(\d+)\:(\d+)\](\S+)', location)                
        self.location_start = int(location_parse.group(1))
        self.location_end = int(location_parse.group(2))
        self.location_direction = location_parse.group(3)

def output_fasta(mongo_db_name):
    client = pymongo.MongoClient()
    db = client[mongo_db_name]

    all_species = db.collection_names(False)

    for species in all_species:
        current_species_collection = db[species]
        with open('hits_fasta/{0}_hits.fna'.format(species), 'w+') as output_handle:
            for record in current_species_collection.find():
                if record['hits']:
                    output_handle.write(
                        '>{0}|{1}|{2}\n{3}\n'.format(
                            record['species'],
                            record['locus_tag'],
                            str(record['_id']),
                            record['dna_seq'],
                            )
                        )

def output_groups(mongo_db_name):
    client = pymongo.MongoClient()
    db = client[mongo_db_name]

    all_species = db.collection_names(False)
    
    all_hits_lists = []
    for species in all_species:
        current_species_collection = db[species]
        all_hits_lists.append(output_islands(db, species))

    merged_hits = []
    for species_queries in all_hits_lists:
        for query_set in species_queries:
            new_set = set()
            for query in query_set:
                query_hits = get_mongo_record(db, query)['hits']
                for hit in query_hits:
                    new_set.update([(hit['hit_species'], hit['hit_id'])])
            query_set.update(new_set)
            merged_hits.append(list(query_set))
    
    groups = map(list, collapse_lists(merged_hits))

    group_no = 0
    with open('kvasir/groups5000.tsv', 'w+') as output_handle:
        output_handle.write(
            '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(
            'groups',
            'species',
            'locus_tag',
            'contig',
            'location',
            'annotation',
            'dna_seq'
            )
        )
        for group in groups:
            group_no += 1
            for entry in group:
                db_handle = get_mongo_record(db, entry)
                output_handle.write(
                    'Group {0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(
                    str(group_no).zfill(2),
                    db_handle['species'],
                    db_handle['locus_tag'],
                    db_handle['contig'],
                    db_handle['location'],
                    db_handle['annotation'],
                    db_handle['dna_seq']
                    )
                )

def get_mongo_record(db, id_tuple):
    '''id_tuple should be (species, _id), where _id is the string 
        representation of an ObjectId''' 
    species, mongo_id = id_tuple
    current_species_collection = db[species]
    return current_species_collection.find_one({'_id':ObjectId(mongo_id)})

def collapse_lists(list_of_lists):
    # example input: [[1,2,3],[3,4],[5,6,7],[1,8,9,10],[11],[11,12],[13],[5,12]]
    result = []
    for d in list_of_lists:
        d = set(d)

        matched = [d]
        unmatched = []
        # first divide into matching and non-matching groups
        for g in result:
            if d & g:
                matched.append(g)
            else:
                unmatched.append(g)
        # then merge matching groups
        result = unmatched + [set().union(*matched)]
    return result

def get_tag_int(locus_tag):
    return int(locus_tag[-5:])

def output_compare_matrix(mongo_db_name):
    client = pymongo.MongoClient()
    db = client[mongo_db_name]

    all_species = db.collection_names(False)
    
def pair_compare(first_species, second_species, db):



# For testing
output_compare_matrix('full_pipe_test')

#if __name__ == '__main__':
#    import sys
#    output_fasta(sys.argv[1])
