#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

import pymongo
import os
import re
from itertools import groupby, combinations
from bson.objectid import ObjectId
import pandas as pd
import KvDataStructures as kv

def output_tsv(mongo_db_name):
    for current_species_collection, species in kv.mongo_iter(mongo_db_name):
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

def output_fasta(mongo_db_name):
    for current_species_collection, species in kv.mongo_iter(mongo_db_name):
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

def output_groups(mongo_db_name, output_file='kvasir/groups5000.tsv'):
    hits_list = []
    for current_species_collection, species in kv.mongo_iter(mongo_db_name):
        current_species_islands = get_islands(current_species_collection)

        for island in current_species_islands:
            new_set = set()
            for gene_id in island:
                gene_hits = get_mongo_record(current_species_collection, gene_id)['hits']
                for hit in gene_hits:
                    new_set.update([(hit['hit_species'], hit['hit_id'])])
            current_id.update(new_set)
            hits_list.append(list(current_id))
    
    groups = map(list, collapse_lists(hits_list))

    group_no = 0
    with open(output_file, 'w+') as output_handle:
        print output_handle
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
                db_handle = get_mongo_record(, entry)
                output_handle.write(
                    '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(
                    str(group_no).zfill(2),
                    db_handle['species'],
                    db_handle['locus_tag'],
                    db_handle['contig'],
                    db_handle['location'],
                    db_handle['annotation'],
                    db_handle['dna_seq']
                    )
                )
    return pair_group_compare(all_species, output_file)

def output_compare_matrix(mongo_db_name):
    client = pymongo.MongoClient()
    db = client[mongo_db_name]

    all_species = db.collection_names(False)
    with open('kvasir/compare_matrix.tsv', 'w+') as output_handle:
        output_handle.write('Species 1\tSpecies 2\tshared CDS\tshared nt\n')
        for species_pair in combinations(all_species, 2):
            comparison = pair_compare(species_pair[0], species_pair[1], db)
            output_handle.write(
                '{0}\t{1}\t{2}\t{3}\n'.format(
                    species_pair[0],
                    species_pair[1],
                    comparison[0],
                    comparison[1],
                    )
                )
    
def pair_compare(first_species, second_species, db):
    first_species_collection = db[first_species]
    shared_CDS = 0
    shared_nt = 0
    shared_group = 0

    for first_species_record in first_species_collection.find():
        if first_species_record['hits']:
            for hit in first_species_record['hits']:
                if hit['hit_species'] == second_species:
                    shared_CDS += 1
                    second_species_record = get_mongo_record(db, (hit['hit_species'], hit['hit_id']))
                    hit_loc = gene_location(second_species_record['location'])
                    shared_nt += (hit_loc.end - hit_loc.start)
    return (shared_CDS, shared_nt)

def pair_group_compare(list_of_species, group_output_file):
    species_combos = list(combinations(list_of_species, 2))
    pair_counts = {n:0 for n in species_combos}

    df = pd.read_csv(group_output_file, sep='\t')
    df2 = df.loc[:,['groups','species']]

    for group in df2.groupby('groups'):
        group_species = list(group[1]['species'])
        for combo in species_combos:
            if combo[0] in group_species:
                if combo[1] in group_species:
                    pair_counts[combo] += 1

    return pair_counts

        #for combo in species_combos:
        #    print combo

"""Getters"""
def get_islands(species_collection):
    islands = []
    species_hits = []
    
    for record in species_collection.find():
        if record['hits']:
            species_hits.append(
                kv.gene(
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
            elif abs(entry_1.location.start - entry_2.location.end) <= 5000:
                entry_recorded = 1
                islands.append([entry_1.identifier, entry_2.identifier])
        if entry_recorded == -1:
            islands.append([entry_1.identifier])

    return collapse_lists(islands)

def get_mongo_record(species_collection, mongo_id): # Check others' use
    return species_collection.find_one({'_id':ObjectId(mongo_id)})

"""Basic Use Functions"""
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

# For testing
output_groups('full_pipe_test', '/Users/KBLaptop/googleDrive/work/duttonLab/working_database/kvasir/groups5000.tsv')
#if __name__ == '__main__':
#    import sys
#    pair_group_compare(sys.argv[1])
