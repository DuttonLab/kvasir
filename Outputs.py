#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

import pymongo
import os
from itertools import groupby, combinations
import pandas as pd
import KvDataStructures as kv

def output_tsv(mongo_db_name):
    for current_species_collection in kv.mongo_iter(mongo_db_name):
        species = current_species_collection.name
        with open('kvasir/{0}_hits.tsv'.format(species),'w+') as output_handle:
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
    for current_species_collection in kv.mongo_iter(mongo_db_name):
        species = current_species_collection.name
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

def get_groups():
    groups_list = []
    for current_species_collection in kv.mongo_iter():
        # Generate list of islands, formated as list of lists of gene ids. Each inner list represents one island eg: 
        # [[gene1, gene2, gene3],[gene4, gene5],[gene6, gene7]]
        current_species_islands = get_islands(current_species_collection)

        # each sublist represents one island...
        for island in current_species_islands:
            hit_set = set() # container for hits 
            for gene_id in island:
                gene_hits = kv.get_mongo_record(current_species_collection, gene_id[1])['hits']
                
                # Pulls each hit and builds id tuple, then appends it to group_set
                for hit in gene_hits:
                    hit_species_collection = kv.get_species_collection(hit['hit_species'])
                    hit_db_record = kv.get_mongo_record(hit_species_collection, hit['hit_id'])
                    hit_set.add(
                        tuple((hit_db_record['species'], str(hit_db_record['_id'])))
                        )
            # add id tuples for hits to island list...
            island.update(hit_set)
            # And add new island (with multiple species) to groups_list
            groups_list.append(list(island))

    # Since each species' islands are built independently, there's a lot of redundancy
    # So... Collapse lists that contain shared elements and deduplicate
    return map(list, collapse_lists(groups_list))

def output_groups(groups_list, output_file='tmp/groups_test.tsv'):
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
        
        group_no= 0
        for group in groups_list:
            group_no += 1
            # Entry is `(species, id)`
            for entry in group:
                species_collection = kv.get_species_collection(entry[0])
                db_handle = kv.get_mongo_record(species_collection, entry[1])
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

def output_compare_matrix(groups_list): #needs fixing
    
    with open('tmp/compare_matrix.tsv', 'w+') as output_handle:
        output_handle.write('Species 1\tSpecies 2\tshared CDS\tshared nt\tshared groups\n')
        for species_pair in combinations(kv.get_collections(), 2):
            comparison = pair_compare(species_pair[0], species_pair[1])
            
            shared_groups = 0
            for group in groups_list:
                if any(species_pair[0] in x for x in group):
                    if any(species_pair[1] in y for y in group):
                        shared_groups += 1


            output_handle.write(
                '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(
                    species_pair[0],
                    species_pair[1],
                    comparison[0],
                    comparison[1],
                    shared_groups
                    )
                )
    
def pair_compare(species_1, species_2):    
    shared_CDS = 0
    shared_nt = 0

    for species_1_record in kv.get_species_collection(species_1).find():
        if species_1_record['hits']:
            for hit in species_1_record['hits']:
                if hit['hit_species'] == species_2:
                    shared_CDS += 1
                    species_2_record = kv.get_mongo_record(
                        kv.get_species_collection(hit['hit_species']),
                        hit['hit_id']
                    )
                    hit_loc = kv.gene_location(species_2_record['location'])
                    shared_nt += (hit_loc.end - hit_loc.start)
    return (shared_CDS, shared_nt)


def get_islands(species_collection):
    islands = []
    species_hits = []
    
    # Add mongo_record for each hit in any gene
    for record in species_collection.find():
        if record['hits']:
            species_hits.append(
                kv.get_mongo_record(
                    kv.get_species_collection(record['species']), str(record['_id'])
                )
            )
    for entry_1 in species_hits:
        entry_recorded = -1
        for entry_2 in species_hits:
            if entry_1 == entry_2:
                pass
            elif entry_1['contig'] != entry_2['contig']:
                pass
            else:
                location_1 = kv.gene_location(entry_1['location'])
                location_2 = kv.gene_location(entry_2['location'])
                if abs(location_1.start - location_2.end) <= 5000:
                    entry_recorded = 1
                    islands.append([
                        (entry_1['species'], str(entry_1['_id'])),
                        (entry_2['species'], str(entry_2['_id']))
                    ])
        if entry_recorded == -1:
            islands.append([(entry_1['species'], str(entry_1['_id']))])

    return collapse_lists(islands)

"""Basic Use Functions"""
def collapse_lists(list_of_lists):
    # example input: [[1,2,3],[3,4],[5,6,7],[1,8,9,10],[11],[11,12],[13],[5,12]]
    # example output: [[1,2,3,4,8,9,10],[5,6,7,11,12],[13]]
    # from stackoverflow user YXD: http://stackoverflow.com/questions/30917226/collapse-list-of-lists-to-eliminate-redundancy
    result = []
    for l in list_of_lists:
        s = set(l)

        matched = [s]
        unmatched = []
        # first divide into matching and non-matching groups
        for g in result:
            if s & g:
                matched.append(g)
            else:
                unmatched.append(g)
        # then merge matching groups
        result = unmatched + [set().union(*matched)]
    return result

def get_tag_int(locus_tag):
    return int(locus_tag[-5:])

# For testing
kv.mongo_init('full_pipe_test')
output_compare_matrix(get_groups())



#if __name__ == '__main__':
#    import sys
#    pair_group_compare(sys.argv[1])
