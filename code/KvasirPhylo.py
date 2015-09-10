#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# Unless otherwise indicated, licensed under GNU Public License (GPLv3)

import os

import KvDataStructures as kv
import pandas as pd

from itertools import combinations_with_replacement
from subprocess import Popen, PIPE
from skbio import DNA
from skbio import DistanceMatrix
from skbio import tree
from skbio.alignment import StripedSmithWaterman

def get_gene_distance(seq_1, seq_2):
    """
    Returns hamming distance between two DNA sequences
    Alignment based on Striped Smith-Waterman algorithm
    """
    query = StripedSmithWaterman(seq_1.upper())
    alignment = query(seq_2.upper())
    q = DNA(alignment.aligned_query_sequence)
    t = DNA(alignment.aligned_target_sequence)
    return q.distance(t)

def get_perc_id(seq_1, seq_2):
    return 1.0 - get_gene_distance(seq_1, seq_2)

def get_16S_distance(species_1, species_2):  
    """
    Returns hamming distance of the 16S seq of two species in MongoDB
    """
    if species_1 == species_2:
        return 0.0
    else:
        s1_ssu = kv.db['16S'].find_one({'species':species_1})
        s2_ssu = kv.db['16S'].find_one({'species':species_2})
        if s1_ssu and s2_ssu:
            if len(s1_ssu['dna_seq']) > 800 and len(s1_ssu['dna_seq']) > 800:
                return get_gene_distance(str(s1_ssu['dna_seq']), str(s2_ssu['dna_seq']))
            else:
                return None
        else:
            return None

def get_distance_matrix(core=False, to_file=True):
    all_species = kv.get_collection('core').distinct('species')
    if core:
        pass
    else:
        all_species.extend(kv.get_collection('other').distinct('species'))

    ssu_species = [n for n  in all_species if kv.db['16S'].find_one({'species':n})]
    distance_matrix = pd.DataFrame(data={n:0.0 for n in ssu_species}, index=ssu_species, columns=ssu_species)
    
    for pair in combinations_with_replacement(ssu_species, 2):
        distance = get_16S_distance(pair[0], pair[1])
        if distance:
            distance_matrix[pair[0]][pair[1]] = distance
            distance_matrix[pair[1]][pair[0]] = distance

    if to_file:
        distance_matrix.to_csv('distance_matrix.csv')
    else:
        return distance_matrix

def get_tree(core=False, newick=False):
    core_collection = kv.get_collection('core')
    all_species = core_collection.distinct('species')
    if core:
        pass
    else:
        other_collection = kv.get_collection('other')
        all_species.extend(other_collection.distinct('species'))
    ssu_species = [n for n  in all_species if kv.db['16S'].find_one({'species':n})]
    
    dm = DistanceMatrix(get_distance_matrix(core=core, to_file=False), ssu_species)
    t = tree.nj(dm)
    print t.ascii_art()
    tips = []
    for node in t.tips():
        print node.name, node.length
        tips.append(node.name.replace(' ', '_'))
    if newick:
        n = tree.nj(dm, result_constructor=str)
        print n
    else:
        return (t, tips)

def ssu_fasta():
    with open('16s.fna', 'w+') as out_handle:
        for species in kv.get_collection('16S').distinct('species'):
            ssu = kv.get_collection('16S').find_one({'species':species})
            if ssu:
                out_handle.write(kv.make_gene_fasta(ssu, to_file=False))
            else:
                print species

def add_ssu(supp_file):
    # df = pd.read_csv(supp_file)
    # print df.columns
    # new_df = pd.DataFrame()
    # # for i in range(len(df['Strain'])):
    #     # print df['Strain'][i].replace(' ', '_').replace('.', '')

    # strain = pd.Series([df['Strain'][i].replace(' ', '_').replace('.', '') for i in range(len(df['Strain']))], name='strain')
    # ssus = df['sequences of the 16s rRNA genes']
    
    # ssu = pd.Series([ssus[i].replace(r'\n', '') if not pd.isnull(ssus[i]) else None for i in range(len(ssus))], name='16S')
    # new_df['strain'] = strain
    # new_df['16S'] = ssu

    # new_df.to_csv('ssu.csv')

    ssu_df = pd.read_csv(supp_file)
    for i in range(len(ssu_df['strain'])):
        print ssu_df['strain'][i], ssu_df['16S'][i]
        if not pd.isnull(ssu_df['16S'][i]):
            gene_record = {
                'species':ssu_df['strain'][i],
                'location':{
                    'contig':None,
                    'start':None,
                    'end':None,
                    'strand':None,
                },
                'annotation':'Small subunit ribosomal RNA',
                'dna_seq':ssu_df['16S'][i],
                'kvtag':None,
                'type':'16S'
                }

            print "adding 16S gene!"
            kv.get_collection('16S').remove({'species':ssu_df['strain'][i]})
            print kv.get_collection('16S').find_one({'species':ssu_df['strain'][i]})
            kv.get_collection('16S').insert_one(gene_record)
            print kv.get_collection('16S').find_one({'species':ssu_df['strain'][i]})

def get_list_from_tree(newick_file):
    pass


if __name__ == '__main__':
    os.chdir('/Users/KBLaptop/computation/kvasir/data/output/reorg/')
    kv.mongo_init('reorg')
