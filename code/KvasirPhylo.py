#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# Unless otherwise indicated, licensed under GNU Public License (GPLv3)

import os

import KvDataStructures as kv
import pandas as pd

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
    query = StripedSmithWaterman(seq_1)
    alignment = query(seq_2)
    q = DNA(alignment.aligned_query_sequence)
    t = DNA(alignment.aligned_target_sequence)
    return q.distance(t)

def get_perc_id(seq_1, seq_2):
    return 1.0 - get_gene_distance(seq_1, seq_2)

def get_16S_distance(species_1, species_2):  
    """
    Returns hamming distance of the 16S seq of two species in MongoDB
    """
    s1_ssu = str(kv.db['16S'].find_one({'species':species_1})['dna_seq'])
    s2_ssu = str(kv.db['16S'].find_one({'species':species_2})['dna_seq'])
    return get_gene_distance(s1_ssu, s2_ssu)

def distance_matrix(core=False, to_file=True):
    all_species = kv.get_collection('core').distinct('species')
    if core:
        pass
    else:
        all_species.extend(kv.get_collection('other').distinct('species'))

    distance_matrix = pd.DataFrame(data={n:0.0 for n in all_species}, index=all_species)
    
    for pair in combinations_with_replacement(all_species, 2):
        distance = get_16S_distance(pair[0], pair[1])
        distance_matrix[pair[0]][pair[1]] = distance
        distance_matrix[pair[1]][pair[0]] = distance

    if to_file:
        distance_matrix.to_csv('distance_matrix.csv')
    else:
        return distance_matrix

def get_tree(core=False, newick=False):
    all_species = kv.get_collection('core').distinct('species')
    if core:
        pass
    else:
        all_species.extend(kv.get_collection('other').distinct('species'))
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

if __name__ == '__main__':
    seq_1 = 'AACAT'
    seq_2 = 'AAGAT'
    print get_gene_distance(seq_1, seq_2)
    print get_perc_id(seq_1, seq_2)