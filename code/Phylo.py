#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

import re
import os
import KvDataStructures as kv
import Outputs as o
from matplotlib import pyplot as plt
from itertools import combinations
from numpy import std
from Bio import SeqIO
from Bio.Blast import NCBIXML
from bson.objectid import ObjectId
from subprocess import Popen, PIPE

def make_species_fasta(species):
    with open('{}.fna'.format(species.replace(' ', '_')), 'w+') as output_handle:
        for record in kv.get_collection(species).find():
            output_handle.write(
                ">{}|{}\n{}\n".format(record['species'].replace(' ', '_'), record['_id'], record['dna_seq'])
            )
    return '{}.fna'.format(species)

def fasta_blast(species_1, species_2):
    s1 = make_species_fasta(species_1)
    s2 = make_species_fasta(species_2)
    
    # Popen(
    #     ['makeblastdb',
    #         '-in', s2,
    #         '-dbtype', 'nucl',
    #         '-out', '{}_blastdb'.format(species_2.replace(' ', '_')),
    #         '-title', species_2,
    #     ]
    # ).wait()
    out = Popen(
        ['blastn',
        '-query', s1,
        '-db', '{}_blastdb'.format(species_2.replace(' ', '_')),
        '-outfmt', '5',
        '-out', 'blast_results.xml',
        ],
    ).communicate()[0]

def get_clocks(species):
    clock_genes = ["dnaG", "frr", "infC", "nusA", "pgk", "pyrG", "rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplK", "rplL", "rplM", "rplN", "rplP", "rplS", "rplT", "rpmA", "rpoB", "rpsB", "rpsC", "rpsE", "rpsI", "rpsJ", "rpsK", "rpsM", "rpsS", "smpB", "tsf"]
    current_species_collection = kv.get_collection(species)
    clock_dict = {}
    for record in current_species_collection.find():
        if any(gene.lower() in record['annotation'].lower() for gene in clock_genes):
            print record['annotation']

def all_by_all_distance(species_1, species_2):
    global fig_counter
    fig_counter += 1
    with open('blast_results.xml', 'r') as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        distance_list = []
        for blast_record in blast_records:
            s1_ssu = None
            s2_ssu = None
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.align_length > 100:
                        pident = float(hsp.positives)/float(hsp.align_length)
                        length = hsp.align_length
                        distance_list.append((length, pident))
                    break
                break

        x, y = zip(*distance_list)
        average = sum(y) / float(len(y))
        dev = std(y)
        plt.figure(fig_counter)
        plt.scatter(x,y)
        plt.ylabel('percent identity')
        plt.xlabel('length')
        plt.xlim(xmin=0)
        plt.axhline(y=average, color='r')
        plt.axhline(y=(average+dev), color='r', linestyle=':')
        plt.axhline(y=(average-dev), color='r', linestyle=':')
        plt.title("{} and {}".format(species_1, species_2))
        plt.savefig('Figure_{}.pdf'.format(fig_counter))
        plt.close()


if __name__ == '__main__':
    import sys
    kv.mongo_init('once_again')
    os.chdir('/Users/KBLaptop/computation/kvasir/data/output/once_again/')
    ls = kv.get_species_collections()
    ls.remove('Arthrobacter_arilaitensis_Re117')
    pairs = combinations(ls, 2)
    fig_counter = 0
    for pair in pairs:
        fasta_blast(pair[0],pair[1])
        all_by_all_distance(pair[0],pair[1])
