#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

import re
import os
import KvDataStructures as kv
import Outputs as o
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
    
    Popen(
        ['makeblastdb',
            '-in', s2,
            '-dbtype', 'nucl',
            '-out', '{}_blastdb'.format(species_2.replace(' ', '_')),
            '-title', species_2,
        ]
    ).wait()
    out = Popen(
        ['blastn',
        '-query', s1,
        '-db', '{}_blastdb'.format(species_2.replace(' ', '_')),
        '-outfmt', '5',
        '-out', 'blast_results.xml',
        '-num_alignments', '1',
        ],
    ).communicate()[0]

def get_clocks(species):
    clock_genes = ["dnaG", "frr", "infC", "nusA", "pgk", "pyrG", "rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplK", "rplL", "rplM", "rplN", "rplP", "rplS", "rplT", "rpmA", "rpoB", "rpsB", "rpsC", "rpsE", "rpsI", "rpsJ", "rpsK", "rpsM", "rpsS", "smpB", "tsf"]
    current_species_collection = kv.get_collection(species)
    clock_dict = {}
    for record in current_species_collection.find():
        if any(gene.lower() in record['annotation'].lower() for gene in clock_genes):
            print record['annotation']

def all_by_all_distance():
    for current_species_collection in kv.mongo_iter():
        pass

if __name__ == '__main__':
    kv.mongo_init('dutton_cheese')
    print kv.get_species_collections()
    os.chdir('/Users/KBLaptop/computation/kvasir/data/output/{}/'.format(kv.db.name))
    fasta_blast('Arthrobacter sp JB182', 'Brevibacterium sp JB5')
