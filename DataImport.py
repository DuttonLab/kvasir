#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''This script is designed to import data from a genbank file and
write it to a mongoDB database.'''

import DataTypes as Data
from Bio import SeqIO
import re
from subprocess import Popen
import os

def import_data(some_genbank, mongo_db_name, path_to_database):
    with open(some_genbank, 'r') as open_file:
        current_species = Data.Species(get_species_name(some_genbank))
        print 'Working on importing {0}'.format(current_species.species_name)

        for record in SeqIO.parse(open_file, 'gb'):
            current_species.add_contig(record.description)
                    
            for feature in record.features:
                if feature.type == 'CDS':
                    cds = Data.Gene(
                        feature.qualifiers['locus_tag'][0],
                        feature.qualifiers['product'][0],
                        feature.qualifiers['translation'][0],
                        )
                    print 'Adding {0} to import'.format(cds.annotation)
                    current_species.add_gene(cds)
        
        write_database(current_species, mongo_db_name, path_to_database)
            
def get_species_name(path_to_genbank):
    name = re.search(r"(\w+)\.gb", path_to_genbank)
    return name.group(1)

def write_database(species, mongo_db_name, path_to_database):
    from pymongo import MongoClient
    print 'opening mongod'
    mongod = Popen(
        ["mongod", "--dbpath", os.path.expanduser(path_to_database)],
    )

    client = MongoClient()
    db = client[mongo_db_name]
    
    species_collection = db[species.species_name]

    for gene in species.genes:
        gene_record = {
            'species':species.species_name,
            'locus_tag':gene.locus_tag,
            'annotation':gene.annotation,
            'translation':gene.aa_seq,
        }
        print 'adding {0} to {1} database'.format(gene_record['annotation'], gene_record['species'])
        species_collection.insert_one(gene_record)
    print 'closing mongod'
    mongod.terminate()

#For testing:
import_data('/Users/KBLaptop/computation/hgt/seqs/genomes/haloFixed.gb', 'mongo_test_again', '~/computation/db/')

#if __name__ == "__main__":
#    import sys
#    try:
#        main(sys.argv[1], sys.argv[2], sys.argv[3])
#    else:
#        print 'Well that didn\'t work... Make you you enter a genbank file, a database name and a path to the database'