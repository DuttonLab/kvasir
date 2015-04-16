#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
This script is designed to import data from a genbank file and
write it to a mongoDB database. Must have Mongod running, in terminal:
`mongod --dbpath path/to/db`
'''

import DataTypes as Data
from Bio import SeqIO
import re
from subprocess import Popen
import os

# Initial function to call. Opens and parses genbank file and creates Species objects
# (see DataTypes.py) filled with Gene objects.
def import_data(some_genbank, mongo_db_name):
    with open(some_genbank, 'r') as open_file:
        current_species = Data.Species(get_species_name(some_genbank))
        print 'Working on importing {0}'.format(current_species.species_name)

        # Each "record" in genbank file is read, corresponds to individual contigs
        for record in SeqIO.parse(open_file, 'gb'):
            current_species.add_contig(record.description)

            # Within each record may be a number of CDS's annotated, this section
            # pulls their descriptions and builds Gene objects        
            for feature in record.features:
                if feature.type == 'CDS':
                    cds = Data.Gene(
                        feature.qualifiers['locus_tag'][0],
                        feature.qualifiers['product'][0],
                        feature.qualifiers['translation'][0],
                        )
                    cds.add_location(feature.location)
                    print 'Adding {0} to import'.format(cds.annotation)
                    
                    # Adds gene objects to Species object created above
                    current_species.add_gene(cds)
        
        # Calls function (below) to build and populate MongoDB
        write_database(current_species, mongo_db_name)
            
def get_species_name(path_to_genbank):
    name = re.search(r"(\w+)\.gb", path_to_genbank)
    return name.group(1)

# Takes Species object built in import_data and adds it to MongoDB. Each database
# may take multiple Species collections, which will be filled with individual gene
# records.
def write_database(species, mongo_db_name):
    from pymongo import MongoClient
    
    # Connects to mongod server and connects to database specified in import_data.
    # Collections are documents that contain individual species. This will append
    # data if species collection already exists
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

#For testing:
#import_data('/Users/KBLaptop/computation/hgt/seqs/genomes/haloFixed.gb', 'location_import_test')

if __name__ == '__main__':
    import sys
    import_data(sys.argv[1], sys.argv[2])