#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

import DataTypes as Data
from Bio import SeqIO
import re

def import_data(some_genbank):
    with open(some_genbank, 'r') as open_file:
        current_species = Data.Species(get_species_name(some_genbank))
        for record in SeqIO.parse(open_file, 'gb'):
            current_species.add_contig(record.description)
            for feature in record.features:
                if feature.type == 'CDS':
                    cds = Data.Gene(
                        feature.qualifiers['locus_tag'][0],
                        feature.qualifiers['product'][0],
                        feature.qualifiers['translation'][0],
                        )
                    current_species.add_gene(cds)
        write_database(current_species)
        #for gene in current_species.genes:
            #print current_species.genes[gene].get_faa()
            
def get_species_name(path_to_genbank):
    name = re.search(r"(\w+)\.gb", path_to_genbank)
    return name.group(1)

def write_database(species):
    from pymongo import MongoClient
    db = MongoClient().demeter_test_db
    
    species_collection = db[species.species_name]

    for gene in species.genes:
        gene_record = {
            'species':species.species_name,
            'locus_tag':gene.locus_tag,
            'annotation':gene.annotation,
            'translation':gene.aa_seq,
        }
        species_collection.insert_one(gene_record)

