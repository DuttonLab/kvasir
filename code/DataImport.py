#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
This script is designed to import data from a genbank file and
write it to a mongoDB database. Must have Mongod running, in terminal:
`mongod --dbpath path/to/db`.
'''

import re, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import KvDataStructures as kv


def import_file(some_genbank):
    with open(some_genbank, 'r') as open_file:
        current_species = get_species_name(some_genbank)
        species_collection = kv.get_collection(current_species)

        print 'Working on importing {0}'.format(current_species)

        # Each "record" in genbank file is read, corresponds to individual contigs
        for record in SeqIO.parse(open_file, 'gb'):
            current_contig = get_contig(record.name)
            print "Importing {}".format(current_contig)
            ssu_gene = get_16S(record)
            if ssu_gene:
                parsed_location = kv.get_gene_location(ssu_gene[0].location)
                gene_record = {
                    'species':current_species,
                    'location':{
                        'contig':current_contig,
                        'start':parsed_location[0],
                        'end':parsed_location[1],
                        'strand':parsed_location[2],
                    },
                    'annotation':ssu_gene[0].qualifiers['product'][0],
                    'dna_seq':ssu_gene[1],
                    'kvtag':ssu_gene[0].qualifiers['kvtag'][0],
                    'type':'16S'
                    }
                print "adding 16S gene!"
                species_collection.insert_one(gene_record)
                kv.get_collection('16S').insert_one(gene_record)

            for feature in record.features:
                if feature.type == 'CDS':
                    parsed_location = kv.get_gene_location(feature.location)
                    gene_record = {
                        'species':current_species,
                        'location':{
                            'contig':current_contig,
                            'start':parsed_location[0],
                            'end':parsed_location[1],
                            'strand':parsed_location[2],
                            'index':None
                        },
                        'kvtag':feature.qualifiers['kvtag'][0],
                        'annotation':feature.qualifiers['product'][0],
                        'dna_seq':get_dna_seq(feature, record),
                        'aa_seq':feature.qualifiers['translation'][0],
                        'type':'gene'
                        }
                    species_collection.insert_one(gene_record)

def get_dna_seq(feature, record):
    if feature.location.strand == 1:
        dna_seq = record.seq[
            feature.location.start:feature.location.end
        ]
    elif feature.location.strand == -1:
        dna_seq = record.seq[
            feature.location.start:feature.location.end
        ].reverse_complement()
    return str(dna_seq)

def get_16S(gbk_record):
    for feature in gbk_record.features:
        if feature.type == 'rRNA':
            if re.search(r"16[sS]|ssu|SSU", str(feature.qualifiers['product'])):
                if len(feature) > 1000:
                    return (feature, get_dna_seq(feature, gbk_record))
                else:
                    return None

def get_species_name(path_to_genbank):
    import re
    name = re.search(r"validated_gbk/(.+)_validated\.gb", path_to_genbank)
    return name.group(1)

# Need to fix search! Only returns "contig"...
def get_contig(record_name):
    import re
    parse_contig = re.search(r'kvc_(\d\d\d)', record_name)
    return parse_contig.group(1)

def import_folder(genbank_folder):
    import os
    print 'we\'re importing a folder now!'
    for a_file in os.listdir(genbank_folder):
        current_file = os.path.join(genbank_folder, a_file)
        print 'importing {0}!'.format(current_file)
        import_file(current_file)

def import_16S(some_genbank, make_fasta=False):
    with open(some_genbank, 'r') as open_file:
        current_species = get_species_name(some_genbank)
        species_collection = kv.get_collection(current_species)
        all_16S = []
        print 'Working on importing 16s from {0}'.format(current_species)

        # Each "record" in genbank file is read, corresponds to individual contigs
        for record in SeqIO.parse(open_file, 'gb'):
            current_contig = get_contig(record.name)
            print "Importing {}".format(current_contig)
            ssu_gene = get_16S(record)
            if ssu_gene:
                gene_record = {
                    'species':current_species,
                    'contig':current_contig,
                    'location':str(ssu_gene[0].location),
                    'annotation':ssu_gene[0].qualifiers['product'][0],
                    'dna_seq':ssu_gene[1],
                    'type':'16S'
                    }
                print "adding 16S gene!"
                species_collection.insert_one(gene_record)
                kv.get_collection('16S').insert_one(gene_record)

                if make_fasta:
                    all_16S.append(gene_record)
        if make_fasta:
            return all_16S

# kv.mongo_init('permissive')
# input_dir = '/Users/KBLaptop/computation/kvasir/data/output/permissive/validated_gbk'
    
# for gb in os.listdir(input_dir):
#     if gb.endswith('.gb'):
#         import_16S(os.path.join(input_dir, gb))

#For testing:
#kv.mongo_init('scratch')
#import_file('/Users/KBLaptop/googleDrive/work/duttonLab/working_database/Microbacterium_gubbenese_published.gb')

#if __name__ == '__main__':
#    import sys
#    import os
#    kv.mongo_init(sys.argv[2])
#    if os.path.isdir(sys.argv[1]) is True:
#        print 'Looks like {0} is a directory!'.format(sys.argv[1])
#        import_folder(sys.argv[1])
#    elif os.path.isfile(sys.argv[1]) is True:
#        import_file(sys.argv[1])
#    else:
#        print 'Something\'s really wrong here...'