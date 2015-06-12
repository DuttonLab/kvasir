#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
This script is designed to import data from a genbank file and
write it to a mongoDB database. Must have Mongod running, in terminal:
`mongod --dbpath path/to/db`
'''

def import_data(some_genbank, mongo_db_name):
    from pymongo import MongoClient
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    
    client = MongoClient()
    db = client[mongo_db_name]

    with open(some_genbank, 'r') as open_file:
        current_species = get_species_name(some_genbank)
        species_collection = db[current_species]
        
        print 'Working on importing {0} into {1}'.format(current_species, mongo_db_name)

        # Each "record" in genbank file is read, corresponds to individual contigs
        for record in SeqIO.parse(open_file, 'gb'):
            current_contig = record.description
       
            for feature in record.features:
                if feature.type == 'CDS':
                    
                    if feature.location.strand == 1:
                        dna_seq = record.seq[
                            feature.location.start:feature.location.end
                        ]
                    elif feature.location.strand == -1:
                        dna_seq = record.seq[
                            feature.location.start:feature.location.end
                        ].reverse_complement()
                    
                    gene_record = {
                        'species':current_species,
                        'contig':current_contig,
                        'location':str(feature.location),
                        'locus_tag':feature.qualifiers['locus_tag'][0],
                        'annotation':feature.qualifiers['product'][0],
                        'dna_seq':str(dna_seq),
                        'aa_seq':feature.qualifiers['translation'][0],
                        }
                    print "adding {0} to {1}".format(gene_record['locus_tag'], mongo_db_name)
                    species_collection.insert_one(gene_record)


def get_species_name(path_to_genbank):
    import re
    name = re.search(r"(\w+)\.gb", path_to_genbank)
    return name.group(1)

def import_folder(genbank_folder, mongo_db_name):
    import os
    print 'we\'re importing a folder now!'
    for a_file in os.listdir(genbank_folder):
        current_file = os.path.join(genbank_folder, a_file)
        print 'importing {0}!'.format(current_file)
        import_data(current_file, mongo_db_name)

#For testing:
#import_folder('/Users/KBLaptop/computation/genomes/', 'prot-nucl')

if __name__ == '__main__':
    import sys
    import os
    if os.path.isdir(sys.argv[1]) is True:
        print 'Looks like {0} is a directory!'.format(sys.argv[1])
        import_folder(sys.argv[1], sys.argv[2])
    elif os.path.isdir(sys.argv[1]) is False:
        import_data(sys.argv[1], sys.argv[2])
    else:
        print 'Something\'s really wrong here...'