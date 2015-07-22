#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

clock_genes = ["dnaG", "frr", "infC", "nusA", "pgk", "pyrG", "rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplK", "rplL", "rplM", "rplN", "rplP", "rplS", "rplT", "rpmA", "rpoB", "rpsB", "rpsC", "rpsE", "rpsI", "rpsJ", "rpsK", "rpsM", "rpsS", "smpB", "tsf"]

import re
import os
import KvDataStructures as kv
from bson.objectid import ObjectId
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from DataImport import get_species_name

def make_fasta(gene):    
    with open('{}.fna'.format(gene), 'w+') as output_fna, open('{}.faa'.format(gene), 'w+') as output_faa:
        for current_species_collection in kv.mongo_iter():
            seq_found = False
            for record in current_species_collection.find():
                if re.search(gene, record['annotation'], re.IGNORECASE):
                    print record['species'], ', ', record['locus_tag']
                    print record['annotation']
                    if seq_found:
                        if len(record['dna_seq']) > len(seq_found['dna_seq']):
                            seq_found = record
                    elif len(record['dna_seq']) > 300:
                        seq_found = record
            if seq_found:
                output_fna.write(
                    '>{0}|{1}|{2}\n{3}\n'.format(
                        seq_found['species'].replace(' ', '_'),
                        seq_found['locus_tag'],
                        str(seq_found['_id']),
                        seq_found['dna_seq'],
                        )
                    )
                output_faa.write(
                        '>{0}|{1}|{2}\n{3}\n'.format(
                            record['species'].replace(' ', '_'),
                            record['locus_tag'],
                            str(record['_id']),
                            record['aa_seq'],
                            )
                        )


# kv.mongo_init('permissive')
# input_dir = '/Users/KBLaptop/computation/kvasir/data/output/permissive/validated_gbk'
    
# for gb in os.listdir(input_dir):
#     if gb.endswith('.gb'):
#         import_16S(os.path.join(input_dir, gb))



                    