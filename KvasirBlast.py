#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
Must have Mongod running, in terminal: `mongod --dbpath path/to/db`
'''
from pymongo import MongoClient
from subprocess import Popen
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import os

def kvasir_blast(mongo_db_name, blast_database):

    client = MongoClient()
    db = client[mongo_db_name]

    all_species = db.collection_names(False)
    print all_species

    for species in all_species:
        current_species_collection = db[species]
        print current_species_collection
        output_faa = './tmp/{0}.fna'.format(mongo_db_name)

        for gene in current_species_collection.find():
            with open(output_faa, 'w+') as output_handle:
                output_handle.write('{0}\n{1}\n'.format(
                    gene['locus_tag'],
                    gene['dna_seq'],
                    )
                )

            blast_handle = NcbiblastnCommandline(
                query=output_faa,
                db=mongo_db_name,
                perc_identity=99,
                outfmt=5,
                out="/tmp/blast_out_tmp.xml",
                max_hsps=20
                )
            print blast_handle
            stdout, stderr = blast_handle()

            blast_records = NCBIXML.parse('./tmp/blast_out_tmp.xml')
            print 'hi there!' 
            print blast_records
            for blast_record in blast_records:
                print blast_record
                for alignment in blast_record.alignments:
                    print alignment
                    for hsp in alignment.hsps:
                        print hsp
                        if hsp.positives / alignment.length > 0.9:
                            print('****Alignment****')
                            print('sequence:', alignment.title)
                            print('length:', alignment.length)
                            print('e value:', hsp.expect)
                            print(hsp.query[0:75] + '...')
                            print(hsp.match[0:75] + '...')
                            print(hsp.sbjct[0:75] + '...')
                        else:
                            print('That one didn\'t match!')

            os.remove(output_faa)
            os.remove('./tmp/blast_out_tmp.xml')
            
#for testing
kvasir_blast('pipe_test', 'pipe_test')



