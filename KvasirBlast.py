#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
Must have Mongod running, in terminal: `mongod --dbpath path/to/db`
'''

from __future__ import print_function
from pymongo import MongoClient
from subprocess import Popen
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
import os

def kvasir_blast(mongo_db_name, blast_database):

    client = MongoClient()
    db = client[mongo_db_name]

    all_species = db.collection_names(False)
    print(all_species)

    for species in all_species:
        current_species_collection = db[species]
        print(current_species_collection)
        output_faa = '{0}.faa'.format(mongo_db_name)

        for gene in current_species_collection.find():
            with open(output_faa, 'w+') as output_handle:
                output_handle.write('{0}\n{1}\n'.format(
                    gene['locus_tag'],
                    gene['translation'],
                    )
                )

            blast_handle = NcbiblastpCommandline(
                query=output_faa,
                db=mongo_db_name,
                evalue='1e-50',
                outfmt=5,
                out="blast_out_tmp.xml",
                max_hsps=20
                )
            print(blast_handle)
            stdout, stderr = blast_handle()

            #blast_records = NCBIXML.parse('blast_out_tmp.xml')
            #print('hi there!', blast_records)
            #for blast_record in blast_records:
            #    for alignment in blast_record.alignments:
            #        for hsp in alignment.hsps:
            #            if hsp.positives / alignment.length > 0.9:
            #                print('****Alignment****')
            #                print('sequence:', alignment.title)
            #                print('length:', alignment.length)
            #                print('e value:', hsp.expect)
            #                print(hsp.query[0:75] + '...')
            #                print(hsp.match[0:75] + '...')
            #                print(hsp.sbjct[0:75] + '...')
            #            else:
            #                print('That one didn\'t match!')

            #os.remove('blast_out_tmp.xml')
            os.remove(output_faa)


def dont_use_this():
    blast_handle = NcbiblastpCommandline(
        query="/Users/KBLaptop/computation/hgt/seqs/genomes/brachy_IMG.faa",
        db='mongo_test_again',
        evalue=0.001,
        outfmt=5,
        out="test.xml"
        )
    print(blast_handle)
    stdout, stderr = blast_handle()

    result_handle = open("test.xml")

    blast_records = NCBIXML.parse(result_handle)

    for blast_record in blast_records:
        print(blast_record.query)
        print(blast_record.database)

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.positives / alignment.length > 0.9:
                    print('****Alignment****')
                    print('sequence:', alignment.title)
                    print('length:', alignment.length)
                    print('e value:', hsp.expect)
                    print(hsp.query[0:75] + '...')
                    print(hsp.match[0:75] + '...')
                    print(hsp.sbjct[0:75] + '...')


