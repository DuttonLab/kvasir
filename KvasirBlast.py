#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
Must have Mongod running, in terminal: `mongod --dbpath path/to/db`
'''

def blast(mongo_db_name, blast_database):
    from pymongo import MongoClient
    from subprocess import Popen
    from Bio.Blast import NCBIXML
    from Bio.Blast.Applications import NcbiblastnCommandline
    import os
    
    client = MongoClient()
    db = client[mongo_db_name]

    all_species = db.collection_names(False)

    for species in all_species:
        current_species_collection = db[species]
        query_fasta = './tmp/{0}.fna'.format(species)

        for gene in current_species_collection.find():
            with open(query_fasta, 'w+') as query_handle:
                query_handle.write('>{0}|{1}\n{2}\n'.format(
                    gene['_id'],
                    gene['locus_tag'],
                    gene['dna_seq'],
                    )
                )

        blast_handle = NcbiblastnCommandline(
            query=query_fasta,
            db=mongo_db_name,
            perc_identity=99,
            outfmt=5,
            out="./tmp/blast_out_tmp.xml",
            max_hsps=20
            )
        print blast_handle
        stdout, stderr = blast_handle()

        with open('./tmp/blast_out_tmp.xml', 'r') as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                print blast_record
                for alignment in blast_record.alignments:
                    print alignment
                    for hsp in alignment.hsps:
                        print 'e value: ' + str(hsp.expect)
                        print(hsp.query[0:75] + '...')
                        print(hsp.match[0:75] + '...')
                        print(hsp.sbjct[0:75] + '...')  

        os.remove(query_fasta)
        os.remove('./tmp/blast_out_tmp.xml')
            
#for testing
#kvasir_blast('pipe_test', 'pipe_test')



