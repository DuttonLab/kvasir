#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
Must have Mongod running, in terminal: `mongod --dbpath path/to/db`
'''

import os, re
from subprocess import Popen
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from bson.objectid import ObjectId
import KvDataStructures as kv


def blast():
    for current_species_collection in kv.mongo_iter():
        print 'Blasting {0}'.format(current_species_collection.name)
        query_fasta = 'kvasir/{0}.fna'.format(current_species_collection.name)

        with open(query_fasta, 'w+') as query_handle:
            for query in current_species_collection.find():
                query_handle.write('>{0}|{1}\n{2}\n'.format(
                    query['species'],
                    query['_id'],
                    query['dna_seq']
                    )
                )

        blast_handle = NcbiblastnCommandline(
            query=query_fasta,
            db='kvasir/{0}'.format(kv.db.name),
            perc_identity=99,
            outfmt=5,
            out="kvasir/blast_out_tmp.xml",
            max_hsps=20
            )
        #print blast_handle
        stdout, stderr = blast_handle()

        with open('kvasir/blast_out_tmp.xml', 'r') as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                query_parse = re.search(r'(\w+)\|(\w+)', blast_record.query)
                query_species = query_parse.group(1)
                query_id = query_parse.group(2)

                db_query = kv.get_mongo_record(current_species_collection, query_id)
                hits_list = []

                for alignment in blast_record.alignments:
                    hit_parse = re.search(r'(\w+)\|(\w+)', alignment.hit_def)
                    hit_species = hit_parse.group(1)
                    hit_id = hit_parse.group(2)
                    
                    if db_query['species'] == hit_species:
                        pass
                    else:
                        print 'hit for {0} detected:\nspecies: {1}'.format(db_query['locus_tag'], hit_species)
                        hits_list.append({
                            'hit_species':hit_species,
                            'hit_id':hit_id,
                            })
                
                print 'Updataing mongoDB with hits'
                current_species_collection.update_one(
                    {'_id':ObjectId(query_id), 'species':query_species},
                    {'$set':{'hits':hits_list}},
                    upsert=True
                    ) 

            os.remove(query_fasta)
            os.remove('kvasir/blast_out_tmp.xml')
            

if __name__ == '__main__':
    import sys
    blast(sys.argv[1], sys.argv[1])

