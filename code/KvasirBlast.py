#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
Must have Mongod running, in terminal: `mongod --dbpath path/to/db`
'''

import os, re
from subprocess import Popen, PIPE
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from bson.objectid import ObjectId
import KvDataStructures as kv


def blast():
    for current_species_collection in kv.mongo_iter():
        query_fasta = 'output/{0}_tmp.fna'.format(current_species_collection.name)

        with open(query_fasta, 'w+') as query_handle:
            print "making {}".format(query_handle)
            for query in current_species_collection.find():
                query_handle.write('>{0}|{1}\n{2}\n'.format(
                    query['species'],
                    query['_id'],
                    query['dna_seq']
                    )
                )
        
        print 'Blasting {0}'.format(current_species_collection.name)
        out = Popen(
            ['blastn',
            '-query', query_fasta,
            '-db', 'output/{0}'.format(kv.db.name),
            '-outfmt', '5',
            '-out', 'output/blast_out_tmp.xml',
            ],
            stdout=PIPE
        ).communicate()[0]

        with open('output/blast_out_tmp.xml', 'r') as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            hits_list = []
            for blast_record in blast_records:
                query_parse = re.search(r'(\w+)([\w\s]+)\|(\w+)', blast_record.query)
                query_genus = query_parse.group(1)
                query_name = '{}{}'.format(query_parse.group(1), query_parse.group(2))
                query_id = query_parse.group(3)

                db_query = kv.get_mongo_record(current_species_collection, query_id)

                for alignment in blast_record.alignments:
                    hit_parse = re.search(r'(\w+)([\w\s]+)\|(\w+)', alignment.hit_def)
                    hit_genus = hit_parse.group(1)
                    hit_name = '{}{}'.format(hit_parse.group(1), hit_parse.group(2))
                    hit_id = hit_parse.group(3)
                    
                    if query_name == hit_name:
                        pass
                    elif query_genus == hit_genus:
                        print "Oops! {} and {} are the same genus, skipping...".format(query_name, hit_name)
                        pass
                    else:
                        print '=======\nhit for {0} detected:\nspecies: {1}\n======='.format(db_query['locus_tag'], hit_name)
                        hits_list.append({
                            'hit_species':hit_name,
                            'hit_id':hit_id,
                            })
                
            print 'Updataing mongoDB with hits'
            hits_collection = kv.get_collection('hits')
            hits_collection.insert_one(
                {
                'species':query_name,
                'hits':hits_list
                }
            ) 

        os.remove(query_fasta)
        os.remove('output/blast_out_tmp.xml')
            
def test_blast():
    for current_species_collection in kv.mongo_iter():
        query_fasta = 'output/tmp/{0}.fna'.format(current_species_collection.name)

        with open(query_fasta, 'w+') as query_handle:
            print "making {}".format(query_handle)
            for query in current_species_collection.find():
                query_handle.write('>{0}|{1}\n{2}\n'.format(
                    query['species'],
                    query['_id'],
                    query['dna_seq']
                    )
                )
        print 'Blasting {0}'.format(current_species_collection.name)
        out = Popen(
            ['blastn',
            '-query', query_fasta,
            '-db', 'output/{0}'.format(kv.db.name),
            '-outfmt', '5',
            '-out', 'output/tmp/{}_out_tmp.xml'.format(current_species_collection.name),
            ],
            stdout=PIPE
        ).communicate()[0]

if __name__ == '__main__':
    import sys
    kv.mongo_init(sys.argv[1])
    blast()

