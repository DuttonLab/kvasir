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
    from bson.objectid import ObjectId
    import os
    import re

    client = MongoClient()
    db = client[mongo_db_name]

    all_species = db.collection_names(False)

    for species in all_species:
        print 'Blasting {0}'.format(species)
        current_species_collection = db[species]
        query_fasta = 'kvasir/{0}.fna'.format(species)

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
            db='kvasir/{0}'.format(mongo_db_name),
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

                db_query = current_species_collection.find_one({'_id':ObjectId(query_id)})
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
                    {'_id':ObjectId(query_id)},
                    {'$set':{'hits':hits_list}},
                    upsert=True
                    )

                            #print 'Query: {0}_{1}\nHit:{2}_{3}'.format(
                            #                                        db_query['species'],
                            #                                        db_query['locus_tag'],
                            #                                        hit_species, 
                            #                                        hit_tag
                            #                                        )
                            #for hsp in alignment.hsps:
                            #    pass
                                #print hsp.query
                                #print 'e value: ' + str(hsp.expect)
                                #print(hsp.query[0:75] + '...')
                                #print(hsp.match[0:75] + '...')
                                #print(hsp.sbjct[0:75] + '...')  

            os.remove(query_fasta)
            os.remove('kvasir/blast_out_tmp.xml')
            

if __name__ == '__main__':
    import sys
    blast(sys.argv[1], sys.argv[1])

