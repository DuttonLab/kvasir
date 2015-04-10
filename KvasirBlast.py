#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

def make_blast_db(mongo_db_name, path_to_database):
    from pymongo import MongoClient
    from subprocess import Popen
    import os
 
    mongod = Popen(
        ["mongod", "--dbpath", os.path.expanduser(path_to_database)],
    )

    client = MongoClient()
    db = client[mongo_db_name]
    
    collection = db['test_collection']

    all_species = db.collection_names(False)
    output_faa = '{0}.faa'.format(mongo_db_name)
    with open(output_faa, 'w+') as output_handle:
        for species in all_species:
            current_species_collection = db[species]
            for gene in current_species_collection.find():
                output_handle.write('>{0} | {1} | {2}\n{3}\n'.format(
                    gene['species'],
                    gene['locus_tag'],
                    gene['annotation'],
                    gene['translation'],
                    )
                )
    
    
        Popen(
            ['makeblastdb',
            '-in', output_faa,
            '-dbtype', 'prot',
            '-out', mongo_db_name,
            '-title', mongo_db_name
            ]
        ).wait()

        Popen(['rm', '{0}'.format(output_faa)]).wait()

    mongod.terminate()


make_blast_db('mongo_test_again', '~/computation/db/')
