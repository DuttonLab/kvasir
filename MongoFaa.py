#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

from pymongo import MongoClient

def make_fasta_from_mongodb(mongo_db_name):
    client = MongoClient()
    db = client[mongo_db_name]
    all_species = db.collection_names(False)

    for species in all_species:
        current_species_collection = db[species]
        with open('./tmp/{0}'.format(species), 'w+') as output_handle:

            for gene in current_species_collection.find():
                output_handle.write('>gnl|{0}|{1}_{2}| {3}\n{4}\n'.format(
                            mongo_db_name,
                            gene['species'],
                            gene['locus_tag'],
                            gene['annotation'],
                            gene['translation'],
                            )
                        )

make_fasta_from_mongodb('refactor_test')