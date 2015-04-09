#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

def make_blast_db(mongo_db_name):
    from pymongo import MongoClient
    import subprocess
    client = MongoClient()
    db = client[mongo_db_name]

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
    
    
        subprocess.call('makeblastdb -in {0} -dbtype "prot" -out {1} -title {1}'.format(
            output_faa,
            mongo_db_name,
            ), shell=True
        )

        subprocess.call('rm {0}'.format(output_faa), shell=True)


make_blast_db('demeter_test_db')
