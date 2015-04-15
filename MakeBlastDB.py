#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
This script is designed to read a mongo database created with DataImport,
create a faa file, build a BLAST database with that file, then delete the 
faa file
'''

def make_blast_db(mongo_db_name, path_to_database):
    from pymongo import MongoClient
    from subprocess import Popen
    import os
 
    # Open and connect to Mongod server. Also opens database (which should have
    # have been created with DataImport.py)
    mongod = Popen(
        ["mongod", "--dbpath", os.path.expanduser(path_to_database)],
    )

    client = MongoClient()
    db = client[mongo_db_name]

    # Reads database and makes list of all collections (representing species)
    all_species = db.collection_names(False)
    # Handle for temporary .faa file that will contain all CDS for all species
    output_faa = '{0}.faa'.format(mongo_db_name)

    # For each collection (species) in the database, reads each gene record and
    # appends the gene and its aa sequence in FASTA format. The .faa file will 
    # contain records for all species stored in the database
    with open(output_faa, 'w+') as output_handle:
        for species in all_species:
            current_species_collection = db[species]
            for gene in current_species_collection.find():
                output_handle.write('>lcl|{0}|{1}_{2} {3}\n{4}\n'.format(
                    mongo_db_name
                    gene['species'],
                    gene['locus_tag'],
                    gene['annotation'],
                    gene['translation'],
                    )
                )
    
        # calls makeblastdb €from shell
        Popen(
            ['makeblastdb',
            '-in', output_faa,
            '-dbtype', 'prot',
            '-out', mongo_db_name,
            '-title', mongo_db_name,
            '-parse_seqids'
            ]
        ).wait() # waits for this operation to terminate before moving on

    # Removes temporary .faa file
    os.remove(output_faa)

    # Always have to close mongod
    mongod.terminate()

if __name__ == '__main__':
    import sys
    import_data(sys.argv[1], sys.argv[2])
