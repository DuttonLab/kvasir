#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
This script is designed to read a mongo database created with DataImport,
create a fasta file, build a BLAST database with that file, then delete the 
fasta file. Must have Mongod running, in terminal: `mongod --dbpath path/to/db`
'''

def make_blast_db(mongo_db_name, seq_type='nucl'):
    from pymongo import MongoClient
    from subprocess import Popen
    import os

    # Previous Mongod instance should be running
    client = MongoClient()
    db = client[mongo_db_name]

    # Reads database and makes list of all collections (representing species)
    all_species = db.collection_names(False)
    # Handle for temporary .fasta file that will contain all CDS for all species
    output_fasta = 'kvasir/{0}.fasta'.format(mongo_db_name)

    # For each collection (species) in the database, reads each gene record and
    # appends the gene and its aa sequence in FASTA format. The .fasta file will 
    # contain records for all species stored in the database
   
    with open(output_fasta, 'w+') as output_handle:

        for species in all_species:
            current_species_collection = db[species]
            for gene in current_species_collection.find():
                if seq_type == 'nucl':
                    seq = gene['dna_seq']
                elif seq_type == 'prot':
                    seq = gene['aa_seq']
                else:
                    print 'That\'s not a valid sequence type, use "nucl" or "prot"'
                    break


                output_handle.write('>{0}_{1}|{2}|{3}\n{4}\n'.format(
                    mongo_db_name,
                    gene['species'],
                    gene['_id'],
                    gene['locus_tag'],
                    seq,
                    )
                )
    
        # calls makeblastdb from shell
        print "making a database!"
        Popen(
            ['makeblastdb',
            '-in', output_fasta,
            '-dbtype', seq_type,
            '-out', 'kvasir/{0}'.format(mongo_db_name),
            '-title', mongo_db_name,
            ]
        ).wait() # waits for this operation to terminate before moving on

    # Removes temporary .fasta file
    os.remove(output_fasta)

#for testing
#make_blast_db(sys.argv[1])

#if __name__ == '__main__':
#    import sys
#    make_blast_db(sys.argv[1])
