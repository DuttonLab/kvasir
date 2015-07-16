#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''
This script is designed to read a mongo database created with DataImport,
create a fasta file, build a BLAST database with that file, then delete the 
fasta file. Must have Mongod running, in terminal: `mongod --dbpath path/to/db`
'''
import os
from subprocess import Popen
import KvDataStructures as kv

def make_blast_db(seq_type='nucl'):
    # Reads database and makes list of all collections (representing species)
    # Handle for temporary .fasta file that will contain all CDS for all species
    database_name = kv.db.name
    output_fasta = 'kvasir/{0}.fasta'.format(database_name)

    # For each collection (species) in the database, reads each gene record and
    # appends the gene and its aa sequence in FASTA format. The .fasta file will 
    # contain records for all species stored in the database
   
    with open(output_fasta, 'w+') as output_handle:

        for current_species_collection in kv.mongo_iter():
            for gene in current_species_collection.find():
                if seq_type == 'nucl':
                    seq = gene['dna_seq']
                elif seq_type == 'prot':
                    seq = gene['aa_seq']
                else:
                    print 'That\'s not a valid sequence type, use "nucl" or "prot"'
                    break

                output_handle.write('>{0}|{1}\n{2}\n'.format(
                    gene['species'],
                    gene['_id'],
                    seq,
                    )
                )
    
        # calls makeblastdb from shell
        print "making a database!"
        Popen(
            ['makeblastdb',
            '-in', output_fasta,
            '-dbtype', seq_type,
            '-out', 'kvasir/{0}'.format(database_name),
            '-title', database_name,
            ]
        ).wait() # waits for this operation to terminate before moving on

    # Removes temporary .fasta file
    os.remove(output_fasta)

#for testing
#kv.mongo_init('restructured')
#make_blast_db()

if __name__ == '__main__':
    import sys
    kv.mongo_init(sys.argv[1])
    print 'collections:', kv.get_collections()
    make_blast_db()