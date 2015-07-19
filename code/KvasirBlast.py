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

def make_blast_db(seq_type='nucl'):
    # Reads database and makes list of all collections (representing species)
    # Handle for temporary .fasta file that will contain all CDS for all species
    output_fasta = 'output/{0}/{0}.fasta'.format(kv.db.name)

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
            '-out', 'output/{0}/{0}_blastdb'.format(database_name),
            '-title', database_name,
            ]
        ).wait() # waits for this operation to terminate before moving on

    # Removes temporary .fasta file
    os.remove(output_fasta)

def blast():
    for current_species_collection in kv.mongo_iter():
        query_fasta = 'output/{}/blast_results/{}_tmp.fna'.format(kv.db.name, current_species_collection.name)

        with open(query_fasta, 'w+') as query_handle:
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
            '-out', 'output/{}/blast_results/{}_blast.xml'.format(kv.db.name, current_species_collection.name),
            '-perc_identity', '99'
            ],
            stdout=PIPE
        ).communicate()[0]

        os.remove(query_fasta)
        # os.remove('output/blast_out_tmp.xml')

def blast_to_db():
    blast_dir = 'output/{}/blast_results/'.format(kv.db.name)
    for f in os.listdir(blast_dir):
        file_handle = os.path.join(blast_dir, f)
        with open(file_handle, 'r') as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            hits_dict = {}
            for blast_record in blast_records:
                query_parse = re.search(r'(\w+)([\w\s]+)\|(\w+)', blast_record.query)
                query_genus = query_parse.group(1)
                query_name = '{}{}'.format(query_parse.group(1), query_parse.group(2))
                query_id = query_parse.group(3)

                hits_dict[query_id] = []
                
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
                        print '=======\nhit for {0} detected:\nspecies: {1}\n======='.format(query_name, hit_name)
                        hits_dict[query_id].append((hit_name, hit_id))
                
            print 'Updataing mongoDB with hits'
            hits_collection = kv.get_collection('hits')
            hits_collection.insert_one(
                {
                'species':query_name,
                'hits':hits_dict
                }
            ) 

def hits_reset():
    kv.remove_collection('hits')
        
if __name__ == '__main__':
    import sys
    kv.mongo_init(sys.argv[1])
    if len(sys.argv) > 2:
        if sys.argv[2] == 'reset':
            hits_reset()
    else:
        
        blast()
        blast_to_db()

