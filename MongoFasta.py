#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

from pymongo import MongoClient

def make_fasta_from_mongodb(mongo_db_name, seq_type='nucl'):
    from pymongo import MongoClient
    client = MongoClient()
    db = client[mongo_db_name]
    all_species = db.collection_names(False)
    #print all_species
    for species in all_species:
        current_species_collection = db[species]
        #print current_species_collection
        with open('./tmp/{0}.fasta'.format(species), 'w+') as output_handle:

            for gene in current_species_collection.find():
                if seq_type == 'nucl':
                    seq = gene['dna_seq']
                elif seq_type == 'prot':
                    seq = seq = gene['aa_seq']
                else:
                    print 'That\'s not a valid sequence type, use "nucl" or "prot"'
                    break
               
                output_handle.write('>gnl|{0}|{1}_{2}| {3}\n{4}\n'.format(mongo_db_name,
                    gene['species'],
                    gene['locus_tag'],
                    gene['annotation'],
                    seq,
                    )
                )

#testing
#make_fasta_from_mongodb('all_genomes')

#if __name__ == '__main__':
#    import sys
#    if sys.argv[2] == 'nucl':
#        make_fasta_from_mongodb(sys.argv[1])
#    elif sys.argv[2] == 'prot'
#        make_fasta_from_mongodb(sys.argv[1], 'prot')
#    else:
#        print 'Did you use "prot" or "nucl" for your sequence type?'