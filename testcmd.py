#import os
#import sys
#
#print 'Here we go!'
#gbk_folder = os.path.abspath(sys.argv[1])
#
#print 'Checking files...'
#for the_file in os.listdir(gbk_folder):
#    if not the_file.startswith('.'):
#        if not os.path.isdir(the_file):
#            print the_file

def testing(mongo_db_name, blast_database):
    from pymongo import MongoClient

    client = MongoClient()
    db = client[mongo_db_name]

    all_species = db.collection_names(False)
    print all_species

    for species in all_species:
        current_species_collection = db[species]
        print current_species_collection
        output_faa = './tmp/{0}.faa'.format(mongo_db_name)

        for gene in current_species_collection.find():
            print gene['locus_tag']
            print gene
            

testing('pipe_test', 'pipe_test')