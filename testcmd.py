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
    #print all_species

    for species in all_species:
        current_species_collection = db[species]
        #print current_species_collection

        for gene in current_species_collection.find():
            print gene['new_field']
            current_species_collection.update_one(
                {'locus_tag':gene['locus_tag']},
                {'$set' : {"new_field":1}},
                upsert=True)

testing('import_test', '/Users/KBLaptop/googleDrive/work/Dutton Lab shared documents/Projects/HGT/genome_assemblies/RAST_annotated/kvasir/import_test')