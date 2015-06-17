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
    import pymongo
    from bson.objectid import ObjectId

    client = pymongo.MongoClient()
    db = client[mongo_db_name]

    all_species = db.collection_names(False)
    #print all_species

    for species in all_species:
        current_species_collection = db[species]
        print current_species_collection
        for record in current_species_collection.find():
            if record['hits']:
                print record['species']
                print record['locus_tag']
                
                for hit in record['hits']:
                    hit_db_record = current_species_collection.find_one({'_id':ObjectId(hit['hit_id'])})
                    print hit_db_record
                    print hit['hit_id']
        

        #for gene in current_species_collection.find():
        #    print gene['_id']
            #print gene['new_field']
            #current_species_collection.update_one(
            #    {'locus_tag':gene['locus_tag']},
            #    {'$set' : {"new_field":1}},
            #    upsert=True)

testing('blast_fix', '/Users/KBLaptop/googleDrive/work/Dutton Lab shared documents/Projects/HGT/genome_assemblies/RAST_annotated/kvasir/blast_fix')