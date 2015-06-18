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

def testing(mongo_db_name):
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
            print record
            #if record['hits']:
            #    print record['species']
            #    print record['_id']
            #    for hit in record['hits']:
            #        #hit_db_record = db[hit['hit_species']].find_one({'_id': ObjectId(hit['hit_id'])})
            #        print hit['hit_id']
        

        #for gene in current_species_collection.find():
        #    print gene['_id']
            #print gene['new_field']
            #current_species_collection.update_one(
            #    {'locus_tag':gene['locus_tag']},
            #    {'$set' : {"new_field":1}},
            #    upsert=True)

def find_contig(some_genbank):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    import re
    
    with open(some_genbank, 'r') as open_file:
        current_species = get_species_name(some_genbank)

        # Each "record" in genbank file is read, corresponds to individual contigs
        for record in SeqIO.parse(open_file, 'gb'):
            parse_contig = re.search(r'contig_\d+|NODE_\d+', record.description)
            current_contig = parse_contig.group(0)
            print current_contig

def get_species_name(path_to_genbank):
    import re
    name = re.search(r"(\w+)_\w+_\w+.gb", path_to_genbank)
    return name.group(1)

def collapse_lists(list_of_lists):
    # example input: [[1,2,3],[3,4],[5,6,7],[1,8,9,10],[11],[11,12],[13],[5,12]]
    result = []
    for d in list_of_lists:
        d = set(d)

        matched = [d]
        unmatched = []
        # first divide into matching and non-matching groups
        for g in result:
            if d & g:
                matched.append(g)
            else:
                unmatched.append(g)
        # then merge matching groups
        result = unmatched + [set().union(*matched)]
    return result

def random_code():
    x, y = (1,2)
    print x



