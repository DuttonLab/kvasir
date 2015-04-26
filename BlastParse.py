#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

from Bio.Blast import NCBIXML
import re
from pymongo import MongoClient

client = MongoClient()

species_db = client['all_genomes']

hits_db = client['hits_test']
hits_collection = hits_db['arthro_hits']

with open('./tmp/arthro_all.xml', 'r') as result_handle:
    blast_records = NCBIXML.parse(result_handle)
    grep = re.compile(r'gnl\|all_genomes\|(\w+IMG|\w+Fixed2?)_(\w+)')

    for blast_record in blast_records:
        if len(blast_record.alignments) > 1:
            query_match = grep.match(str(blast_record.query))
            query_species = query_match.group(1)
            query_locus = query_match.group(2)
            query_collection = species_db[str(query_match.group(1))]
            query_entry = species_db[query_species].find({'locus_tag':query_match.group(2)})[0]
            print 'checking hits for {0} at {1}'.format(blast_record.query, query_entry['_id'])
            
            for hit in blast_record.alignments:
                hit_match = grep.match(str(hit))
                
                if hit_match is None:
                    print 'couldn\'t parse {0}'.format(str(hit))
                elif query_match.group(2) == hit_match.group(2):
                    pass
                else:
                    if hit.length < 50:
                        pass
                    else:
                        hit_collection = species_db[str(hit_match.group(1))]
                        for entry in hit_collection.find({'locus_tag':hit_match.group(2)}):
                            hits_collection.insert_one({
                                'query_id':query_entry['_id'],
                                'hit_id':entry['_id']
                            })

                            print entry['_id']
