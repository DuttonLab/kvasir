#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# Unless otherwise indicated, licensed under GNU Public License (GPLv3)

'''
Must have Mongod running, in terminal: `mongod --dbpath path/to/db`
'''

import os, re
from subprocess import Popen, PIPE
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from bson.objectid import ObjectId
from pymongo.cursor import Cursor
from KvasirHGT import core_hgt_groups
import KvDataStructures as kv

def make_blast_db(source, name=None, remove_source=True):
    """
    Produces BLAST database from `source`
    Optional - provide name (defaults to `source`)
    Set remove_source=False to keep fasta file (if created)
    
    Source types:
    - fasta file (use path, must end with `.fna`)
    - Mongo collection (use name of collection)
    - list of dicts containing at least keys `species`, `_id`, `dna_seq`
    - Mongo cursor eg. `collection.find({'key':value})`
    """
    # If there's no directory for blast db's, create one
    if not os.path.isdir('blast_databases/'):
        os.makedirs('blast_databases')
    
    output_fasta = None
    
    if os.path.isfile(source):
        # Input is fasta file?
        if source.endswith('.fna'):
            output_fasta = source
            if not name:
                name = os.path.basename(source)[:-4]
            remove_source = False
        else:
            print "Not a valid file type, use .fna"

    else:
        output_fasta = '{0}_all.fasta'.format(kv.db.name)     
        genes = None
        with open(output_fasta, 'w+') as output_handle:
            if source in kv.get_collections():
                genes = kv.get_collection(source).find()
                if not name:
                    name = source
            elif type(source) == list:
                genes = source
            elif type(source) == Cursor:
                genes = source
        
            for gene in genes:
                output_handle.write('>{0}|{1}\n{2}\n'.format(
                    gene['species'],
                    gene['_id'],
                    gene['dna_seq'],
                    )
                )

    while not name:
        name = str(raw_input("enter name for BLAST database: "))

    # calls makeblastdb from shell
    print "making a database!"
    Popen(
        ['makeblastdb',
        '-in', output_fasta,
        '-dbtype', 'nucl',
        '-out', 'blast_databases/{0}'.format(name),
        '-title', name,
        ]
    ).wait() # waits for this operation to terminate before moving on

    if remove_source:
        os.remove(output_fasta)

def blast_vs_fasta(query, subject):
    """
    Blast `query` against `subject`. Both must be paths to fasta file
    Returns list of lists, each `[sseqid, qseqid, pident, length]`   
    """
    out = Popen(
        ['blastn',
        '-query', query,
        '-subject', subject,
        '-outfmt', '10 sseqid qseqid pident length',
        ], stdout=PIPE
    ).communicate()[0]
    if out:
        result = re.split(r'[,\n]', out.rstrip())[0]
        return [result[i:i+4] for i in range(len(result))[0::4]]

def blast_vs_db(query, db):
    """
    Blast `subject` (fasta file) against `db` (blast db). 
    Returns list of lists, each `[qseqid, sseqid, pident, length]`   
    """
    out = Popen(
        ['blastn',
        '-query', query,
        '-db', db,
        '-outfmt', '10 qseqid sseqid pident length',
        ], stdout=PIPE
    ).communicate()[0]
    if out:
        result = re.split(r'[,\n]', out.rstrip())
        return [result[i:i+4] for i in range(len(result))[0::4]]

def core_hgt_blast(perc_identity='99'):
    """
    Blasts all core genomes against core db
    - Set `perc_identity` if desired (default = 99)
    """
    if not os.path.isdir('blast_results/core/'):
        os.makedirs('blast_results/core/')
    for species in kv.get_collection('core').distinct('species'):
        query_fasta = 'blast_results/core/{}_tmp.fna'.format(species)
        
        with open(query_fasta, 'w+') as query_handle:
            for query in kv.get_collection('core').find({'species':species}):
                if query['type'] == 'gene':
                    query_handle.write('>{0}|{1}\n{2}\n'.format(
                        query['species'],
                        query['_id'],
                        query['dna_seq']
                        )
                    )
        
        print 'Blasting {0}'.format(species)
        out = Popen(
            ['blastn',
            '-query', query_fasta,
            '-db', 'blast_databases/core',
            '-outfmt', '5',
            '-out', 'blast_results/core/{}_{}_blast.xml'.format(species, perc_identity),
            '-perc_identity', perc_identity
            ],
            stdout=PIPE
        ).communicate()[0]

        os.remove(query_fasta)

def blast_to_db(db='core', perc_identity='99'):
    blast_dir = 'blast_results/{}/'.format(db)
    for f in os.listdir(blast_dir):
        if f.endswith('{}_blast.xml'.format(perc_identity)):
            file_handle = 'blast_results/{}/{}'.format(db,f)
            with open(file_handle, 'r') as result_handle:
                blast_records = NCBIXML.parse(result_handle)
                hits_dict = {}
                for blast_record in blast_records:
                    query_parse = re.search(r'(\w+)\|(\w+)', blast_record.query)
                    query_genus_parse = re.match(r'([A-Za-z]+)_', blast_record.query)
                    query_genus = query_genus_parse.group(1)
                    query_name = query_parse.group(1)
                    query_id = query_parse.group(2)

                    hits_dict[query_id] = []

                    for alignment in blast_record.alignments:
                        hit_parse = re.search(r'(\w+)\|(\w+)', alignment.hit_def)
                        hit_genus_parse = re.match(r'([A-Za-z]+)_', alignment.hit_def)
                        hit_genus = hit_genus_parse.group(1)

                        hit_name = hit_parse.group(1)
                        hit_id = hit_parse.group(2)
                        if query_name == hit_name:
                            pass
                        elif query_genus == hit_genus:
                            print "Oops! {} and {} are the same genus, skipping...".format(query_name, hit_name)
                            pass
                        elif kv.get_mongo_record(hit_name, hit_id)['type'] == '16S':
                            print 'Skipping 16S hit'
                        else:
                            print '=======\nhit for {0} detected:\nspecies: {1}\n======='.format(query_name, hit_name)
                            hits_dict[query_id].append((hit_name, hit_id))
                    
                print 'Updataing mongoDB with hits'
                hits_collection = kv.get_collection('hits')
                hits_collection.update_one(
                    {'species':query_name},
                    {'$set':{'{}_hits_{}'.format(db, perc_identity):{x:hits_dict[x] for x in hits_dict if hits_dict[x]}}},
                    upsert=True
                    ) 

def hits_reset():
    kv.remove_collection('hits')

def other_blast():
    groups_list = core_hgt_groups()
    groups_list.sort(key=len, reverse=True)
    
    for i in range(len(groups_list)):
        group_hits = []
        kv.make_id_list_fasta(groups_list[i], 'core')
        results = blast_vs_db('tmp.fna', 'blast_databases/other')

        hits_collection = kv.get_collection('hits')
        if results:
            for j in range(len(results)):
                group_hits.append(results[j])
        hits_collection.insert_one({'group':(i+1), 'group_hits':group_hits})


if __name__ == '__main__':
    import sys
    kv.mongo_init('pacbio2')
    os.chdir('/Users/KBLaptop/computation/kvasir/data/output/pacbio2/')
    # kv.mongo_init(sys.argv[1])
    # os.chdir('output/{}/'.format(sys.argv[1]))
    make_blast_db('core')
    make_blast_db('other')
    hits_reset()
    hgt_blast(perc_identity='90')
    hgt_blast(perc_identity='95')
    hgt_blast(perc_identity='99')
    blast_to_db(perc_identity='90')
    blast_to_db(perc_identity='95')
    blast_to_db(perc_identity='99')
    
    other_blast()
