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

def blast_vs_fasta(subject, query):
    """
    Blast `subject` against `query`. Both must be paths to fasta file
    Returns list of lists, each `[sseqid, qseqid, pident, length]`   
    """
    out = Popen(
        ['blastn',
        '-query', subject,
        '-subject', query,
        '-outfmt', '10 sseqid qseqid pident length',
        ], stdout=PIPE
    ).communicate()[0]
    if out:
        result = re.split(r'[,\n]', out.rstrip())[0]
        return [result[i:i+4] for i in range(len(result))[0::4]]

def blast_vs_db(subject, db):
    """
    Blast `subject` (fasta file) against `db` (blast db). 
    Returns list of lists, each `[sseqid, qseqid, pident, length]`   
    """
    out = Popen(
        ['blastn',
        '-query', subject,
        '-db', db,
        '-outfmt', '10 sseqid qseqid pident length',
        ], stdout=PIPE
    ).communicate()[0]
    if out:
        result = re.split(r'[,\n]', out.rstrip())[0]
        return [result[i:i+4] for i in range(len(result))[0::4]]

def core_hgt_blast(perc_identity='99'):
    """
    Blasts all core genomes against core db
    - Set `perc_identity` if desired (default = 99)
    """
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
            '-out', 'blast_results/core/{}_blast.xml'.format(species),
            '-perc_identity', perc_identity
            ],
            stdout=PIPE
        ).communicate()[0]

        os.remove(query_fasta)

def core_blast_to_db():
    blast_dir = 'blast_results/core/'
    for f in os.listdir(blast_dir):
        if f.endswith('.xml'):
            file_handle = 'blast_results/core/{}'.format(f)
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
                hits_collection.insert_one(
                    {
                    'species':query_name,
                    'core_hits':hits_dict
                    }
                ) 

def hits_reset():
    kv.remove_collection('hits')

if __name__ == '__main__':
    import sys
    kv.mongo_init('reorg')
    kv.mongo_init(sys.argv[1])
    os.chdir('output/{}/'.format(sys.argv[1]))
    make_blast_db('core')
    make_blast_db('other')
    hits_reset()
    core_hgt_blast()
    core_blast_to_db()

