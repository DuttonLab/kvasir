#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# Unless otherwise indicated, licensed under GNU Public License (GPLv3)

import os
from Bio import SeqIO
import numpy as np
import KvDataStructures as kv
import pandas as pd

def get_karyotype(some_gbk):
    """
    Convert Genbank file into Karyotype file for Circos
    - Each contig is a "chromosome"
    - format: 'chr - ID LABEL START END COLOR'
    """
    with open (some_gbk, 'r') as file_handle:
        counter = 1
        contigs = []
        sp_name = None
        for record in SeqIO.parse(file_handle, 'gb'):
            sp_name = kv.parse_genbank_name(some_gbk)
            contigs.append((record.name, len(record)))
        sp_strain = sp_name[2]
        with open('circos/karyotypes/karyotype_{}.txt'.format(os.path.basename(some_gbk)[:-13]), 'w+') as karyotype:
            for contig in contigs:
                if contig[1] > 1000:
                    karyotype.write('chr - {0}{1} {2} {3} {4} {5},{6},{7}\n'.format(
                        sp_strain,
                        contig[0],
                        contig[0],
                        '1',
                        contig[1],
                        np.random.randint(0,255),
                        np.random.randint(0,255),
                        np.random.randint(0,255),
                        )
                    )
                    counter +=1
                else:
                    break

def get_links():
    hits_collection = kv.get_collection('hits')
    
    with open('circos/links.txt', 'w+') as out_handle:
        for species in hits_collection.find():
            all_hits = species['core_hits']
            for gene in all_hits:
                if all_hits[gene]:
                    s1_record = kv.get_mongo_record(species['species'], gene)
                    s1_strain = kv.parse_species_name(species['species'])
                    for hit in all_hits[gene]:
                        s2_record = kv.get_mongo_record(hit[0], hit[1])
                        s2_strain = kv.parse_species_name(hit[0])
                        out_handle.write('{0}kvc_{1} {2} {3} {4}kvc_{5} {6} {7}\n'.format(
                            s1_strain[2],
                            s1_record['location']['contig'],
                            s1_record['location']['start'],
                            s1_record['location']['end'],
                            s2_strain[2],
                            s2_record['location']['contig'],
                            s2_record['location']['start'],
                            s2_record['location']['end'],
                            )
                        )

def get_gc(folder):
    """
    Takes folder of csv files from Geneious GC content export
    
    - csv format eg: `Position,GC Content - kvc_001,GC Content - kvc_002,...`
    - file names eg: `JB182_gc.csv
    """
    for f in os.listdir(folder):
        strain = f[:-7]
        print strain


if __name__ == '__main__':
    os.chdir('/Users/KBLaptop/computation/kvasir/data/output/reorg/')
    kv.mongo_init('reorg')

    get_gc('circos/GC')

    # for entry in kv.get_collection('core').distinct('species'):
    #     print 'karyotypes/karyotype_{}.txt;'.format(entry)
    
    # get_links()
    
    # if not os.path.isdir('circos'):
    #     os.makedirs('circos')
    # for f in os.listdir('validated_gbk/'):
    #     if f.endswith('.gb'):
    #         print f
    #         get_karyotype('validated_gbk/{}'.format(f))
