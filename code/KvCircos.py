#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# Unless otherwise indicated, licensed under GNU Public License (GPLv3)

import os
import re
from Bio import SeqIO, SeqUtils
from KvasirHGT import core_hgt_groups
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
        contigs = []
        sp_name = None
        for record in SeqIO.parse(file_handle, 'gb'):
            sp_name = kv.parse_genbank_name(some_gbk)
            contigs.append((record.name, len(record)))
        sp_strain = sp_name[2]
        if not os.path.isdir('circos/karyotypes/'):
            os.makedirs('circos/karyotypes/')
        with open('circos/karyotypes/karyotype_{}.txt'.format(os.path.basename(some_gbk)[:-13]), 'w+') as karyotype:
            color = [np.random.randint(0,255), np.random.randint(0,255), np.random.randint(0,255)]
            for contig in contigs:
                if contig[1] > 1000:
                    karyotype.write('chr - {0}{1} {2} {3} {4} {5},{6},{7}\n'.format(
                        sp_strain,
                        contig[0],
                        contig[0],
                        '1',
                        contig[1],
                        *color
                        )
                    )
                else:
                    break

def get_links(group=None, perc_identity='99'):
    hits_collection = kv.get_collection('hits')
    group_hits = None
    if not os.path.isdir('circos/links/'):
            os.makedirs('circos/links/')
    out_name = 'circos/links/all_links_{}.txt'.format(perc_identity)
    if group:
        groups = core_hgt_groups()
        group_hits = sorted(groups, key=len, reverse=True)[group - 1]
        out_name = 'circos/links/group{}_links_{}.txt'.format(group, perc_identity)
    
    with open(out_name, 'w+') as out_handle:
        for species in hits_collection.find():
            print species
            try:
                all_hits = species['core_hits_{}'.format(perc_identity)]
                hits_to_write = None
                if group:
                    hits_to_write = {gene:all_hits[gene] for gene in all_hits if (species['species'], gene) in group_hits}
                else:
                    hits_to_write = all_hits
                for gene in hits_to_write:
                    if hits_to_write[gene]:
                        s1_record = kv.get_mongo_record(species['species'], gene)
                        s1_strain = kv.parse_species_name(species['species'])
                        for hit in hits_to_write[gene]:
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
            except KeyError:
                pass


def get_gc(some_gbk):
    """
    Get GC content for contigs in genbank file, output file for circs line plot

    """    
    with open (some_gbk, 'r') as file_handle:
        gc_points = []
        sp_name = None
        for record in SeqIO.parse(file_handle, 'gb'):
            for b in range(len(record))[500::1000]:
                gc_cont = SeqUtils.GC(record.seq[b-500:b+499])
                gc_points.append((record.name, b-500, b+499, gc_cont))
            
            sp_name = kv.parse_genbank_name(some_gbk)
        sp_strain = sp_name[2]
        with open('circos/GC/gc_{}.txt'.format(os.path.basename(some_gbk)[:-13]), 'w+') as out_handle:
            values = [x[3] for x in gc_points]
            stats = (min(values), max(values), np.average(values), np.std(values))
            out_handle.write("# Min: {}\n# Max: {}\n# Avg, Std: {}, {}\n".format(
                stats[0], stats[1], stats[2], stats[3]
                )
            )
            for point in gc_points:
                out_handle.write('{0}{1} {2} {3} {4}\n'.format(
                    sp_strain,
                    point[0],
                    point[1],
                    point[2],
                    point[3]
                    )
                )
            return stats

def gc_circos(gb_folder):
    if not os.path.isdir('circos/GC/'):
            os.makedirs('circos/GC/')
    with open('circos/GC/plots.conf', 'w+') as out_handle:
        out_handle.write('<plots>\n')
        for f in os.listdir(gb_folder):
            if f.endswith('.gb'):
                print f
                stats = get_gc('{}{}'.format(gb_folder, f))
                out_handle.write("""    <plot>
    show  = yes
    type  = line

    file  = GC/gc_{}.txt
    r1    = 0.95r
    r0    = 1.05r
    max   = {}
    min   = {}

    color     = grey
    thickness = 1

    <rules>
    <rule>
    condition    = var(value) > {}
    color = dgreen
    </rule>
    <rule>
    condition    = var(value) < {}
    color = dred
    </rule>
    </rules>
    </plot>\n
    """.format(f[:-13], stats[1], stats[0], stats[2]+stats[3], stats[2]-stats[3])
                )
        out_handle.write('</plots>\n')

if __name__ == '__main__':
    os.chdir('/Users/KBLaptop/computation/kvasir/data/output/pacbio2/')
    kv.mongo_init('pacbio2')

    # if not os.path.isdir('circos'):
    #     os.makedirs('circos')
    # for f in os.listdir('validated_gbk/'):
    #     if f.endswith('.gb'):
    #         print f
    #         get_karyotype('validated_gbk/{}'.format(f))

    get_links(perc_identity='95')
    get_links(perc_identity='90')

    # gc_circos('validated_gbk/')

    # for entry in kv.get_collection('core').distinct('species'):
    #     print 'karyotypes/karyotype_{}.txt;'.format(entry)
        
    