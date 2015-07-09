#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

import os, re
from Bio import SeqIO
from Bio.Seq import Seq

def add_locus_tag(some_genbank):
    with open(some_genbank, 'r') as open_file:
        
        lt_counter = 0
        genome = SeqIO.parse(open_file, 'gb')
        new_genome = []
        for record in genome:
            record.name = record.name[0:15]
            for feature in record.features:
                if feature.type == 'CDS':
                    if 'locus_tag' in feature.qualifiers:
                        print 'locus_tag {0} is good!'.format(feature.qualifiers['locus_tag'])
                        pass
                    else:
                        lt_counter += 1
                        feature.qualifiers['locus_tag'] = 'kb_{0}'.format(str(lt_counter).zfill(5))
                        'No locus_tag found for CDS, adding {0}'.format(feature.qualifiers['locus_tag'])
            new_genome.append(record)

        with open(os.path.abspath('kvasir/{0}_validated.gbk'.format(some_genbank[:-3])), 'w+') as output_handle:
            SeqIO.write(new_genome, output_handle, 'gb')
            return os.path.abspath('kvasir/{0}_validated.gbk'.format(some_genbank[:-3]))

def check_dupe_locus_tags(some_genbank):
    with open(some_genbank, 'r') as open_file:
        genome = SeqIO.parse(open_file, 'gb')
        
        locus_tag_list = []
        new_genome = []

        duplicate_flag = -1
        lt_counter = 0

        for record in genome:
            for feature in record.features:
                if feature.type == 'CDS':
                    tag = feature.qualifiers['locus_tag']
                    if tag in locus_tag_list:
                        print 'duplicate locus_tag detected: {0}'.format(tag)
                        
                        if duplicate_flag == -1:
                            duplicate_flag == 1
                    else:
                        locus_tag_list.append(tag)

        if duplicate_flag == -1:
            print 'No duplicate locus_tag found!'
            
#if __name__ == '__main__':
#    import sys
#    try:
#        add_locus_tag(sys.argv[1])
#        print sys.stderr
#for testing
#add_locus_tag("/Users/KBLaptop/googleDrive/work/Dutton Lab shared documents/Projects/HGT/genome_assemblies/RAST_annotated/1663.18_Arthrobacter.gbk")
#check_dupe_locus_tags("/Users/KBLaptop/googleDrive/work/Dutton Lab shared documents/Projects/HGT/genome_assemblies/RAST_annotated/1663.18_Arthrobacter.gbk")