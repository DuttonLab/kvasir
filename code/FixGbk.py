#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

import os, re
from Bio import SeqIO
from Bio.Seq import Seq

def validate_gbk(some_genbank):
    file_name = os.path.basename(some_genbank)
    new_file_name = 'validated_gbk/{0}_validated.gb'.format(file_name[:-3])
    if not os.path.isfile(new_file_name):
        with open(some_genbank, 'r') as open_file:
            
            lt_counter = 0
            contig_counter = 0
            genome = SeqIO.parse(open_file, 'gb')
            new_genome = []

            for record in genome:
                contig_counter += 1
                record.name = 'kvc_{0}_{1}'.format(str(contig_counter).zfill(3), record.name)[:15]
                for feature in record.features:
                    if feature.type == 'CDS':
                        if 'locus_tag' in feature.qualifiers:
                            print 'locus_tag {0} is good!'.format(feature.qualifiers['locus_tag'])
                            pass
                        else:
                            lt_counter += 1
                            feature.qualifiers['locus_tag'] = 'kvtag_{0}'.format(str(lt_counter).zfill(5))
                            'No locus_tag found for CDS, adding {0}'.format(feature.qualifiers['locus_tag'])
                new_genome.append(record)
            if not os.path.isdir('validated_gbk/'):
                os.makedirs('validated_gbk/')
            with open(new_file_name, 'w+') as output_handle:
                SeqIO.write(new_genome, output_handle, 'gb')
                return os.path.abspath(new_file_name)
    else:
        print 'already fixed, skipping validation!'
        return False

def has_dupe_locus_tags(parsed_gbk):
    locus_tag_list = []
    new_genome = []

    for record in parsed_gbk:
        for feature in record.features:
            if feature.type == 'CDS':
                tag = feature.qualifiers['locus_tag']
                if tag in locus_tag_list:
                    print 'duplicate locus_tag detected: {0}'.format(tag)
                    return True
                else:
                    locus_tag_list.append(tag)

    return False

#def fix_locus_tags(parsed_gbk, dupes=False):
#    if dupes:
#        print 'Providing unique locus_tags... just kidding, this function doesn\'t work!'
#        lt_counter = 0
#        break


            
#if __name__ == '__main__':
#    import sys
#    try:
#        add_locus_tag(sys.argv[1])
#        print sys.stderr
#for testing
#add_locus_tag("/Users/KBLaptop/googleDrive/work/Dutton Lab shared documents/Projects/HGT/genome_assemblies/RAST_annotated/1663.18_Arthrobacter.gbk")
#check_dupe_locus_tags("/Users/KBLaptop/googleDrive/work/Dutton Lab shared documents/Projects/HGT/genome_assemblies/RAST_annotated/1663.18_Arthrobacter.gbk")