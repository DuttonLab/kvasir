#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# Unless otherwise indicated, licensed under GNU Public License (GPLv3)

import os, re
from Bio import SeqIO

def confirm_species_name(file_name):
    user_validate = str(raw_input("Is {} the correct species name? [y/n]".format(file_name[:-3])))
    if user_validate == 'y':
        return file_name[:-3].replace(' ', '_').replace('.', '')
    elif user_validate == 'n':
        print "Please use only letters for the following"
        genus = raw_input("Please enter organism genus: ")
        species = raw_input("Please enter organism species (use \"sp\" if unknown):")
        print "Please use only letters, numbers and '_' for the following"
        strain = raw_input("Please enter organism strain (use \"undef\" if unknown):")
        if re.match(r'^[A-Za-z]+$', genus) and re.match(r'^[A-Za-z]+$', species) and re.match(r'^\w+$', strain):
            return '{}_{}_{}'.format(genus, species, strain)
        else:
            "Something's wrong, try again"
            confirm_species_name(file_name)
    else:
        print "Please enter \"y\" or \"n\"... try again"
        confirm_species_name(file_name)

def validate_gbk(some_genbank, user_validate=True):
    file_name = os.path.basename(some_genbank)
    sp_name = None
    if user_validate:
        sp_name = confirm_species_name(file_name)
    else:
        sp_name = file_name[:-3].replace(' ', '_').replace('.', '')
    
    if sp_name:
        new_file_name = 'validated_gbk/{0}_validated.gb'.format(sp_name)
        
        if not os.path.isfile(new_file_name):
            with open(some_genbank, 'r') as open_file:
                
                lt_counter = 0
                contig_counter = 0
                genome = SeqIO.parse(open_file, 'gb')
                new_genome = []

                for record in genome:
                    new_genome.append(record)
                new_genome.sort(key=len, reverse=True)
                
                for record in new_genome:
                    contig_counter += 1
                    try:
                        source = record.annotations['source']
                    except KeyError:
                        record.annotations['source'] = 'None'
                    record.annotations['comment'] = """Validated gb file for use with Kvasir HGT finder www.github.com/kescobo/kvasir
                    Original Def: {}
                    Original source: {}
                    """.format(record.name, record.annotations['source'])
                    record.name = 'kvc_{0}'.format(str(contig_counter).zfill(3))
                    record.annotations['source'] = sp_name
                    for feature in record.features:
                        lt_counter += 1
                        feature.qualifiers['kvtag'] = str(lt_counter).zfill(5)
                
                if not os.path.isdir('validated_gbk/'):
                    os.makedirs('validated_gbk/')
                with open(new_file_name, 'w+') as output_handle:
                    SeqIO.write(new_genome, output_handle, 'gb')
                    return os.path.abspath(new_file_name)
        else:
            print 'already fixed, skipping validation!'
            return False

if __name__ == '__main__':
    os.chdir('/Users/KBLaptop/computation/tmp')
    validate_gbk('/Users/KBLaptop/computation/genomes/arthro_IMG.gb', False)