#!/usr/bin/env python
from __future__ import division

class SequenceRecord(object):
    def __init__(self, sequence, gene_name, species_name):
        self.sequence = sequence
        self.gene_name = gene_name
        self.species_name = species_name
    def get_fasta(self):
        safe_species_name = self.species_name.replace(' ','_')
        header = '>' + self.gene_name + '_' + safe_species_name
        return header + '\n' + self.sequence + '\n'

class DNARecord(SequenceRecord):

    def complement(self):
        replacement1 = self.sequence.replace('A', 't')
        replacement2 = replacement1.replace('T', 'a')
        replacement3 = replacement2.replace('C', 'g')
        replacement4 = replacement3.replace('G', 'c')
        return replacement4.upper()
    
    def get_gc_content(self):
        length = len(self.sequence)
        g_count = self.sequence.count('G')
        c_count = self.sequence.count('C')
        gc_content = (g_count + c_count) / length
        return gc_content
    

class ProteinRecord(SequenceRecord):

    def get_hydrophobic(self):
        aa_list=['A','I','L','M','F','W','Y','V']
        protein_length = len(self.sequence)
        total = 0
        for aa in aa_list:
            aa = aa.upper()
            aa_count = self.sequence.count(aa)
            total = total + aa_count
        percentage = total * 100 / protein_length
        return percentage

class RNARecord(SequenceRecord):

    def convert_DNA(self):
        replacement = self.sequence.replace('U', 't')
        return replacement.upper()
