#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

class Species(object):
    def __init__(self, species_name, genes=[], contigs=[]):
        self.species_name = species_name
        self.genes = genes
        self.contigs = contigs

    def add_contig(self, contig_name):
        self.contigs.append(contig_name)

    def add_gene(self, Gene):
        self.genes.append(Gene)

class Gene(object):
    def __init__(self, locus_tag, annotation, aa_seq):
        self.locus_tag = locus_tag
        self.annotation = annotation
        #self.dna_seq = dna_seq
        self.aa_seq = aa_seq
        
    def add_description(self, description):
        self.description = description

    def add_location(self, gene_location):
        self.location = gene_location

    #def get_fna(self):
    #    try:
    #        species_name = self.Species.species_name
    #    except:
    #        print 'No Species Defined!'
    #        pass
    #    safe_species_name = self.species_name.replace(' ','_')
    #    header = '>' + self.locus_tag + '_' + safe_species_name
    #    return header + '\n' + self.dna_seq + '\n'

    def get_faa(self):
        #try:
        #    species_name = self.Species.species_name
        #except:
        #    print 'No Species Defined!'
        #    pass
        #safe_species_name = self.species_name.replace(' ','_')
        header = '>' + self.locus_tag# + '_' + safe_species_name
        return header + '\n' + self.aa_seq + '\n'
