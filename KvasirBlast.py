#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

from pymongo import MongoClient
from subprocess import Popen
import os

#from Bio.Blast import NCBIWWW
#result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
#
#print result_handle.read()

#save_file = open("my_blast.xml", "w")
#save_file.write(result_handle.read())
#save_file.close()
#result_handle.close()

from Bio.Blast.Applications import NcbiblastxCommandline

#blast_handle = NcbiblastxCommandline(
#    query="/Users/KBLaptop/computation/hgt/seqs/genomes/breviFixed.faa",
#    db='mongo_test_again',
#    evalue=0.001,
#    outfmt=5,
#    out="test2.xml"
#    )
#print blast_handle
#stdout, stderr = blast_handle()

result_handle = open("test2.xml")

from Bio.Blast import NCBIXML
blast_records = NCBIXML.parse(result_handle)

#for blast_record in blast_records:
#    print blast_record.query
#    print blast_record.database

for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            #if hsp.expect < 0.04:
                print('****Alignment****')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('e value:', hsp.expect)
                print(hsp.query[0:75] + '...')
                print(hsp.match[0:75] + '...')
                print(hsp.sbjct[0:75] + '...')

#blastx_cline = NcbiblastxCommandline(query="opuntia.fasta", db="nr", evalue=0.001, outfmt=5, out="opuntia.xml")
#blastx_cline
#NcbiblastxCommandline(cmd='blastx', out='opuntia.xml', outfmt=5, query='opuntia.fasta',
#db='nr', evalue=0.001)
#>>> print(blastx_cline)
#blastx -out opuntia.xml -outfmt 5 -query opuntia.fasta -db nr -evalue 0.001
#>>> stdout, stderr = blastx_cline()