#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

from Bio.Blast import NCBIXML
import re

with open('./tmp/arthro_all.xml', 'r') as result_handle:
    blast_records = NCBIXML.parse(result_handle)
    grep = re.compile(r'gnl\|all_genomes\|(\w+_IMG|w+Fixed)_[a-zA-Z]+_(\w+)')

    for blast_record in blast_records:
        if len(blast_record.alignments) > 1:
            query_match = grep.match(str(blast_record.query))
            query_species = query_match.group(1)
            query_locus = query_match.group(2)
            print 'hits:'
            
            for hit in blast_record.alignments:
                hit_match = grep.match(str(hit))
                if query_match.group(2) == hit_match.group(2):
                    print 'none'
                    pass
                else:
                    print hit_match.group(1)
                    print hit_match.group(2)