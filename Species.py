#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

import os
from itertools import combinations
from subprocess import Popen, PIPE
import KvDataStructures as kv

def species_compare_16S():
    ssu = kv.get_collection('16S')
    fna_file = 'kvasir/16s.fna'
    results = []
    pairs = {}
    with open(fna_file, 'w+') as file_handle:
        for entry in ssu.find():
            print entry
            file_handle.write(
                '>{0}|{1}\n{2}\n'.format(
                    entry['species'],
                    entry['_id'],
                    entry['dna_seq']
                    )
                )
    
        out = Popen(
            ['blastn',
            '-query', fna_file,
            '-subject', fna_file,
            '-outfmt', '10 qseqid sseqid pident'
            ],
            stdout=PIPE
        ).communicate()[0]

        for result in out.split('\n'):
            results.append(
                tuple(result.split(','))
                )
    os.remove(fna_file)
    
    for pair in results:
        if len(pair) == 3:
            if float(pair[2]) > 95:
                species_1 = parse_seqid(pair[0])[0]
                species_2 = parse_seqid(pair[1])[0]
                if species_1 == species_2:
                    pass
                else:
                    if species_1 in pairs:
                        pairs[species_1].update([species_2])
                    else:
                        pairs[species_1] = set([species_2])
                    if species_2 in pairs:
                        pairs[species_2].update([species_1])
                    else:
                        pairs[species_2] = set([species_1])
    for key in pairs:
        pairs[key] = list(pairs[key])                        
    ssu.insert_one({'same_genera':pairs})


def parse_seqid(seqid):
    return (tuple(seqid.split('|')))

if __name__ == '__main__':
    import sys
    kv.mongo_init(sys.argv[1])
    species_compare_16S()


