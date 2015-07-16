#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

import fileinput

def split_faa(input_file, split_number):
    line_count = sum(1 for line in open(input_file))
    lines_per_file = line_count / (split_number - 1)

    i = 0
    out_file = open("split0.faa", 'w+')

    for line in fileinput.FileInput(input_file):
        out_file.write(line)
        i+=1
        if i % lines_per_file == 0:
            out_file.close()
            out_file = open("split{0}.faa".format(i/lines_per_file), 'w+')

    out_file.close()  

# for testing
split_faa('./tmp/arthro_IMG.faa', 100)