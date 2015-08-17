 #!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# Unless otherwise indicated, licensed under GNU Public License (GPLv3)

import os
import sys
import KvasirHGT as kh
from KvDataStructures import mongo_init

def get_outputs():
    if not os.path.isdir('hgt_data/'):
       os.makedirs('hgt_data/')
    os.chdir('hgt_data/')

    kh.output_hits_csv()
    kh.output_groups(min_group_size=4)
    kh.output_compare_matrix()
    kh.output_distance_matrix()

if __name__ == '__main__':
    import sys
    os.chdir('/Users/KBLaptop/computation/kvasir/data/output/dutton_genomes2/')
    mongo_init('dutton_genomes2')
    get_outputs()