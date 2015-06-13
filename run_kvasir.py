#!/usr/bin/env 
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

import DataImport
import FixGbk
import MakeBlastDB
import KvasirBlast
import os
import sys

print 'Here we go!'
gbk_folder = sys.argv[1]
mongo_db_name = sys.argv[2]

print 'Checking and importing files...'
for the_file in os.listdir(gbk_folder):
    path_to_file = '{0}{1}'.format(gbk_folder, the_file)
    if not the_file.startswith('.'):
        if not os.path.isdir(path_to_file):
            print 'Checking {0}'.format(the_file)
            FixGbk.check_dupe_locus_tags(path_to_file)

            print 'Importing {0}'.format(the_file)
            DataImport.import_file('{0}{1}'.format(gbk_folder, the_file), mongo_db_name)

MakeBlastDB.make_blast_db(mongo_db_name)

