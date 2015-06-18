#!/usr/bin/env 
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

'''Must have Mongod running, in terminal: `mongod --dbpath path/to/db`'''

import DataImport
import FixGbk
import MakeBlastDB
import KvasirBlast
import os
import sys

print 'Here we go!'
gbk_folder = sys.argv[1]
mongo_db_name = sys.argv[2]

new_folder = '{0}/kvasir/'.format(os.path.abspath(gbk_folder))
if not os.path.isdir(new_folder):
   os.makedirs(new_folder)

print 'Checking and importing files...'
for the_file in os.listdir(gbk_folder):
    path_to_file = '{0}/{1}'.format(gbk_folder, the_file)
    if the_file.endswith('.gb'):
        print 'Checking {0}'.format(the_file)
        validated_file = FixGbk.add_locus_tag(path_to_file)

        print 'Importing {0}'.format(the_file)
        DataImport.import_file(validated_file, mongo_db_name)

MakeBlastDB.make_blast_db(mongo_db_name)
KvasirBlast.blast(mongo_db_name, mongo_db_name)