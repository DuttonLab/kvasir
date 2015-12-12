#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# Unless otherwise indicated, licensed under GNU Public License (GPLv3)

'''
- Create directories `path/input/` and `path/output/`
- Put genbank files into `path/input/core` or `path/input/other`
- Must have Mongod running, in terminal: `mongod --dbpath path/to/db`
- run `python run_kvasir.py [run_name]` from `path/`
'''

import os
import sys
import datetime
import DataImport, FixGbk, KvasirBlast
from KvasirHGT import output_groups
from KvDataStructures import mongo_init, reset_database

print 'Here we go!'

core_gbk_folder = os.path.abspath('data/input/core/')
other_gbk_folder =os.path.abspath('data/input/other/')
exp_name = datetime.datetime.now().strftime('%Y%m%d%t%H%M')

mongo_init(exp_name)
reset_database(exp_name)    

new_folder = 'data/output/{0}/'.format(exp_name)
if not os.path.isdir(new_folder):
   os.makedirs(new_folder)
os.chdir(new_folder)

print 'Checking and importing core genomes...'

for the_file in os.listdir(core_gbk_folder):
    path_to_file = '{0}/{1}'.format(core_gbk_folder, the_file)
    DataImport.import_file(path_to_file, 'core')

print 'Checking and importing other genomes...'

for the_file in os.listdir(other_gbk_folder):
    path_to_file = '{0}/{1}'.format(other_gbk_folder, the_file)
    DataImport.import_file(path_to_file, 'other')
    

KvasirBlast.make_blast_db('core')
KvasirBlast.make_blast_db('other')
KvasirBlast.core_hgt_blast(perc_identity='90')
# KvasirBlast.core_hgt_blast(perc_identity='95')
# KvasirBlast.core_hgt_blast(perc_identity='99')
KvasirBlast.blast_to_db(perc_identity='90')
# KvasirBlast.blast_to_db(perc_identity='95')
# KvasirBlast.blast_to_db(perc_identity='99')
output_groups()

