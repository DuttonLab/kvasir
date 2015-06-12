#!/usr/bin/env 
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

import DataImport
import FixGbk
import MongoFasta
import MakeBlastDB
import KvasirBlast
import os

gbk_folder = sys.argv[1]
mongo_db_name = sys.argv[2]

try:
    DataImport.import_folder(gbk_folder, mongo_db_name)
    pass
except KeyError, e:
    print "Genbank files don't have locus_tag, fixing"
    raise e