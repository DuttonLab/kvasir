import user_settings
import os
import sys
from subprocess import Popen

from DataImport.gb_parse import parse_genbank

in_dir = user_settings.INPUT
out_dir = user_settings.OUTPUT

in_file = os.path.join(in_dir, sys.argv[1])

mongod = Popen(['mongod', '--dbpath', user_settings.MONGOPATH])

mongo_import(parse_genbank(in_file), 'test_collection')

Popen.kill(mongod)