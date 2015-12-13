import user_settings
import os
import sys
import DataImport

in_dir = user_settings.INPUT
out_dir = user_settings.OUTPUT

in_file = os.path.join(in_dir, sys.argv[1])

print in_file
print sys.path
# DataImport.gb_parse.parse_genbank(os.path.join(in_dir, 'pacbio2_IMG.gbk'))