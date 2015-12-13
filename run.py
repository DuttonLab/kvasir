import settings
import os
import sys

in_dir = settings.INPUT
out_dir = settings.OUTPUT

in_file = os.path.join(in_dir, sys.argv[1])

print in_file

# parse_genbank(os.path.join(in_dir, 'pacbio2_IMG.gbk'))