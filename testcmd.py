import os
import sys

print 'Here we go!'
gbk_folder = os.path.abspath(sys.argv[1])

print 'Checking files...'
for the_file in os.listdir(gbk_folder):
    if not the_file.startswith('.'):
        if not os.path.isdir(the_file):
            print the_file