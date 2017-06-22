#!/usr/bin/env python

import pymongo
import argparse
import os
import logging
from kvasir.mongo_import import mongo_import_genbank

parser = argparse.ArgumentParser(description='Import genbank files')

parser.add_argument("mongodb", help="The name of MongoDB database")
parser.add_argument("-i", "--input", help="File or directory to import", required=True)

parser.add_argument("-v", "--verbose", help="Display info status messages", action="store_true")
parser.add_argument("-q", "--quiet", help="Suppress most output", action="store_true")
parser.add_argument("--debug", help="Set logging to debug", action="store_true")

parser.add_argument("-l", "--log",
    help="File path for log file")

args = parser.parse_args()

logpath = None
if args.log:
    logpath = os.path.abspath(args.log)
    if os.path.isdir(logpath):
        logpath = os.path.join(logpath, "kvasir.log")
else:
    logpath = ("kvasir.log")

logger = logging.getLogger("Analysis") # create logger
sh = logging.StreamHandler()
fh = logging.FileHandler(logpath)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d,%H:%M:%S')
sh.setFormatter(formatter)
fh.setFormatter(formatter)

# set level based on args
if args.debug:
    sh.setLevel(logging.DEBUG)
    fh.setLevel(logging.DEBUG)
elif args.verbose:
    sh.setLevel(logging.INFO)
    fh.setLevel(logging.INFO)
elif args.quiet:
    sh.setLevel(logging.ERROR)
    fh.setLevel(logging.ERROR)
else:
    sh.setLevel(logging.WARNING)
    fh.setLevel(logging.WARNING)

logger.addHandler(sh) # add handler to logger
logger.addHandler(fh)

success = 0

INPUT = os.path.abspath(args.input)
DB = pymongo.MongoClient()[args.mongodb]

if os.path.isdir(INPUT):
    for f in os.listdir(INPUT):
        if f.endswith(".gb") or f.endswith(".gbk"):
            logging.info("** Importing {} **".format(f))
            mongo_import_genbank(os.path.join(INPUT, f), DB, "genes")
            success = 1
elif os.path.isfile(INPUT):
    if INPUT.endswith(".gb") or f.endswith(".gbk"):
        logging.info("** Importing {} **".format(INPUT))
        mongo_import_genbank(INPUT, DB, "genes")
        success = 1

if not success:
    raise ValueError("Input must be genbank file or folder containing genbank files")
