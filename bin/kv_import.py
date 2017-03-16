#!/usr/bin/env python

import pymongo
import argparse
import os
import logging
from kvasir.mongo_import import mongo_import_genbank

parser = argparse.ArgumentParser(description='Import genbank files')

parser.add_argument("mongodb", help="The name of MongoDB database")
parser.add_argument("-i", "--input", help="File or directory to import", required=True)

parser.add_argument("-v", "--verbose", help="Display debug status messages", action="store_true")
parser.add_argument("-q", "--quiet", help="Suppress most output", action="store_true")
parser.add_argument("--debug", help="set logging to debug", action="store_true")

parser.add_argument("-l", "--log",
    help="File path for log file")

args = parser.parse_args()

logpath = None
if args.log:
    logpath = os.path.abspath(args.log)
    if os.path.isdir(logpath):
        logpath = os.path.join(logpath, "kvasir.log")

if args.debug:
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s || %(levelname)s: %(message)s", filename=logpath)
elif args.verbose:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s || %(levelname)s: %(message)s", filename=logpath)
elif args.quiet:
    logging.basicConfig(level=logging.ERROR, format="%(asctime)s || %(levelname)s: %(message)s", filename=logpath)
else:
    logging.basicConfig(level=logging.WARNING, format="%(asctime)s || %(levelname)s: %(message)s", filename=logpath)


success = 0

INPUT = os.path.abspath(args.input)
DB = pymongo.MongoClient()[args.mongodb]

if os.path.isdir(INPUT):
    for f in os.listdir(INPUT):
        # TODO: figure out way to throw an error if it looks like the same gene/genome is being imported
        if f.endswith(".gb") or f.endswith(".gbk"):
            logging.info("** Importing {} **".format(f))
            mongo_import_genbank(os.path.join(INPUT, f), DB, "genes")
            success = 1
elif os.path.isfile:
    if f.endswith(".gb") or f.endswith(".gbk"):
        logging.inf("** Importing {} **".format(f))
        mongo_import_genbank(os.path.join(INPUT, f), DB, "genes")
        success = 1

if not success:
    raise ValueError("Input must be genbank file or folder containing genbank files")
