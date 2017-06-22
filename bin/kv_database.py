#!/usr/bin/env python

import pymongo
import argparse
import os
import logging
from kvasir import database

parser = argparse.ArgumentParser(description='Import genbank files')

parser.add_argument("mongodb", help="The name of MongoDB database")
parser.add_argument("command", help="Which database command to run",
    choices=["delete", "list_species", "list_contigs"])

parser.add_argument("-s", "--species", help="Species on which to perform operations", nargs="+")

parser.add_argument("-c", "--collection", help="Which collection to perform operations on",
    choices=["genes", "blast_results", "species_distance", "all"], default="genes")

parser.add_argument("-v", "--verbose", help="Display info status messages", action="store_true")
parser.add_argument("-q", "--quiet", help="Suppress most output", action="store_true")
parser.add_argument("--debug", help="Set logging to debug", action="store_true")

parser.add_argument("-l", "--log",
    help="File path for log file")

args = parser.parse_args()
DB = pymongo.MongoClient()[args.mongodb]

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

if args.command == "list_species":
    database.list_species(DB, args.collection, args.species)
elif args.command == "delete":
    database.delete_species(DB, args.collection, args.species)
elif args.command == "list_contigs":
    database.list_contigs(DB, args.species)
