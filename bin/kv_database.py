#!/usr/bin/env python

import pymongo
import argparse
import os
import logging
from kvasir import database

parser = argparse.ArgumentParser(description='Import genbank files')

parser.add_argument("mongodb", help="The name of MongoDB database")
parser.add_argument("command", help="Which database command to run",
    choices=["dedupe", "delete", "list_species"])

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

if args.debug:
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s || %(levelname)s: %(message)s", filename=logpath)
elif args.verbose:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s || %(levelname)s: %(message)s", filename=logpath)
elif args.quiet:
    logging.basicConfig(level=logging.ERROR, format="%(asctime)s || %(levelname)s: %(message)s", filename=logpath)
else:
    logging.basicConfig(level=logging.WARNING, format="%(asctime)s || %(levelname)s: %(message)s", filename=logpath)

if args.command == "list_species":
    database.list_species(DB, args.collection, args.species)
elif args.command == "delete":
    database.delete_species(DB, args.collection, args.species)
elif args.command == "dedupe":
    database.dedupe(DB, args.collection)
