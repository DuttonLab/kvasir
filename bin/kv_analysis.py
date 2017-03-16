#!/usr/bin/env python

import pymongo
import argparse
import os
import logging
from kvasir.output import hgt_groups, output_groups

parser = argparse.ArgumentParser(description='Kvasir Analysis commands')

parser.add_argument("mongodb", help="The name of MongoDB database")

parser.add_argument("command",
    help="which analysis command to run (groups, species_hits)",
    choices=["groups"])

parser.add_argument("-v", "--verbose", help="Display debug status messages", action="store_true")
parser.add_argument("-q", "--quiet", help="Suppress most output", action="store_true")

parser.add_argument("--identity",
    help="minimum identity for BLAST hits",
    default="0.99", required=True)

parser.add_argument("-s", "--spacer",
    help="Maximum distance between genes to be considered in the same group",
    default="5000")
parser.add_argument("--length",
    help="Minimum length (nt) for each protein coding gene",
    default="500")
parser.add_argument("-d", "--species-distance",
    help="Minimum distance between species (0-1)",
    default=0)
parser.add_argument("-g", "--group-size",
    help="Minimum number of genes in a group",
    default=2)
parser.add_argument("-t", "--distance-type",
    help="Type of distance to get from db (usable with species-distance)",
    default="ani")
parser.add_argument("-o", "--output",
    help="File path for output (usable with distance_matrix)", default="./")
parser.add_argument("-l", "--log",
    help="File path for log file")

parser.add_argument("--debug",
    help="set logging to debug", action="store_true")

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



DB = pymongo.MongoClient()[args.mongodb]

if args.command == "groups":
    g = hgt_groups(float(args.identity), DB, int(args.length), int(args.spacer),
        float(args.species_distance), args.distance_type)

    if os.path.isdir(args.output):
        out = os.path.join(args.output, "{}-{}-{}_groups.csv".format(
                                            int(float(args.identity)*100),
                                            args.length,
                                            args.spacer
                                            ))
    else:
        out = os.path.abspath(args.output)

    output_groups(g, out, DB, args.group_size)
