import pymongo
import argparse
import os
import logging
from kvasir.output import hgt_groups, output_groups

parser = argparse.ArgumentParser(description='Kvasir Analysis commands')

parser.add_argument("mongodb", help="The name of MongoDB database")

parser.add_argument("-c", "--command",
    help="which analysis command to run (groups, species_hits)",
    choices=["groups"], required=True)

parser.add_argument("-v", "--verbose", help="Display debug status messages", action="store_true")
parser.add_argument("-q", "--quiet", help="Suppress most output", action="store_true")

parser.add_argument("--identity",
    help="minimum identity for BLAST hits",
    default="0.99", required=True)

parser.add_argument("-s", "spacer",
    help="Maximum distance between genes to be considered in the same group",
    default="5000")
parser.add_argument("-l", "length",
    help="Minimum length (nt) for each protein coding gene",
    default="500")
parser.add_argument("-d", "species-distance",
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

args = parser.parse_args()

if args.verbose:
    logging.basicConfig(level=logging.DEBUG)
elif args.quiet:
    logging.basicConfig(level=logging.WARNING)
else:
    logging.basicConfig(level=logging.INFO)

DB = pymongo.MongoClient()[args.mongodb]

if args.command == "groups":
    g = hgt_groups(float(args.identity), DB, args.length, args.spacer,
        args.species_distance, args.distance_type)

    if os.path.isdir(args.output):
        out = os.path.join(args.output, "{}-{}-{}_groups.csv".format(
                                            int(float(args.identity)*100),
                                            args.length,
                                            args.spacer
                                            ))
    else:
        out = os.path.abspath(args.output)

    output_groups(g, out, db, args.group_size)
