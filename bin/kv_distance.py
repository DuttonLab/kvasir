import pymongo
import argparse
import os
import logging
import pandas as pd
from itertools import combinations, product
from kvasir.mongo_import import mongo_import_distance, mongo_import_distance_matrix
from kvasir.distance import get_ani, get_distance_matrix

logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser(description='Kvasir Analysis commands')

parser.add_argument("mongodb", help="The name of MongoDB database")
parser.add_argument("-c", "--command", help="which analysis command to run (ani, distance_matrix)",
    choices=["ani", "distance_matrix"], required=True)

parser.add_argument("-o", "--output", help="File path for output (usable with distance_matrix)", default="./")
parser.add_argument("-i", "--input", help="File path for input (usable with distance_matrix)")
parser.add_argument("-t", "--distance-type", help="Type of distance to get from db (usable with distance_matrix)", default="ani")

parser.add_argument("-f", "--force", help="Overwrite duplicate records", action="store_true")

args = parser.parse_args()

DB = pymongo.MongoClient()[args.mongodb]

if args.command == "ani":
    species = DB["genes"].distinct("species")
    for sp1, sp2 in combinations(species, 2):
        logging.debug
        if sp1.split().split("_")[0].split("_")[0] == sp2.split()[0].split("_")[0]:
            record_exists = DB["species_distance"].find_one(
                    {"$or":[{"species_1":sp1, "species_2":sp2}, {"species_1":sp2, "species_2":sp1}]})
            if record_exists:
                if args.force:
                    DB["species_distance"].delete_many(
                            {"$or":[{"species_1":sp1, "species_2":sp2}, {"species_1":sp2, "species_2":sp1}]})
                    record_exists = None

            if not record_exists:
                d = 1 - get_ani(sp1, sp2, DB)
                logging.info("    ANI = {} | importing".format(d))
                mongo_import_distance(sp1, sp2, d, DB, dtype="ani")
            else:
                logging.error("{} and {} already have an ANI in the database, stopping. Use -f to overwrite".format(sp1, sp2))
                raise Exception("Careful! You may be trying to duplicate data")
elif args.command == "distance_matrix":
    if args.input:
        logging.info("Importing distance matrix")
        p = os.path.abspath(args.input)
        dm = pd.read_csv(p, header=0, index_col=0)
        mongo_import_distance_matrix(dm, DB, args.distance_type)

    elif args.output:
        logging.info("Extracting distance matrix")
        dm = get_distance_matrix(DB, args.distance_type)
        p = os.path.abspath(args.output)
        if os.path.isdir(p):
            p = os.path.join(p, "distance_matrix.csv")

        dm.to_csv(p)
