import pymongo
import argparse
import os
import logging
import pandas as pd
from itertools import combinations
from kvasir.mongo_import import mongo_import_distance
from kvasir.distance import get_ani

logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser(description='Kvasir Analysis commands')

parser.add_argument("mongodb", help="The name of MongoDB database")
parser.add_argument("-c", "--command", help="which analysis command to run (ani, distance_matrix)",
    choices=["ani", "distance_matrix"], required=True)
parser.add_argument("-f", "--force", help="Overwrite duplicate records", action="store_true")

args = parser.parse_args()

DB = pymongo.MongoClient()[args.mongodb]


if args.command == "ani":
    species = DB["genes"].distinct("species")
    for sp1, sp2 in combinations(species, 2):
        logging.debug
        if sp1.split()[0] == sp2.split()[0]:
            record_exists = DB["species_distance"].find_one(
                    {"$or":[{"species_1":sp1, "species_2":sp2}, {"species_1":sp2, "species_2":sp1}]})
            if record_exists:
                if args.force:
                    DB["species_distance"].delete_many(
                            {"$or":[{"species_1":sp1, "species_2":sp2}, {"species_1":sp2, "species_2":sp1}]})
                    record_exists = None

            if not record_exists:
                d = get_ani(sp1, sp2, DB)
                logging.info("    ANI = {} | importing".format(d))
                mongo_import_distance(sp1, sp2, d, DB, dtype="ani")
            else:
                logging.error("{} and {} already have an ANI in the database, stopping. Use -f to overwrite".format(sp1, sp2))
                raise Exception("Careful! You may be trying to duplicate data")
