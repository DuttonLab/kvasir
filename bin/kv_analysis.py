import pymongo
import argparse
import os
import logging
from kvasir.output import get_matches

logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser(description='Kvasir Analysis commands')

parser.add_argument("mongodb", help="The name of MongoDB database")
parser.add_argument("-c", "--command", help="which analysis command to run (groups, species_hits)",
    choices=["groups", "species_hits"], required=True)

args = parser.parse_args()

DB = pymongo.MongoClient()[args.mongodb]

if args.command == groups:
