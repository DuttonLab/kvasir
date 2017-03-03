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
parser.add_argument("-s", "spacer", help="Maximum distance between genes to be considered in the same group",
    default="5000")
parser.add_argument("-l", "length", help="Minimum length (nt) for each protein coding gene",
    default="500")
parser.add_argument("-d", "species-distance", help="Minimum distance between species (0-1)",
    default=0)
parser.add_argument("-t", "--distance-type", help="Type of distance to get from db (usable with species-distance)", default="ani")

args = parser.parse_args()

DB = pymongo.MongoClient()[args.mongodb]

if args.command == groups:
    
