import pymongo
import argparse
import os
import sys

from kvasir.mongo_import import mongo_import_genbank

parser = argparse.ArgumentParser(description='Import genbank files')

parser.add_argument("db", help="The name of MongoDB database")
parser.add_argument("-i", "--input", help="File or directory to import", required=True)

args = parser.parse_args()

INPUT = args.input
db = pymongo.MongoClient()[args.db]

success = 0
if os.path.isdir(INPUT):
    for f in os.listdir(INPUT):
        if f.endswith(".gb") or f.endswith(".gbk"):
            print("** Importing {}".format(f))
            mongo_import_genbank(os.path.join(INPUT, f), "genes")
            success = 1
elif os.path.isfile:
    if f.endswith(".gb") or f.endswith(".gbk"):
        print("** Importing {}".format(f))
        mongo_import_genbank(os.path.join(INPUT, f), "genes")
        success = 1

if not success:
    raise ValueError("Input must be genbank file or folder containing genbank files")
