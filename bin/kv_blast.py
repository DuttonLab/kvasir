#!/usr/bin/env python

import pymongo
import argparse
import os
import logging
from kvasir.mongo_import import mongo_import_record
from kvasir.make_blast_db import db_to_fna, make_blast_db
from kvasir.run_blast import blast_all, parse_blast_results_xml
from tempfile import NamedTemporaryFile


parser = argparse.ArgumentParser(description='Kvasir BLAST commands')

parser.add_argument("mongodb", help="The name of MongoDB database")
parser.add_argument("-b", "--blastpath", help="path to BLAST database", default="./")
parser.add_argument("command", help="which blast command to run (makedb, blastall, blastone)",
    choices=["makedb", "blastall", "blastone"])
parser.add_argument("-f", "--force", help="Override errors (eg re-importing blast for same species)", action="store_true")

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


DB = pymongo.MongoClient()[args.mongodb]
BLASTPATH = os.path.abspath(args.blastpath)

if args.command == "makedb":
    logging.info("making BLAST database with protein coding sequences from {} at {}".format(args.mongodb, BLASTPATH))
    fasta = db_to_fna(DB, "genes")

    if not os.path.isdir(BLASTPATH):
        os.makedirs(BLASTPATH)

    make_blast_db(fasta.name, "nucl", os.path.join(BLASTPATH, "{}_CDS".format(args.mongodb)))

elif args.command == "blastall":
    blast_db = os.path.join(BLASTPATH, "{}_CDS".format(args.mongodb))

    for species in DB["genes"].distinct("species"):
        if DB["blast_results"].find_one({"query_species":species}):
            logging.warn("Blast results for {} already present in database".format(species))
            if args.force:
                logging.warning("Adding results for {} anyway".format(species))
            else:
                logging.warning("Skipping {} -- use '-f' option to override".format(species))
        logging.info("blasting {}".format(species))

        tmp_file = NamedTemporaryFile(mode="w+")
        for record in DB["genes"].find({"type": "CDS", "species":species}):
            tmp_file.write(">{}\n{}\n".format(
                record["_id"],
                record["dna_seq"]
                )
            )
        tmp_file.seek(0)
        blast_results = blast_all(tmp_file, blast_db, perc_identity=0.50)
        for result in parse_blast_results_xml(DB, blast_results):
            mongo_import_record(result, DB, "blast_results")


elif args.command == "blastone":
    # TODO: Add command for just searching a single genome
    logging.error("Blastone doesn't work yet. Come back later!")
