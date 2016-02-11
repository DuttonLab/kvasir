from settings import *
from DataImport.mongo_import import mongo_import_genbank
from FindHGT.make_blast_db import make_blast_db, db_cds_to_fna
import os


def import_data():
    mongo_import_genbank(INPUT, "genes")  # Perhaps settings.py should include option for collection name?


def blast_db():
    fasta = db_cds_to_fna('genes')  # Collection name option? (see ln9 above)

    # Make separate directory in output for Blast databases. Will probably do this for multiple outputs, might be good
    # to have a function in `Analysis.output`
    db_path = os.path.join(OUTPUT, MONGODB.name, "blast_db")
    if not os.path.isdir(db_path):
        os.makedirs(db_path)

    make_blast_db(fasta.name, "nucl", os.path.join(db_path, "genes"))  # Collection name option? (see ln9 above)


def blast_all():
    pass
