from settings import *
from DataImport.mongo_import import mongo_import_genbank
from FindHGT.make_blast_db import make_blast_db, db_cds_to_fna
from FindHGT.run_blast import blast_all, parse_blast_results_xml
from Analysis import output
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


def blast():
    fasta = db_cds_to_fna('genes')  # Collection name option? (see ln9 above)
    db_path = os.path.join(OUTPUT, MONGODB.name, "blast_db", "genes")

    blast_results = blast_all(fasta, db_path)
    parse_blast_results_xml(blast_results)


def analyze():
    groups = output.hgt_groups(0.99, minimum_length=500,)
    output.output_groups(groups, "/Users/KBLaptop/Desktop/99-500-5000.csv")


def run_circos():
    # circos.get_karyotypes()
    # circos.get_links(0.99, 400)
    circos.get_conf_file("/Users/KBLaptop/computation/kv_data/output/database_name/circos/links/99-400-links.txt")
