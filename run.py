from settings import *
from FindHGT.create_fasta import db_cds_to_fna
from Analysis import output


def import_data():
    parse_genbank_and_insert(INPUT, "test_collection")


def run_blast():
    db_cds_to_fna('test_collection')


def analyze():
    output.get_islands(0.99, 100)
