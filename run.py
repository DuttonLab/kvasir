from settings import *
from DataImport.gb_parse import parse_genbank_and_insert
from FindHGT.create_fasta import db_cds_to_fna


def import_data():

    parse_genbank_and_insert(INPUT, "test_collection")

def run_blast():

    db_cds_to_fna('test_collection')