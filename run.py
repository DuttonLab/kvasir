from settings import *
from DataImport.mongo_import import mongo_import_genbank
from FindHGT.create_fasta import db_cds_to_fna


def import_data():

    mongo_import_genbank(INPUT, "collection")

def run_blast():

    db_cds_to_fna('test_collection')
