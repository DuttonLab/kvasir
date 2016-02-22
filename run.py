from settings import *
from DataImport.mongo_import import mongo_import_genbank
from Analysis import output


def import_data():

    mongo_import_genbank(INPUT, "collection")


def run_blast():

    db_cds_to_fna('test_collection')


def analyze():
    groups = output.hgt_groups(0.99)
    output.output_groups(groups, OUTPUT)

