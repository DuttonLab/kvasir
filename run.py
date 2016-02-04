import user_settings
import os
import sys

def import_data():
    from DataImport.gb_parse import parse_genbank_and_insert
    # from DataImport.mongo_import import mongo_import_list

    # mongo_import_list(parse_genbank(user_settings.INPUT), 'test_collection')

    parse_genbank_and_insert(user_settings.INPUT, "test_collection")

def run_blast():
    from FindHGT.create_fasta import db_cds_to_fna

    db_cds_to_fna('test_collection')