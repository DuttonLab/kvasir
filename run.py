import user_settings
import os
import sys

def import_data():
    from DataImport.gb_parse import parse_genbank
    from DataImport.mongo_import import mongo_import

    in_dir = user_settings.INPUT
    out_dir = user_settings.OUTPUT

    in_file = os.path.join(in_dir, sys.argv[1])

    mongo_import(parse_genbank(in_file), 'test_collection')


def run_blast():
    