from settings import MONGODB as db
from gb_parse import parse_genbank


def mongo_import_record(record, collection):
    """ Insert a single record (dict) into collection

    :param record: gene record
    :param collection: Name of collection in MongoDB
    """
    db[collection].insert_one(record)

def mongo_import_genbank(genbank_file, collection):
    """ Parse genbank file and import into MongoDB

    :param genbank_file: genbank file containing genomic information
    :param collection:
    """
    for record in parse_genbank(genbank_file):
        mongo_import_record(record, collection)