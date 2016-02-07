import pymongo
from settings import MONGODB

client = pymongo.MongoClient()
db = client[MONGODB]


def mongo_import_list(list_of_records, collection):
    """
    Takes a list of dicts and inserts them into MongoDB
    :param list_of_records: list of dicts containing gene records
    :param collection
    """
    for record in list_of_records:
        db[collection].insert_one(record)

    print db.collection_names(False)

def mongo_import_record(record, collection):
    """
    Insert a single `record` into `collection`
    :param record: gene record
    :param collection: MongoDB collection name
    """
    db[collection].insert_one(record)