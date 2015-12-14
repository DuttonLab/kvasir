import pymongo
from user_settings import MONGODB

client = pymongo.MongoClient()
db = client[MONGODB]


def mongo_import(list_of_records, collection):
    """
    Takes a list of dicts and inserts them into MongoDB
    :param list_of_records: list of dicts containing gene records
    """
    for record in list_of_records:
        db[collection].insert_one(record)

    print db.collection_names(False)
