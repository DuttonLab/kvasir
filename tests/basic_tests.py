import pytest
import pymongo

def test_record_import():
    from kvasir import mongo_import

    db = pymongo.MongoClient()["test_db"]
    db["collection"].drop()

    mongo_import.mongo_import_record(
        {"import_key1":"import_value1",
         "import_key2":"import_value2"}, db, "collection")

    record = db["collection"].find_one({"import_key1":"import_value1"})
    assert record["import_key2"] == "import_value2"
