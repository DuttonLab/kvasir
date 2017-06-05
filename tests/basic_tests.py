import pytest
import pymongo

def test_mongo():
    db = pymongo.MongoClient()["db"]
    db["collection"].insert_one({"key1":"value1", "key2":"value2"})

    record = db["collection"].find_one({"key1":"value1"})
    assert record["key2"] == "value2"
