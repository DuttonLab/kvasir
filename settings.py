# Once set up with user variables, change name to `user_settings.py`
import pymongo

INPUT = "/path/to/input/"
OUTPUT = "/path/to/output/"
MONGODB = pymongo.MongoClient()["database_name"]
