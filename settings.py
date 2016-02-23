import pymongo

INPUT = "/path/to/input/"  # path to folder with genbank files
OUTPUT = "path/to/output/"  # path to output folder
MONGODB = pymongo.MongoClient()["database_name"]
