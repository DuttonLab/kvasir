import pymongo

INPUT = "/path/to/input/"  # path to folder with genbank files
OUTPUT = "/Users/KBLaptop/computation/kv_data/output/database_name"  # path to output folder
MONGODB = pymongo.MongoClient()["database_name"]
CIRCOS = "/Users/KBLaptop/computation/circos/"  # path to circos binary
