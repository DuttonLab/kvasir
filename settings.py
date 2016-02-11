# Once set up with user variables, change name to `user_settings.py`
import pymongo

INPUT = "/Users/KBLaptop/computation/kv_data/input/img_core/pacbio2_IMG.gbk"
OUTPUT = "/Users/KBLaptop/computation/kv_data/output"
MONGODB = pymongo.MongoClient()["database_name"]
