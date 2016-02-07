import pymongo
from user_settings import MONGODB
from tempfile import NamedTemporaryFile

filepath = "/Users/KBLaptop/Desktop/test.fna"

client = pymongo.MongoClient()
db = client[MONGODB]


def db_cds_to_fna(collection):
    # fna_file = NamedTemporaryFile()
    # print fna_file.name
    #
    with open(filepath, 'w+') as fna_file:
        for record in db[collection].find({"type": "CDS"}):
            fna_file.write(">{}|{}|{}\n{}\n".format(
                collection,
                record["_id"],
                record["annotation"],
                record["dna_seq"]
                )
            )
        lines = 0

        for line in fna_file:
            if lines < 5:
                print line
                lines += 1
            else:
                break