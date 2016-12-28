import os
from bson.objectid import ObjectId
from subprocess import Popen, PIPE
from settings import *
from tempfile import NamedTemporaryFile
from DataImport.mongo_import import mongo_import_record
from FindHGT.run_blast import blast_all, parse_blast_results_xml
from FindHGT.make_blast_db import db_to_fna, make_blast_db
import pandas as pd

def get_ssu(collection="genes"):
    species = []
    last_species = None
    ssu = None
    for s in MONGODB[collection].find({"type":"16s"}):
        this_species = s["species"]
        seq = s["dna_seq"]
        if not this_species in species:
            species.append(this_species)
            if last_species:
                yield ssu
            else:
                pass

            last_species = this_species
            ssu = (s["_id"], seq)

        elif ssu and len(ssu[1]) < len(seq):
            ssu = (s["_id"], seq)

    yield ssu


def write_ssu(collection="genes"):
    tmp_file = NamedTemporaryFile()
    for ssu in get_ssu(collection):
        if ssu:
            tmp_file.write(">{}\n{}\n".format(ssu[0], ssu[1]))
    tmp_file.seek(0, 2)
    return tmp_file


def ssu_blast_db():
    fasta = write_ssu('genes')

    db_path = os.path.join(OUTPUT, MONGODB.name, "blast_db")
    if not os.path.isdir(db_path):
        os.makedirs(db_path)
    make_blast_db(fasta.name, "nucl", os.path.join(db_path, "ssu"))

def blast_ssu():
    fasta = write_ssu('genes')
    # print(fasta.read())

    db_path = os.path.join(OUTPUT, MONGODB.name, "blast_db", "ssu")

    blast_results = blast_all(fasta, db_path, "50")

    for result in parse_blast_results_xml(blast_results, "16s"):
        query = MONGODB["genes"].find_one({"_id":ObjectId(result["query"])})
        subject = MONGODB["genes"].find_one({"_id":ObjectId(result["subject"])})
        dist = {"query":query["species"],
                "subject":subject["species"],
                "perc_identity":result["perc_identity"],
                "type":"ssu_result"}
        mongo_import_record(dist, "blast_results")

def check_ssu_pair(query, subject, ssu_max=0, distance_db="blast_results"):
    query_species = MONGODB["genes"].find_one({"_id":ObjectId(query)})["species"]
    subject_species = MONGODB["genes"].find_one({"_id":ObjectId(subject)})["species"]

    q16 = MONGODB["genes"].find_one({"species":query_species, "type":"16s"})
    s16 = MONGODB["genes"].find_one({"species":subject_species, "type":"16s"})

    if not q16 and not s16:
        return False
    # since we don't know order of insert, check both
    pair = {"type": "ssu_result", "query": query_species, "subject": subject_species}
    reciprocal = {"type": "ssu_result", "query": subject_species, "subject": query_species}

    blast_pair = MONGODB[distance_db].find_one({"$or": [pair, reciprocal], "perc_identity":{'$gte':ssu_max}})  # will evaluate None if no pair is found

    if blast_pair:
        return True
    else:
        return False

def ssu_import_distance_matrix(filename):
    with open(filename, "r") as f:
        dm = pd.read_csv(filename, index_col=0)
        sp = db["genes"].distinct("species")
        for s1 in sp:
            for s2 in sp:
                dist = dm.loc[s1, s2]
                record = {
                    "type": "blast_result",
                    "seqtype":"ssu",
                    "query": s1,
                    "subject": s2,
                    "perc_identity": dist,
                    "length": None,
                    "bit_score": None,
                    "e-value": None,
                }
                mongo_import_record(record, "blast_results")

def get_distance_matrix(collection="ssu"):
    species = MONGODB[collection].distinct("query")
    l = len(species)
    data = {s:[0 for x in range(l)] for s in species}
    dm = pd.DataFrame(data, index=species)

    for s1 in species:
        for s2 in species:
            dist = MONGODB[collection].find_one({"query":s1, "subject":s2})["perc_identity"]
            dm.loc[s1, s2] = dist # / 100.0

    return dm
