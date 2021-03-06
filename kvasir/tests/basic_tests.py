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


def test_genbank_import():
    from kvasir import mongo_import

    db = pymongo.MongoClient()["test_db"]
    db["genes"].drop()

    mongo_import.mongo_import_genbank("kvasir/tests/testdata/jb4_partial.gb", db, "genes")

    record = db["genes"].find_one({"locus_tag":"Ga0099666_132","type":"CDS"})
    assert record["annotation"] == "hypothetical protein"
    assert record["aa_seq"] == "MLKEQDRDLESRVWDAQRIYDKAERSLSWAQSTQQSRGVFKKLLHGASDQTKVDELTWELEQAQQVYEDISQQHDALAAELEAQKDTVERLRERDGELEDRMGFLDSMLFLTENGTSPFAEDQPDSFSSSADSYLGYDDSYDMSSSPVFGLSDDSSDGMEL"
    assert record["dna_seq"].lower() == "atgctaaaagagcaagaccgagacctagaatctcgcgtgtgggatgcgcaacgcatctacgacaaagcggagcgcagcctttcgtgggcgcagtcgacgcagcaatcgcgcggggtgtttaagaaacttctccacggagcatctgatcagacaaaggttgatgagcttacctgggaactcgagcaagcgcagcaggtgtatgaagatataagccagcagcatgacgcacttgcagctgaactcgaggcacagaaggacaccgtggagagactgcgtgagcgcgatggcgagctggaagacagaatgggcttcttggattccatgctgtttttgacagaaaacgggacgagtccatttgctgaagatcaaccagattccttttcctcatcagctgatagctacctcggatatgacgatagctatgacatgagtagctctcctgtttttgggcttagcgatgattcatctgatggaatggaactctaa"


def test_blast_db_creation():
    from kvasir import make_blast_db
    import os

    db = pymongo.MongoClient()["test_db"]
    db["genes"].drop()

    db["genes"].insert({
        "type":"CDS",
        "dna_seq":"aaatttggggcccc",
        })
    db["genes"].insert({
        "type":"16s",
        "dna_seq":"acgtactgatgcagct"
        })

    f1 = make_blast_db.db_to_fna(db, "genes")
    lns = f1.readlines()
    f1.seek(0)
    assert lns[1] == "aaatttggggcccc\n"

    f2 = make_blast_db.db_to_fna(db, "genes", "16s")
    make_blast_db.make_blast_db(f2.name, "nucl", "./tmp")

    assert os.path.isfile("./tmp.nhr")
    assert os.path.isfile("./tmp.nin")
    assert os.path.isfile("./tmp.nsq")

    os.remove("./tmp.nhr")
    os.remove("./tmp.nin")
    os.remove("./tmp.nsq")


def test_blast_run():
    from tempfile import NamedTemporaryFile
    from kvasir import mongo_import
    from kvasir import make_blast_db
    from kvasir import run_blast
    from kvasir.tests.blast_results import genes
    import os

    db = pymongo.MongoClient()["test_db"]
    db["genes"].drop()

    db["genes"].insert_many(genes)

    f1 = make_blast_db.db_to_fna(db, "genes")
    make_blast_db.make_blast_db(f1.name, "nucl", "./tmp")

    fasta = NamedTemporaryFile(mode="w+")
    for record in db["genes"].find({"type": "CDS", "species":"Corynebacterium sp.  JB4"}):
        fasta.write(">{}\n{}\n".format(
            record["_id"],
            record["dna_seq"]
            )
        )
    fasta.seek(0) == 0


    blast_results = run_blast.blast_all(fasta, "tmp", perc_identity=0.99)
    blast_results.seek(0)
    assert blast_results.readline() == '<?xml version="1.0"?>\n'
    blast_results.seek(0)
    assert 998 <= len(list(blast_results)) <= 999

    blast_results.seek(0)
    r = list(run_blast.parse_blast_results_xml(db, blast_results))
    assert len(r) == 13
    r1 = r[0]
    assert r1["type"] == "blast_result"
    assert "query" in r1 and "subject" in r1 and "perc_identity" in r1
    assert "query_species" in r1 and "subject_species" in r1

    for rec in r:
        assert run_blast.check_species(db, rec["query"], rec["subject"])[0] == False

    os.remove("./tmp.nhr")
    os.remove("./tmp.nin")
    os.remove("./tmp.nsq")


def test_analysis():
    from kvasir.tests.blast_results import blast_results, genes
    from kvasir.output import hgt_groups, output_groups
    from bson.objectid import ObjectId
    import os


    db = pymongo.MongoClient()["test_db"]
    db["genes"].drop()
    db["blast_results"].drop()

    db["genes"].insert_many(genes)
    db["blast_results"].insert_many(blast_results)

    groups = hgt_groups(0.99, db)
    assert len(groups) == 1
    assert len(groups[0]) == 25

    s = sorted(groups[0])
    assert s[0] == ObjectId('5939bbad00947043daa06288')

    output_groups(groups, "tmp.csv", db)

    with open("tmp.csv", "r") as infile:
        l = list(infile)
        assert l[4][0:46] == "0,Arthrobacter sp.  JB182,001,Ga0099663_102740"

    os.remove("./tmp.csv")


def test_database_operations():
    from kvasir.tests.blast_results import blast_results, genes
    from bson.objectid import ObjectId
    from kvasir import database

    db = pymongo.MongoClient()["test_db"]
    db["genes"].drop()
    db["blast_results"].drop()

    db["genes"].insert_many(genes)
    db["blast_results"].insert_many(blast_results)

    assert database.list_species(db)[0] == 'Arthrobacter sp.  JB182'

    database.delete_species(db, "genes", ["Arthrobacter sp.  JB182"])
    assert database.list_species(db)[0] == 'Corynebacterium sp.  JB4'


    assert db["blast_results"].find_one({
        "query":"5939bbad00947043daa062a4"})["subject"] == '5939bbb300947043daa06896'
    database.delete_species(db, "blast_results", ["Arthrobacter sp.  JB182"])
    assert db["blast_results"].find_one({
        "query":"5939bbad00947043daa062a4"}) == None
