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

    mongo_import.mongo_import_genbank("tests/testdata/jb4_partial.gb", db, "genes")

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
    db["genes"].insert(
        {
        "type":"16s",
        "dna_seq":"acgtactgatgcagct"},
        {
        "type":"16s",
        "dna_seq":"aaatttggggcccc"
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
