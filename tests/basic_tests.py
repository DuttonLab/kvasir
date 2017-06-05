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
