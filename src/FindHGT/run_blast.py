from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE
from DataImport.mongo_import import mongo_import_record
from bson.objectid import ObjectId
from Bio.Blast import NCBIXML
from settings import MONGODB as db
import os

def blast_all(query_fasta, blast_db):
    """ Take all records of "type":"CDS" from a collection and blast against a database
    :param query_fasta: open fasta file object
    :param blast_db: string - path to blast database
    :return: temporary xml file with blast results
    """
    tmp_results = NamedTemporaryFile()

    print(Popen(
        ['blastn',
         '-query', query_fasta.name,
         '-db', blast_db,
         '-out', tmp_results.name,
         '-outfmt', '5'],
         stdout=PIPE  # xml output
    ).wait())  # waits for return code before proceeding

    return tmp_results


def parse_blast_results_xml(results_file):
    """ Parse and insert results of BLAST

    :param results_file: blast results in xml format
    """
    counter = 0
    results_file.seek(0)
    for blast_record in NCBIXML.parse(results_file):
        counter += 1
        query_id = blast_record.query

        for alignment in blast_record.alignments:
            hit_id = alignment.hit_def

            if query_id != hit_id:
                hsp = alignment.hsps[0]  # there may be multiple hsps - first one is typically best match
                perc_identity = float(hsp.identities) / float(hsp.align_length)

                if not check_blast_pair(query_id, hit_id):
                    mongo_import_record(
                        {
                            "type": "blast_result",
                            "query": query_id,
                            "subject": hit_id,
                            "perc_identity": perc_identity,
                            "length": hsp.align_length,
                            "bit_score": hsp.bits,
                            "e-value": hsp.expect
                        },
                        "blast_results"
                    )
        if counter % 500 == 0:
            print("---> {} blast records imported".format(counter))


def check_blast_pair(query, subject):
    """ Check database for the presence of a blast result for given pair of record '_id's

    Since blasting seq1 vs seq2 should return the same results as seq2 vs seq1 and we don't want to duplicate data, this
    checks in order to determine whether to import or not. Can also be used to get database entry for any gene pair when
    you know the "_id" value of each.

    :param query:
    :param subject:
    :rtype: MongoDB record (dict) or None
    """
    collection = db["blast_results"]
    query_id, subject_id = ObjectId(query), ObjectId(subject)

    # since we don't know order of insert, check both
    pair = {"type": "blast_result", "query": query_id, "subject": subject_id}
    reciprocal = {"type": "blast_result", "query": subject_id, "subject": query_id}

    blast_pair = collection.find_one({"$or": [pair, reciprocal]})  # will evaluate None if no pair is found

    return blast_pair
