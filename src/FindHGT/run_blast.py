from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE
from make_blast_db import db_cds_to_fna
from DataImport.mongo_import import mongo_import_record
from bson.objectid import ObjectId
from Bio.Blast import NCBIXML
from settings import MONGODB as db


def blast_all(blast_db, collection="gene_info"):
    """
    Take all records of "type":"CDS" from a collection and blast against a database
    :param collection: string - MongoDB collection name
    :param blast_db: string - path to blast database
    :return: temporary xml file with blast results
    """
    tmp_results = NamedTemporaryFile()

    # Because I'm doing this in pieces I put this line in this function, but at least for now, this should be the same
    # fna file that is created when making the blast database. That said, in the future we may want to analyze other
    # genomes against the same blastdb and would therefore want to do this again.
    db_fna = db_cds_to_fna(collection)

    Popen(
        ['blastn',
         '-query', db_fna.name,
         '-db', blast_db,
         '-out', tmp_results.name,
         '-outfmt', '5'  # xml output
         ], stdout=PIPE  # not sure if this is necessary unless we want to print/record the output
    )

    return tmp_results


def parse_blast_results_xml(results_file):
    """

    :param results_file:
    :return:
    """

    results_file.seek(0)
    for blast_record in NCBIXML.parse(results_file):
        query_id = blast_record.query

        for alignment in blast_record.alignments:
            hit_id = alignment.hit_def

            if query_id != hit_id:
                hsp = alignment.hsps[0]  # there may be multiple hsps - first one is typically best match
                perc_identity = hsp.identities / hsp.align_length

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


def check_blast_pair(query, subject):
    """Check database for the presence of a blast result for given pair of record '_id's

    Since blasting seq1 vs seq2 should return the same results as seq2 vs seq1 and we don't want to duplicate data, this
    checks in order to determine whether to import or not

    :param query:
    :param subject:
    :rtype: Bool
    """
    collection = db["blast_results"]
    query_id, subject_id = ObjectId(query_id), ObjectId(subject_id)

    # since we don't know order of insert, check both
    pair = {"type": "blast_result", "query": query_id, "subject": subject_id}
    reciprocal = {"type": "blast_result", "query": subject_id, "subject": query_id}

    if collection.find({"$or": [pair, reciprocal]}):
        return True
    else:
        return False
