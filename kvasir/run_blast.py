from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE
from bson.objectid import ObjectId
from Bio.Blast import NCBIXML
import logging
import os


def blast_all(query_fasta, blast_db, perc_identity=.95):
    """BLAST `query_fasta` against a database
    :param query_fasta: open fasta file object
    :param blast_db: string - path to blast database
    :return: temporary xml file with blast results
    """
    tmp_results = NamedTemporaryFile(mode="w+")
    logging.info("Blasting all by all")
    Popen(
        ['blastn',
         '-query', query_fasta.name,
         '-db', blast_db,
         '-out', tmp_results.name,
         '-perc_identity', str(perc_identity*100),
         '-outfmt', '5'],
         stdout=PIPE  # xml output
    ).wait()  # waits for return code before proceeding

    return tmp_results


def parse_blast_results_xml(db, results_file, seqtype="CDS"):
    """ Parse and insert results of BLAST

    :param results_file: blast results in xml format
    """
    counter = 1
    logging.info("Getting Blast Records")
    try:
        for blast_record in NCBIXML.parse(results_file):

            query_id = blast_record.query

            for alignment in blast_record.alignments:
                hit_id = alignment.hit_def
                same, sp1, sp2 = check_species(db, query_id, hit_id)

                if query_id != hit_id and not same and not check_blast_pair(db, query_id, hit_id):
                    counter += 1
                    if counter % 500 == 0: logging.info("---> {} blast records retrieved".format(counter))

                    hsp = alignment.hsps[0]  # there may be multiple hsps - first one is typically best match
                    perc_identity = float(hsp.identities) / float(hsp.align_length)

                    yield {
                        "type": "blast_result",
                        "seqtype":seqtype,
                        "query": query_id,
                        "query_species": sp1,
                        "subject": hit_id,
                        "subject_species":sp2,
                        "perc_identity": perc_identity,
                        "length": hsp.align_length,
                        "bit_score": hsp.bits,
                        "e-value": hsp.expect
                        }

    except ValueError:
        logging.warning("No hits for {}".format("species"))


def check_blast_pair(db, query, subject):
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


def check_species(db, query, subject):
    """ Check to ensure species are not the same
    - Defaults to checking names (hits from same species are discarded).
    """
    collection = db["genes"]
    query_id, subject_id = ObjectId(query), ObjectId(subject)

    species1 = collection.find_one({"_id":query_id})["species"]
    species2 = collection.find_one({"_id":subject_id})["species"]

    try:
        l1 = " ".join(species1.lower().split()[0:2])
        l2 = " ".join(species2.lower().split()[0:2])
        if l1 == l2:
            return True, species1, species2
        else:
            return False, species1, species2
    except:
        logging.warning("Something went wrong when checking species")
        return False, "", ""
