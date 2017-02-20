from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE
import logging


def db_to_fna(db, collection, seqtype="CDS"):  # See comment on `run.py` ln9
    """
    Takes records of ``"type":seqtype` (like "CDS" or "16s"), writes them to file
    :param db: string - MongoDB database name
    :param collection: string - MongoDB collection name
    :param open_file: file object open for writing
    :return: tmp file IO object
    """
    tmp_file = NamedTemporaryFile(mode="w+")
    for record in db[collection].find({"type": seqtype}):
        tmp_file.write(">{}\n{}\n".format(
            record["_id"],
            record["dna_seq"]
            )
        )
    tmp_file.seek(0)
    return tmp_file


def make_blast_db(fasta_file, record_type, output_path):
    """
    Create blast db from MongoDB collection
    :param fna_file: string - path to fasta file
    :param record_type: 'nucl' or 'prot'
    :param output_path: string - path and name of output blast database (folder must exist)
    :return:
    """
    Popen(
        ['makeblastdb',
         '-in', fasta_file,
         '-dbtype', record_type,
         '-out', output_path,
         ], stdout=PIPE
    ).wait()
    logging.info("BLAST db created at {}".format(output_path))
