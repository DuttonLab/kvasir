from settings import MONGODB as db
from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE


def db_to_fna(collection="genes", seqtype="CDS"):  # See comment on `run.py` ln9
    """
    Takes records of ``"type":seqtype` (like "CDS" or "16s"), writes them to file
    :param collection: string - MongoDB collection name
    :param open_file: file object open for writing
    :return: string - path of file
    """
    tmp_file = NamedTemporaryFile()

    for record in db[collection].find({"type": seqtype}):
        tmp_file.write(">{}\n{}\n".format(
            record["_id"],
            record["dna_seq"]
            )
        )

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
