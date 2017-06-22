from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE
import logging

logger = logging.getLogger(__name__)
logger.propagate = True # passes up to parent logger

def db_to_fna(db, collection, seqtype="CDS"):
    """Takes records of ``"type":seqtype` (like "CDS" or "16s"), writes them to file

    Args:
      db: string - MongoDB database name
      collection: string - MongoDB collection name
      open_file: file object open for writing
      seqtype:  (Default value = "CDS")

    Returns:
      tmp file IO object

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
    """Create blast db from MongoDB collection

    Args:
      fna_file: string - path to fasta file
      record_type: nucl' or 'prot'
      output_path: string - path and name of output blast database (folder must exist)
      fasta_file:

    Returns:

    """
    Popen(
        ['makeblastdb',
         '-in', fasta_file,
         '-dbtype', record_type,
         '-out', output_path,
         ], stdout=PIPE
    ).wait()
    logger.info("BLAST db created at {}".format(output_path))
