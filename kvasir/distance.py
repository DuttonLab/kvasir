import subprocess
import os
from tempfile import NamedTemporaryFile
import logging


def get_ani(species_1, species_2, db):
    modpath = os.path.dirname(os.path.realpath(__file__))
    logging.info("Getting ANI for\n    {}\n    {}".format(species_1, species_2))
    f1 = NamedTemporaryFile(mode="w+")
    f2 = NamedTemporaryFile(mode="w+")
    for contig in db["genes"].find({"type":"contig", "species":species_1}):
        f1.write(">{}\n{}\n".format(contig["contig_id"], contig["dna_seq"]))
    for contig in db["genes"].find({"type":"contig", "species":species_2}):
        f2.write(">{}\n{}\n".format(contig["contig_id"], contig["dna_seq"]))

    return float(subprocess.Popen([
                "ruby", "{}/ani.rb".format(modpath), "--auto", "--quiet",
                "-1", f1.name,
                "-2", f2.name,
                ]).communicate()[0])


def get_distance_matrix(db, dtype):
