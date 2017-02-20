from settings import *
import sys
sys.path.append('src/')
from DataImport.mongo_import import mongo_import_genbank
from FindHGT.make_blast_db import make_blast_db, db_to_fna
from FindHGT.run_blast import blast_all, parse_blast_results_xml
from DataImport.mongo_import import mongo_import_record
from Analysis import output
from tempfile import NamedTemporaryFile
import os


def import_data():
    for f in os.listdir(INPUT):
        if f.endswith(".gb") or f.endswith(".gbk"):
            print("** Importing {}".format(f))
            mongo_import_genbank(os.path.join(INPUT, f), "genes")  # Perhaps settings.py should include option for collection name?

def blast_db():
    fasta = db_to_fna('genes')  # Collection name option? (see ln15 above)

    # Make separate directory in output for Blast databases. Will probably do this for multiple outputs, might be good
    # to have a function in `Analysis.output`
    db_path = os.path.join(OUTPUT, MONGODB.name, "blast_db")
    if not os.path.isdir(db_path):
        os.makedirs(db_path)

    make_blast_db(fasta.name, "nucl", os.path.join(db_path, "genes"))  # Collection name option? (see ln15 above)


def blast():
    db_path = os.path.join(OUTPUT, MONGODB.name, "blast_db", "genes")

    for species in MONGODB["genes"].distinct("species"):
        print("blasting {}".format(species))
        tmp_file = NamedTemporaryFile(mode="w+")
        for record in MONGODB["genes"].find({"type": "CDS", "species":species}):
            tmp_file.write(">{}\n{}\n".format(
                record["_id"],
                record["dna_seq"]
                )
            )
        tmp_file.seek(0)
        blast_results = blast_all(tmp_file, db_path, perc_identity=0.50)
        for result in parse_blast_results_xml(blast_results):
                mongo_import_record(result, "blast_results")


def analyze(minimum_identity, minimum_length=500, dist_between_hits=5000, ssu_max=0.9):
    groups = output.hgt_groups(minimum_identity, minimum_length, dist_between_hits, ssu_max)
    output.output_groups(
        groups, os.path.join(
            OUTPUT, MONGODB.name, "{}-{}-{}-groups.csv".format(
                int((minimum_identity)*100), minimum_length, dist_between_hits
                )
            )
        )

"""
WIP
"""
def run_circos():
    pass
    # circos.get_karyotypes()
    # circos.get_links(0.99, 500, group_no=2, min_dist=3000, annotation="transposase")
    # circos.get_links(0.99, 500, group_no=1, min_dist=3000, annotation="transposase")
    # circos.get_links(0.90, 400, 0.95)
    # circos.get_gc_conf()
    # circos.get_conf_file("path/to/output/database_name/circos/links/99-400-links.txt",
    #                      gc="path/to/output/database_name/circos/GC/gc_plots.conf")
    # circos.get_conf_file("path/to/output/database_name/circos/links/95-400-links.txt",
    #                      name="ciros_95",
    #                      gc="path/to/output/database_name/circos/GC/gc_plots.conf")
    # circos.get_conf_file("path/to/output/database_name/circos/links/90-400-links.txt",
    #                      name="circos_90",
    #                      gc="path/to/output/database_name/circos/GC/gc_plots.conf")
