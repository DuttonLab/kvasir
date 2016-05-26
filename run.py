from settings import *
from DataImport.mongo_import import mongo_import_genbank
from FindHGT.make_blast_db import make_blast_db, db_cds_to_fna
from FindHGT.run_blast import blast_all, parse_blast_results_xml
from Analysis import output, circos
import os


def import_data():
    mongo_import_genbank(INPUT, "genes")  # Perhaps settings.py should include option for collection name?

def blast_db():
    fasta = db_cds_to_fna('genes')  # Collection name option? (see ln9 above)

    # Make separate directory in output for Blast databases. Will probably do this for multiple outputs, might be good
    # to have a function in `Analysis.output`
    db_path = os.path.join(OUTPUT, MONGODB.name, "blast_db")
    if not os.path.isdir(db_path):
        os.makedirs(db_path)

    make_blast_db(fasta.name, "nucl", os.path.join(db_path, "genes"))  # Collection name option? (see ln9 above)


def blast():
    fasta = db_cds_to_fna('genes')  # Collection name option? (see ln9 above)
    db_path = os.path.join(OUTPUT, MONGODB.name, "blast_db", "genes")

    blast_results = blast_all(fasta, db_path)
    parse_blast_results_xml(blast_results)


def analyze():
    groups = output.hgt_groups(0.99, minimum_length=500, dist_between_hits=5000)
    output.output_groups(groups, os.path.join(OUTPUT, "99-500-5000-groups-test.csv"))
    # groups = output.hgt_groups(0.99, minimum_length=400, dist_between_hits=5000)
    # output.output_groups(groups, os.path.join(OUTPUT, "99-400-5000-groups.csv"))
    # groups = output.hgt_groups(0.99, minimum_length=400, dist_between_hits=5000)
    # output.output_groups(groups, os.path.join(OUTPUT, "99-00-5000-groups.csv"))
    # groups = output.hgt_groups(0.99, minimum_length=300, dist_between_hits=5000)
    # output.output_groups(groups, os.path.join(OUTPUT, "99-300-5000-groups.csv"))
    # groups = output.hgt_groups(0.99, minimum_length=400, dist_between_hits=5000)
    # output.output_groups(groups, os.path.join(OUTPUT, "99-400-5000-groups.csv"))
    # groups = output.hgt_groups(0.99, minimum_length=400, dist_between_hits=5000)
    # output.output_groups(groups, os.path.join(OUTPUT, "99-400-5000-groups.csv"))


def run_circos():
    # circos.get_karyotypes()
    circos.get_links(0.99, 500, group_no=2, min_dist=3000, annotation="transposase")
    circos.get_links(0.99, 500, group_no=1, min_dist=3000, annotation="transposase")
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
