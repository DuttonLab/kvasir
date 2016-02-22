from settings import *
from FindHGT.create_fasta import db_cds_to_fna
from Analysis import output, circos


def import_data():
    parse_genbank_and_insert(INPUT, "test_collection")


def run_blast():
    db_cds_to_fna('test_collection')


def analyze():
    groups = output.hgt_groups(0.99, minimum_length=500,)
    output.output_groups(groups, "/Users/KBLaptop/Desktop/99-500-5000.csv")

def run_circos():
    # circos.get_karyotypes()
    # circos.get_links(0.99, 400)
    circos.get_conf_file("/Users/KBLaptop/computation/kv_data/output/database_name/circos/links/99-400-links.txt")