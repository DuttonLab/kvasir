import subprocess
import os
import logging
import pandas as pd
from itertools import combinations
from tempfile import NamedTemporaryFile


def get_ani(species_1, species_2, db):
    """ Get ANI for pair of species in database
    Note: this ANI calculation is totally unreliable for distantly related
    species. Use only for pairs of genomes in the same genus.

    :param species_1: string name of a species in MongoDB
    :param species_2: string name of a species in MongoDB
    :param db: MongoDB pointer
    """
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
                ]).communicate()[0]) / 100


def get_distance_matrix(db, dtype="ani"):
    species = db["genes"].distinct("species")
    dm = pd.DataFrame([[1 for s in species] for s in species], index=species, columns=species)
    for s in species:
        dm.loc[s, s] = 0

    for s1, s2 in combinations(species, 2):
        dist = db["species_distance"].find_one({"$or":[
                                        {"species_1":s1, "species_2":s2},
                                        {"species_1":s2, "species_2":s1}
                                    ]})
        if dist:
            dist = dist["distance"]
            dm.loc[s1, s2] = dist
            dm.loc[s2, s1] = dist

    return dm
