from subprocess import Popen, PIPE
import os
import logging
import pandas as pd
from itertools import combinations
from tempfile import NamedTemporaryFile

logger = logging.getLogger(__name__)
logger.propagate = True # passes up to parent logger

def get_ani(species_1, species_2, db):
    """Get ANI for pair of species in database
    Note: this ANI calculation is totally unreliable for distantly related
    species. Use only for pairs of genomes in the same genus.

    Args:
      species_1: string name of a species in MongoDB
      species_2: string name of a species in MongoDB
      db: MongoDB pointer

    Returns:
      Nothing
    """
    modpath = os.path.dirname(os.path.realpath(__file__))
    logger.info("Getting ANI for\n    {}\n    {}".format(species_1, species_2))
    f1 = NamedTemporaryFile(mode="w+")
    f2 = NamedTemporaryFile(mode="w+")
    for contig in db["genes"].find({"type":"contig", "species":species_1}):
        f1.write(">{}\n{}\n".format(contig["contig_id"], contig["dna_seq"]))
    for contig in db["genes"].find({"type":"contig", "species":species_2}):
        f2.write(">{}\n{}\n".format(contig["contig_id"], contig["dna_seq"]))

    ani = Popen([
                "ruby", "{}/ani.rb".format(modpath), "--auto", "--quiet",
                "-1", f1.name,
                "-2", f2.name
                ], stdout=PIPE).communicate()[0]

    logger.debug("ani for {} and {}".format(species_1, species_2))
    logger.debug("ani variable: {} | variable type: {}".format(ani, type(ani)))

    if ani:
        logger.info(
            "Compared {} and {}, ANI = {}".format(species_1, species_2, ani))
        return float(ani) / 100
    else:
        logger.warning(
            "Tried to compare {} and {}, but got no ANI. Returning 0".format(species_1, species_2))
        return 0


def get_distance_matrix(db, dtype="ani"):
    """Get DataFrame containing specified distance matrix

    Args:
      db: Pointer to MongoDB
      dtype (Str): Type of distance to export (Default value = "ani")

    Returns:
      Distance matrix (DataFrame)
    """
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
