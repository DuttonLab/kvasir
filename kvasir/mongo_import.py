from kvasir.gb_parse import parse_genbank
from itertools import product
import logging

logger = logging.getLogger(__name__)
logger.propagate = True # passes up to parent logger

def mongo_import_record(record, db, collection):
    """Insert a single record (dict) into collection

    Args:
      db: pointer to MongoDB database
      record: gene record
      collection: Name of collection in MongoDB

    Returns:

    """
    db[collection].insert_one(record)

def mongo_import_genbank(genbank_file, db, collection):
    """Parse genbank file and import into MongoDB

    Args:
      db: pointer to MongoDB database
      genbank_file: genbank file containing genomic information
      collection: Name of collection in MongoDB

    Returns:

    """
    for record in parse_genbank(genbank_file):
        logger.debug(record)
        if record["type"] == "ssu":
            mongo_import_record(record, db, "ssu")
        else:
            mongo_import_record(record, db, collection)

def mongo_import_distance(sp1, sp2, distance, db, dtype="ssu"):
    """Add record to MongoDB for species distance (eg 16S or ANI)

    Value for `distance` should be a float between 0:1, where identical species
    have distance == 0. For example, if ANI is calculated, `distance` should be
    1 - ANI.

    Args:
      sp1: name of one of the species
      sp2: name of other species
      distance: species distance calculated, between 0 and 1
      db: pointer to MongoDB database
      dtype: type of distance calculation. Default = "ssu"

    Returns:

    """
    mongo_import_record({
        "type": dtype,
        "species_1":sp1,
        "species_2":sp2,
        "distance": distance,
        }, db, "species_distance")


def mongo_import_distance_matrix(dm, db, dtype="ssu"):
    """Add records to MongoDB for species distance (eg 16S or ANI) from distance matrix

    Value for `distance` should be a float between 0:1, where identical species
    have distance == 0. For example, if ANI is calculated, `distance` should be
    1 - ANI.

    NOTE: species names in index and columns should match species in database.
    If they don't, they will sill be imported, but won't be matched to existing
    records.

    Args:
      dm: pandas dataframe containing pairwise species distances
      db: pointer to MongoDB database
      dtype: type of distance calculation. Default = "ssu"

    Returns:

    """
    for s1, s2 in product(dm.index, dm.columns):
        if s1 == s2:
            pass
        else:
            d = dm.loc[s1, s2]
            logger.debug(d)
            logger.debug(type(d))
            mongo_import_distance(s1, s2, d, db, dtype)
