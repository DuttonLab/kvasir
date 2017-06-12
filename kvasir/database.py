import logging


def delete_species(db, collection, species):
    """Delete all records in specified collection for given loci

    Args:
      db: pointer to MongoDB
      collection (Str): MongoDB collection name
      species (Str): name of species to delete

    Returns:
      Nothing
    """


    if collection == "genes":
        logging.info("Deleting all gene records for {}".format(species))
        db[collection].delete_many({"species":{"$in":species}})
    elif collection == "blast_results":
        logging.info("Deleting all blast results for {}".format(species))
        db[collection].delete_many({"$or":[
                        {"query_species":{"$in":species}},
                        {"subject_species":{"$in":species}}]})
    elif collection == "species_distance":
        logging.info("Deleting all distance records for {}".format(species))
        db[collection].delete_many({"$or":[
                        {"species_1":{"$in":species}},
                        {"species_2":{"$in":species}}]})
    elif collection == "all":
        logging.info("Deleting records of all types for {}".format(species))
        db["genes"].delete_many({"species":{"$in":species}})
        db["blast_results"].delete_many({"$or":[
                        {"query_species":{"$in":species}},
                        {"subject_species":{"$in":species}}]})
        db["species_distance"].delete_many({"$or":[
                        {"species_1":{"$in":species}},
                        {"species_2":{"$in":species}}]})
    else:
        logging.error("{} is not a valid collection, no action taken".format(collection))
        raise Exception()


def list_species(db, species=[]):
    """Print name of species in "genes" collection to logger
    Can be useful to get exact strings for "species" arguments in other database
    operations.

    Args:
      db: Pointer to MongoDB
      species: optional list of species, restricts output to species in list (Default value = [])

    Returns:
      Nothing
    """

    species_list = db["genes"].distinct("species")
    if species:
        species_list = set(species).intersection(set(species_list))
        print("The \"{}\" database contains:".format(db.name))
        print("\n".join(sorted(species_list)))
        print("but not:")
        print("\n".join(sorted(set(species).difference(set(species_list)))))
    else:
        print("The \"{}\" database contains:".format(db.name))
        print("\n".join(sorted(species_list)))

    return(sorted(species_list))


def list_contigs(db, species=[]):
    """Print all contigs for each species in "genes" collection to logger

    Args:
      db: Pointer to MongoDB
      species: optional list of species, restricts output to species in list (Default value = [])

    Returns:
      Nothing
    """
    species_list = db["genes"].distinct("species")
    not_found = None
    if species:
        species_list = set(species).intersection(set(species_list))
        not_found = set(species).difference(set(species_list))
    for sp in sorted(species_list):
        print(sp)
        for record in db["genes"].find({"species":sp, "type":"contig"}):
            print("    {}".format(record["contig_id"]))

    if not_found:
        logging.info(
            "Did not find entries for {}:\n     ".format("\n    {}".join(not_found)))


def delete_loci(db, locus_tag_list):
    """Delete all records in genes and blast_results for given loci

    Args:
      db: Pointer to MongoDB
      locus_tag_list (List): locus tags to be removed (as strings)

    Returns:
      Nothing
    """
    for record in db["genes"].find({"locus_tag":{"$in":Set(locus_tag_list)}}):
        db["genes"].delete_one({"_id":record["_id"]})
        db["blast_results"].delete_many({"subject":str(record["_id"])})
        db["blast_results"].delete_many({"query":str(record["_id"])})
