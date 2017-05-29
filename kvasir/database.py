import logging


def delete_species(db, collection, species):
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


def list_species(db, collection, species=[]):
    species_list = db[collection].distinct("species")
    if species:
        species_list = set(species).intersection(set(species_list))
        print("The \"{}\" database contains:".format(db.name))
        print("\n".join(sorted(species_list)))
        print("but not:")
        print("\n".join(sorted(set(species).difference(set(species_list)))))
    else:
        print("The \"{}\" database contains:".format(db.name))
        print("\n".join(sorted(species_list)))


def dedupe(db, collection, species=[]):
    if species:
        if collection == "genes":
            records = db[collection].find({"species":{"$in":species})
        elif collection = "blast_results":
            records = db[collection].find({"$or":[
                {"query_species":{"$in":species}},
                {"subject_species":{"$in":species}}
                }])
        else:
            logging.error("Can only dedupe 'genes' or 'blast_results' collection")
    else:
        records = db[collection].find()

    for record in records:
        if collection == "genes":
            if "locus_tag" in record:
                search = {
                    "species": record["species"]
                    "locus_tag": record["locus_tag"]
                    "dna_seq": record["dna_seq"]
                    }
            elif record["type"] == "contig":
                search = {
                    "species": record["species"]
                    "contig_id": record["locus_tag"]
                    "dna_seq": record["dna_seq"]
                    }
        elif collection = "blast_result":
            search = {
                "$or":[{
                    "query_species":record["query_species"],
                    "query":record["query"],
                    "subject_species":record["subject_species"],
                    "subject":record["subject"]},
                    {"query_species":record["subject_species"],
                    "query":record["subject"],
                    "subject_species":record["query_species"],
                    "subject":record["query"]}],
                "perc_identity": record["perc_identity"],
                "length": record["length"]
            }

        search["_id"] = {"$ne": record["_id"]}

        db[collection].delete_many(search)
