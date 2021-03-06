import pandas as pd
import logging
from bson.objectid import ObjectId
from kvasir.distance import get_distance_matrix

logger = logging.getLogger(__name__)
logger.propagate = True # passes up to parent logger

def get_hits_from_species(species, minimum_identity, db, minimum_length=100, minimum_species_distance=0, dm=None):
    """Generator yielding blast hits for a speices >= certain length and >= certain identity

    If a gene has multiple hits, there will be duplicates. Use `set()` to remove
    duplicates.

    Args:
      species (Str): name of species
      minimum_identity (Float): lowest value of perc_identity to consider (0.0 : 1.0)
      db: pointer to MongoDB
      minimum_length (Int): minimum length of blast hit in bp (Default value = 100)
      minimum_species_distance (Float, optional): Set to constrain comparisons. Eg. a value of  0.1 will exclude species that have greater than 90\% identity (Default value = 0). Requires distance matrix (`dm` argument)
      dm (DataFrame): distance matrix containing

    Returns:
      Generator - each round yields `_id` value (in ObjectId format) for a blast_results hit

    """
    logger.info("=> Getting blast hits for {}".format(species))
    if minimum_species_distance:
        logger.debug("   > using distance matrix".format(species))

    for record in db['blast_results'].find(
            {"$or":[{"query_species": species},{"subject_species":species}],
             'perc_identity': {'$gte': minimum_identity},
             'length'       : {'$gte': minimum_length},
            }):
        logger.debug("     > Found hit in {}:".format(species))
        if not minimum_species_distance or dm.loc[record["query_species"], record["subject_species"]] > minimum_species_distance:
            if record["query_species"] == species:
                logger.debug("       > {}".format(record["query"]))
                yield ObjectId(record["query"])
            elif record["subject_species"] == species:
                logger.debug("       > {}".format(record["subject"]))
                yield ObjectId(record["subject"])


def get_islands(species, minimum_identity, db, minimum_length=100, dist_between_hits=3000, minimum_species_distance=0, dm=None):
    """Get blast hits within species that are within x base pairs of each other

    Args:
      species (Str): name of species
      minimum_identity (Float): lowest value of perc_identity to consider (0.0 : 1.0)
      db: pointer to MongoDB
      minimum_length (Int): minimum length of blast hit in bp (Default value = 100)
      dist_between_hits (Int): maximum distance (in bp) between hits to be considered part of the same island (Default value = 3000)
      minimum_species_distance (Float, optional): Set to constrain comparisons. Eg. a value of  0.1 will exclude species that have greater than 90\% identity (Default value = 0). Requires species distance collection
      dm (Str, optional): type of distance value to use (Default value = None)

    Returns:
        List of lists - each sublist contains `_id` values (as ObjectId) for a gene in `genes` collection from a single island.

    """
    logger.debug("=> Getting islands from {}".format(species))
    hits = get_hits_from_species(
        species, minimum_identity, db,
        minimum_length, minimum_species_distance, dm
        ) # generator

    hit_list = [db['genes'].find_one({'_id': x}) for x in hits]
    hit_list = sorted(hit_list, key=lambda x: (x["location"]["contig"], x["location"]["start"]))
    islands = []

    last = {"contig":None, "end":0}

    logger.debug("   > got hits, building islands".format(species))
    for hit in hit_list:
        loc = hit["location"]
        if loc["contig"] == last["contig"] and loc["start"] - last["end"] <= dist_between_hits:
            islands[-1].append(hit["_id"])
        else:
            islands.append([hit["_id"]])
        last = loc

    logger.debug("    > got islands, delivering".format(species))
    return islands


def collapse(list_of_iterables):
    """Reduces list of any iterable that can be converted to a set to non-redundant list of lists

    Combines iterables with identical elements and returns list of lists.
    **Example input**: [[1,2,3],[3,4],[5,6,7],[1,8,9,10],[11],[11,12],[13],[5,12]]
    **Example output**: [[1,2,3,4,8,9,10],[5,6,7,11,12],[13]]

    from stackoverflow user YXD: http://stackoverflow.com/questions/30917226/collapse-list-of-lists-to-eliminate-redundancy

    Args:
      list_of_iterables: list containing any iterables (tuples, lists etc) that can be converted to a set

    Returns:
      List of lists without overlap in contents.
    """
    result = []
    for l in list_of_iterables:
        s = set(l)

        matched = [s]
        unmatched = []
        # first divide into matching and non-matching groups
        for g in result:
            if s & g:
                matched.append(g)
            else:
                unmatched.append(g)
        # then merge matching groups
        result = unmatched + [set().union(*matched)]

    return sorted(result, key=len, reverse=True)


def find_all_hits(some_id, db, minimum_identity, minimum_length=100):
    """Generator yielding blast hits

    Args:
      some_id: ObjectId
      db: pointer to MongoDB
      minimum_identity (Float): lowest value of perc_identity to consider (0.0 : 1.0)
      minimum_length (Int): minimum length of blast hit in bp (Default value = 100)
      minimum_species_distance (Float, optional): Set to constrain comparisons. Eg. a value of  0.1 will exclude species that have greater than 90\% identity (Default value = 0). Requires species distance collection
    Returns:
      generator - each round returns an `ObjectId` for a gene in `genes` collection

    """
    # since we don't know order of insert, check both
    as_query = {"type": "blast_result", "query": some_id}
    as_subject = {"type": "blast_result", "subject": some_id}

    blast_hits = db["blast_results"].find({"$or": [as_query, as_subject],
                             "perc_identity": {"$gte": minimum_identity},
                             "length":        {"$gte": minimum_length}
                              }
                         )  # will evaluate None if no pair is found

    # we only need to get the hit back
    if blast_hits:
        for blast_hit in blast_hits:
            if blast_hit['query'] == some_id:
                yield ObjectId(blast_hit['subject'])
            elif blast_hit['subject'] == some_id:
                yield ObjectId(blast_hit['query'])


def hgt_groups(minimum_identity, db, minimum_length=100, dist_between_hits=3000, minimum_species_distance=0, dtype="ani"):
    """Returns mutilspecies groups of genes as list of lists.

    Args:
      minimum_identity (Float): lowest value of perc_identity to consider (0.0 : 1.0)
      db: pointer to MongoDB
      minimum_length (Int): minimum length of blast hit in bp (Default value = 100)
      dist_between_hits (Int): maximum distance (in bp) between hits to be considered part of the same island (Default value = 3000)
      minimum_species_distance (Float, optional): Set to constrain comparisons. Eg. a value of  0.1 will exclude species that have greater than 90\% identity (Default value = 0). Requires species distance collection
      dtype (Str, optional): type of distance value to use (Default value = ani). Can use `None` if minimum_species_distance == 0

    Returns:
      list of lists - each sublist is a group containing gene `_id`s
    """
    if minimum_species_distance:
        logger.debug("    > getting distance matrix")
        dm = get_distance_matrix(db, dtype)
    else:
        dm = None

    groups = []
    for species in db["genes"].distinct("species"):
        islands = get_islands(species, minimum_identity, db, minimum_length, dist_between_hits, minimum_species_distance, dm)
        recorded = []
        logger.debug("     > got islands, updating groups")
        for gindx, g in enumerate(groups):
            groupids=list(map(str, g))
            for iindx, i in enumerate(islands):
                islids = list(map(str, i))
                if db["blast_results"].find_one({"$or":[
                        {"query":{"$in":groupids},"subject":{"$in":islids}},
                        {"subject":{"$in":groupids},"query":{"$in":islids}}]}):
                    logger.debug("     > hit found for group {} and island {}".format(gindx, iindx))
                    groups[gindx].update(i)
                    recorded.append(iindx)
        for iindx in range(len(islands)):
            if not iindx in recorded:
                groups.append(set(islands[iindx]))
        groups = collapse(groups)
    logger.info("=> Found {} groups".format(len(groups)))
    return groups

def output_groups(groups_list, output_file, db, min_group_size=2):
    """Returns .csv file with information for each CDS in an HGT group

    - Optional: set minimum number of CDS to be considered a group

    Args:
      groups_list: param output_file:
      output_file (Str): Path to output_file
      db: pointer to MongoDB
      min_group_size (Int): (Default value = 2)

    Returns:
      None. Makes a file at path `output_file`
    """
    df_columns = ['species','group','locus_tag','contig','start','end','strand','annotation','id', 'dna_seq']
    df = pd.DataFrame(columns=df_columns)
    group_no= 0
    groups_list.sort(key=len, reverse=True)

    for group in groups_list:
        if len(group) >= min_group_size:
            group_no += 1
            for entry in group:
                db_handle = db['genes'].find_one({'_id': ObjectId(entry)})
                if db_handle:
                    annotation = db_handle['annotation'].replace(',', '')  # prevents CSV screw-up
                    newrow = pd.DataFrame([
                        [db_handle["species"],
                        str(group_no).zfill(3),
                        db_handle['locus_tag'],
                        db_handle['location']['contig'],
                        db_handle['location']['start'],
                        db_handle['location']['end'],
                        db_handle['location']['strand'],
                        annotation,
                        str(db_handle['_id']),
                        db_handle['dna_seq']
                        ]], columns = df_columns)
                    df=df.append(newrow)
                else:
                    logger.warning("No database entry found for {}".format(entry))

    df = df.sort_values(["group", "species", "contig", "start"])
    logger.info("Creating file at {}".format(output_file))
    df.to_csv(output_file, columns=df_columns)
