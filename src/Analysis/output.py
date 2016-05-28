from settings import MONGODB as db
from bson.objectid import ObjectId
import pandas as pd


def get_hgt(minimum_identity, minimum_length=100, maximum_identity=1.0):
    """ Generator yielding blast hits >= certain length and >= certain identity

    :param minimum_identity: lowest value of perc_identity to consider (0.0 : 1.0)
    :type minimum_identity: Float
    :param minimum_length:
    :type minimum_length: Int
    :rtype generator: query id and subject id for each hit
    """
    hgt = db['blast_results'].find(
        {'perc_identity': {'$gte': minimum_identity, '$lte': maximum_identity},
         'length'       : {'$gte': minimum_length}}
    )

    for record in hgt:
        q, s = ObjectId(record['query']), ObjectId(record['subject'])
        if db['genes'].find_one({'_id': q}) and db['genes'].find_one({'_id': s}):  # all hits should refer to a record, but they don't...
            if not db['genes'].find_one({'_id': q})['species'] == db['genes'].find_one({'_id': s})['species']:
                yield q, s, record

def get_islands(minimum_identity, minimum_length=100, dist_between_hits=3000):
    """ Get blast hits within species that are within x base pairs of each other
    :param minimum_identity:
    :param minimum_length:
    :param dist_between_hits: number of base-pairs between hits considered significant
    :return:
    """
    hit_list = []

    for id1, id2, hit in get_hgt(minimum_identity, minimum_length):
        hit_list.extend([id1, id2])

    ids = list(set(hit_list))  # set removes duplicates, but mongo needs list for query
    # Yields list of islands for each species. There's probably a way more efficient way to do this.
    for species_id_list in get_hits_from_species(ids):
        species_hits_list = [db['genes'].find_one({'_id': x}) for x in species_id_list]
        islands = []
        for entry_1 in species_hits_list:
            entry_recorded = False
            for entry_2 in species_hits_list:
                if entry_1 == entry_2:
                    pass
                else:
                    location_1 = entry_1['location']
                    location_2 = entry_2['location']
                    if location_1['contig'] == location_2['contig']:
                        if abs(location_1['end'] - location_2['start']) <= dist_between_hits:
                            entry_recorded = True
                            islands.append([entry_1['_id'], entry_2['_id']])

            if not entry_recorded:
                islands.append([entry_1['_id']])
        yield collapse(islands)


def get_hits_from_species(hits_list):
    """ iterator returning records from each species in list that matches any in list of `_id`s,

    """
    for species in db['genes'].distinct('species'):
        print("---> Getting blast hits for {}".format(species))
        yield [record['_id'] for record in db['genes'].find({'species': species, 'type': 'CDS', '_id': {'$in': hits_list}})]


def collapse(list_of_iterables):
    """ Reduces list of any iterable that can be converted to a set to non-redundant list of lists

    Combines iterables with identical elements and returns list of lists.
    **Example input**: [[1,2,3],[3,4],[5,6,7],[1,8,9,10],[11],[11,12],[13],[5,12]]
    **Example output**: [[1,2,3,4,8,9,10],[5,6,7,11,12],[13]]

    from stackoverflow user YXD: http://stackoverflow.com/questions/30917226/collapse-list-of-lists-to-eliminate-redundancy

    :param list_of_iterables: list containing any iterables (tuples, lists etc) that can be converted to a set
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


def find_all_hits(some_id, minimum_identity, minimum_length=100):
    """ Generator yielding blast hits

    :param some_id: ObjectID
    :return: generator
    """
    collection = db["blast_results"]

    # since we don't know order of insert, check both
    as_query = {"type": "blast_result", "query": str(some_id)}  # maybe should use ObjectId on import, see issue#10
    as_subject = {"type": "blast_result", "subject": str(some_id)}

    blast_hits = collection.find({"$or": [as_query, as_subject],
                             "perc_identity": {"$gte": minimum_identity},
                             "length":        {"$gte": minimum_length}
                              }
                         )  # will evaluate None if no pair is found

    # we only need to get the hit back
    if blast_hits:
        for blast_hit in blast_hits:
            if blast_hit['query'] == str(some_id):
                yield ObjectId(blast_hit['subject'])
            elif blast_hit['subject'] == str(some_id):
                yield ObjectId(blast_hit['query'])


def hgt_groups(minimum_identity, minimum_length=100, dist_between_hits=3000):
    """
    Returns mutilspecies groups of genes as list of lists.

    """
    groups_list = []
    for islands in get_islands(minimum_identity, minimum_length, dist_between_hits):

        # each sublist represents one island...
        for island in islands:
            hit_set = set() # container for hits
            for gene_id in island:
                gene_hits = find_all_hits(gene_id, minimum_identity, minimum_length)

                # Pulls each hit id, then appends it to group_set
                for hit in gene_hits:
                    hit_set.add(hit)
            # add id for hits to island list...
            island.update(hit_set)
            # And add new island (with multiple species) to groups_list
            groups_list.append(list(island))

    # Since each species' islands are built independently, there's a lot of redundancy
    # So... Collapse lists that contain shared elements and deduplicate
    return map(list, collapse(groups_list))


def output_groups(groups_list, output_file, min_group_size=2):
    """
    Returns .csv file with information for each CDS in an HGT group

    - Optional: set minimum number of CDS to be considered a group
    """
    df_index = ['group', 'locus_tag','contig','start','end','strand','annotation', 'id', 'dna_seq']
    df = pd.DataFrame()
    group_no= 0
    groups_list.sort(key=len, reverse=True)

    for group in groups_list:
        if len(group) >= min_group_size:
            group_no += 1
            for entry in group:
                db_handle = db['genes'].find_one({'_id': ObjectId(entry)})
                annotation = db_handle['annotation'].replace(',', '')  # prevents CSV screw-up
                series = pd.Series(
                    [str(group_no).zfill(3),
                    db_handle['locus_tag'],
                    db_handle['location']['contig'],
                    db_handle['location']['start'],
                    db_handle['location']['end'],
                    db_handle['location']['strand'],
                    annotation,
                    str(db_handle['_id']),
                    db_handle['dna_seq']
                    ],
                    index=df_index,
                    name=db_handle['species']
                )
                df=df.append(series)
    print("Creating file at {}".format(output_file))
    df.to_csv(output_file, columns=df_index)
