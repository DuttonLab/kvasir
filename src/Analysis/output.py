from settings import MONGODB as db
from bson.objectid import ObjectId

def get_hgt(minimum_identity, minimum_length=100):
    """ Generator yielding blast hits >= certain length and >= certain identity

    :param minimum_identity: lowest value of perc_identity to consider (0.0 : 1.0)
    :type minimum_identity: Float
    :param minimum_length:
    :type minimum_length: Int
    :rtype generator: query id and subject id for each hit
    """
    hgt = db['blast_results'].find(
        {'perc_identity': {'$gte': minimum_identity},
         'length'       : {'$gte': minimum_length}}
    )

    for record in hgt:
        yield record['query'], record['subject']

def get_islands(minimum_identity, minimum_length=100, dist_between_hits=3000):
    """ Get blast hits within species that are within x base pairs of eachother
    :param minimum_identity:
    :param minimum_length:
    :param dist_between_hits: number of base-pairs between hits considered significant
    :return:
    """
    hit_list = []

    for id1, id2 in get_hgt(minimum_identity, minimum_length):
        hit_list.extend((ObjectId(id1), ObjectId(id2)))

    ids = list(set(hit_list))  # set removes duplicate, but mongo needs list for query

    for species in db['genes'].distinct('species'):
        for record in db['genes'].find({'species': species, '_id': {'$in': ids}}):
            print "huh?"


def collapse_lists(list_of_iterables):
    """ Reduces list of any iterable that can be converted to a set to non-redundant list of lists

    Combines iterables with identical elements and returns list of lists.
    **Example input**: [[1,2,3],[3,4],[5,6,7],[1,8,9,10],[11],[11,12],[13],[5,12]]
    **Example output: [[1,2,3,4,8,9,10],[5,6,7,11,12],[13]]

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
