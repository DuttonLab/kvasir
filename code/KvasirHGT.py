#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# Unless otherwise indicated, licensed under GNU Public License (GPLv3)

import KvDataStructures as kv
import pandas as pd

def core_hgt_stats():
    """
    Returns stats of HGT (number of events etc)
    """
    collection = kv.get_collection('core')
    df_index = ['Total_CDS', 'HGT_CDS', 'Islands']
    df = pd.DataFrame()
    for species in collection.distinct('species'):
        hits = kv.get_collection('hits').find_one({'species':species})['core_hits']
        series = pd.Series([
            sum([1 for x in collection.find({'species':species})]),
            sum([1 for x in hits if hits[x]]),
            len(get_islands(species))
        ], name=species, index=df_index)
        df = df.append(series)

    df.to_csv('stats.csv',  columns=df_index)

def collapse_lists(list_of_lists):
    """
    Reduces list of lists by combining lists with identical elements
    **Example input**: [[1,2,3],[3,4],[5,6,7],[1,8,9,10],[11],[11,12],[13],[5,12]]
    **Example output: [[1,2,3,4,8,9,10],[5,6,7,11,12],[13]]

    from stackoverflow user YXD: http://stackoverflow.com/questions/30917226/collapse-list-of-lists-to-eliminate-redundancy
    """
    result = []
    for l in list_of_lists:
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
    return result

def get_islands(species_name):
    """
    For each species, combines HGT hits co-occurring within 5kb of eachother
    Returns list of lists of `(species, _id)` tuples
    """
        
    islands = []
    species_hits_list = []
    
    # Add mongo_record for each hit in any gene
    all_hits = kv.get_collection('hits')
    species_hits = all_hits.find_one({'species':species_name})['core_hits']

    
    for query_id in species_hits:
        if species_hits[query_id]:
            species_hits_list.append(
                kv.get_mongo_record(species_name, query_id)
                )

    for entry_1 in species_hits_list:
        entry_recorded = False
        for entry_2 in species_hits_list:
            if entry_1 == entry_2:
                pass
            elif entry_1['location']['contig'] != entry_2['location']['contig']:
                pass
            else:
                location_1 = entry_1['location']
                location_2 = entry_2['location']
                if abs(location_1['end'] - location_2['start']) <= 5000:
                    entry_recorded = True
                    islands.append([
                        (entry_1['species'], str(entry_1['_id'])),
                        (entry_2['species'], str(entry_2['_id']))
                    ])
        if not entry_recorded:
            islands.append([(entry_1['species'], str(entry_1['_id']))])

    return collapse_lists(islands)

def core_hgt_groups():
    """
    Returns mutilspecies groups of genes as list of lists.
    
    - Starts with species islands `get_islands()`
    - For each island, if a hit from that island is in another island,
      group them together
    """
    all_hits = kv.get_collection('hits')
    groups_list = []
    for h in all_hits.find():
        current_species_islands = get_islands(h['species'])
        
        # each sublist represents one island...
        for island in current_species_islands:
            if island: # many lists are empty, skip those
                hit_set = set() # container for hits 
                for gene_id in island:
                    gene_hits = h['core_hits'][gene_id[1]]
                    
                    # Pulls each hit id tuple, then appends it to group_set
                    for hit in gene_hits:
                        hit_set.add((hit[0], hit[1]))
                # add id tuples for hits to island list...
                island.update(hit_set)
                # And add new island (with multiple species) to groups_list
                groups_list.append(list(island))

    # Since each species' islands are built independently, there's a lot of redundancy
    # So... Collapse lists that contain shared elements and deduplicate
    return map(list, collapse_lists(groups_list))

def output_groups(min_group_size=2):
    """
    Returns .csv file with information for each CDS in an HGT group

    - Optional: set minimum number of CDS to be considered a group
    """ 
    output_file = '{}_groups.csv'.format(kv.db.name)
    df_index = ['group','kvtag','contig','start','end','strand','annotation','dna_seq']
    df = pd.DataFrame()
    group_no= 0
    groups_list = core_hgt_groups()
    groups_list.sort(key=len, reverse=True)

    for group in groups_list:
        if len(group) >= min_group_size:
            group.sort(key=lambda entry:entry[0])
            group_no += 1
            for entry in group: # Entry is `(species, id)`
                
                db_handle = kv.get_mongo_record(*entry)
                annotation = db_handle['annotation'].replace(',','') # prevents CSV screw-up
                series = pd.Series(
                    [str(group_no).zfill(3),
                    db_handle['kvtag'],
                    db_handle['location']['contig'],
                    db_handle['location']['start'],
                    db_handle['location']['end'],
                    db_handle['location']['strand'],
                    annotation,
                    db_handle['dna_seq']
                    ],
                    index=df_index,
                    name=db_handle['species']
                )
                df=df.append(series)
    df.to_csv(output_file, columns=df_index)

if __name__ == '__main__':
    import os
    kv.mongo_init('reorg')
    os.chdir('/Users/KBLaptop/computation/kvasir/data/output/reorg/')
    output_groups()