# !/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# CC-BY

import re
import os
import KvDataStructures as kv
import Outputs as o
from itertools import combinations, permutations
from matplotlib import pyplot as plt
from Bio import SeqIO
from Bio.Blast import NCBIXML
from bson.objectid import ObjectId
from subprocess import Popen, PIPE
import numpy as np
import pymongo
import plotly.plotly as py

def ssu_perc_id(species_1, species_2):
    s1_ssu = kv.db['16S'].find_one({'species':species_1})
    s2_ssu = kv.db['16S'].find_one({'species':species_2})

    s1_fasta = make_gene_fasta(s1_ssu)
    s2_fasta = make_gene_fasta(s2_ssu)

    out = Popen(
        ['blastn',
        '-query', s1_fasta,
        '-subject', s2_fasta,
        '-outfmt', '10 pident',
        ], stdout=PIPE
    ).communicate()[0]
    if out:
        result = re.split(r'[,\n]', out.rstrip())[0]
        return result

def make_species_fasta(species):
    if not os.path.isdir('fastas'):
        os.makedirs('fastas')    
    fasta = 'fastas/{}.fna'.format(species)
    if not os.path.isfile(fasta):
        with open(fasta, 'w+') as output_handle:
            for record in kv.get_collection(species).find():
                output_handle.write(
                    ">{}|{}\n{}\n".format(record['species'].replace(' ', '_'), record['_id'], record['dna_seq'])
                )
    return fasta

def make_gene_fasta(gene_record):
    if not os.path.isdir('fastas/genes'):
        os.makedirs('fastas/genes')    
    fasta = 'fastas/genes/{}_{}.fna'.format(gene_record['species'], gene_record['kvtag'])
    if not os.path.isfile(fasta):
        with open(fasta, 'w+') as output_handle:
            output_handle.write(
                ">{}|{}\n{}\n".format(gene_record['species'].replace(' ', '_'), gene_record['_id'], gene_record['dna_seq'])
            )
    return fasta

def make_indexed_fasta(species):
    if not os.path.isdir('fastas/'):
        os.makedirs('fastas/')
    id_list =[]
    indexed_species = kv.index_contigs(kv.get_collection(species))
    indexed_species2 = kv.index_contigs(kv.get_collection(species))
    fasta = 'fastas/{}_indexed.fna'.format(species)
    for record in indexed_species:
        id_list.append(record['_id'])
    if not os.path.isfile(fasta):
        with open(fasta, 'w+') as output_handle:
            for record in indexed_species2:
                output_handle.write(
                    ">{}|{}\n{}\n".format(record['species'].replace(' ', '_'), record['_id'], record['dna_seq'])
                )
                
    return (fasta, id_list)

def fasta_blast(species_1, species_2, word_size=28):
    if not os.path.isdir('pairwise_blast'):
        os.makedirs('pairwise_blast')

    s1 = make_species_fasta(species_1)
    s2 = make_species_fasta(species_2)
    results = 'pairwise_blast/{}_{}-blast_results.xml'.format(species_1, species_2)
    
    if not os.path.isfile('pairwise_blast/{}_blastdb.nhr'.format(species_2)):
        Popen(
            ['makeblastdb',
                '-in', s2,
                '-dbtype', 'nucl',
                '-out', 'pairwise_blast/{}_blastdb'.format(species_2),
                '-title', os.path.basename(species_2),
            ]
        ).wait()

    if word_size > 7:    
        out = Popen(
            ['blastn',
            '-query', s1,
            '-db', 'pairwise_blast/{}_blastdb'.format(species_2),
            '-outfmt', '5',
            '-word_size', str(word_size),
            '-out', results,
            ],
        ).communicate()[0]
    else:
        print "Word size is too small, skipping"
        return False

    blast_hits = count_records(results)
    print "Got {} hits".format(blast_hits)
    if blast_hits < 800:
        if word_size > 7:
            print "Not enough hits, trying with word size = {}".format(word_size-2)
            os.remove(results)
            return fasta_blast(species_1, species_2, word_size-2)
        else:
            print "Still not enough hits, but we're giving up trying to get to 800"
            return results
    elif blast_hits >= 800:
        print "Got enough hits!"
        return results
    else:
        print "something went wrong"
        os.remove(results)
        return False

def output_loc_hist(species_1, species_2, ax):
    s1, id_list = make_indexed_fasta(species_1)
    s2 = make_species_fasta(species_2)
    if not os.path.isfile('pairwise_blast/{}_blastdb.nhr'.format(species_2)):
        Popen(
            ['makeblastdb',
                '-in', s2,
                '-dbtype', 'nucl',
                '-out', 'pairwise_blast/{}_blastdb'.format(species_2),
                '-title', os.path.basename(species_2),
            ]
        ).wait()

    indexed_blast = blast_one(s1, 'pairwise_blast/{}_blastdb'.format(species_2))
    concatenated_subject = kv.concat_contigs(kv.get_collection(species_1))
    
    xys = []
    last_end = 0
    
    for i in range(len(indexed_blast))[0::4]:
        # print indexed_blast[i:i+4]

        subject = concatenated_subject[ObjectId(kv.fasta_id_parse(indexed_blast[i])[1])]
        query = kv.get_mongo_record(*kv.fasta_id_parse(indexed_blast[i+1]))
        
        x1 = subject['location']['start']
        if x1 <= last_end:
            x1 = last_end + 1

        x2 = subject['location']['end']
        last_end = x2

        y = float(indexed_blast[i+2])
        print [(x1-0.1, 0), (x1, y), (x2, y), (x2+0.1, 0)]
        xys.extend([(x1-0.1, 0), (x1, y), (x2, y), (x2+0.1, 0)])

    xys.sort(key=lambda xy:xy[0])
    x, y = zip(*xys)

    if len(x) > 200:
        cmap = plt.get_cmap('jet')
        color = cmap(np.random.rand())

        ax.plot(x, y, marker=None, color=color, linestyle='-', label='{}'.format(species_2))
        print "Plotting {} vs {}!".format(species_1, species_2)
        
        ssu_id = ssu_perc_id(species_1, species_2)
        xmin, xmax = plt.xlim()
        ax.plot([xmin, xmax], [ssu_id, ssu_id], color=color, label='{} 16S'.format(species_2))
       
    else:
        print "Not that many hits between {} and {}, skipping!".format(species_1, species_2)

    # plt.fill_between(x, 0, y, color='#B0384D', alpha=0.5)

def blast_one(query, database, word_size=28):
    out, err = Popen(
            ['blastn',
            '-query', query,
            '-db', database,
            '-num_alignments', '1',
            '-outfmt', '10 qseqid sseqid pident length',
            '-word_size', str(word_size),
            ], stdout=PIPE
        ).communicate()
    if out:
        result = re.split(r'[,\n]', out.rstrip())
        return result
    else:
        return None

def get_clocks(species):
    clock_genes = ["dnaG", "frr", "infC", "nusA", "pgk", "pyrG", "rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplK", "rplL", "rplM", "rplN", "rplP", "rplS", "rplT", "rpmA", "rpoB", "rpsB", "rpsC", "rpsE", "rpsI", "rpsJ", "rpsK", "rpsM", "rpsS", "smpB", "tsf"]
    current_species_collection = kv.get_collection(species)
    clock_dict = {}
    for record in current_species_collection.find():
        if any(gene.lower() in record['annotation'].lower() for gene in clock_genes):
            print record['annotation']

def all_by_all(species_1, species_2):
    # results = fasta_blast(species_1, species_2)
    results = 'pairwise_blast/{}_{}-blast_results.xml'.format(species_1, species_2)
    if results:
        with open(results, 'r') as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            hits_list = []
            for blast_record in blast_records:
                qsp, qid = kv.fasta_id_parse(blast_record.query)
                query_record = kv.get_mongo_record(qsp, qid)
                for alignment in blast_record.alignments:
                    asp, aid = kv.fasta_id_parse(alignment.hit_def)
                    alignment_record = kv.get_mongo_record(asp, aid)
                    for hsp in alignment.hsps:
                        if hsp.align_length > 100:
                            pident = float(hsp.positives)/float(hsp.align_length)
                            length = hsp.align_length
                            hits_list.append((query_record, alignment_record))
                        break
                    break
            return hits_list
    else:
        print "Blast didn't work for some reason"

def global_distance(species_1, species_2):
    hits_list = all_by_all(species_1, species_2)
    hits_distance = []
    ssu_distance = None
    for hit in hits_list:
        distance = o.get_gene_distance(str(hit[0]['dna_seq']), str(hit[1]['dna_seq']))
        hits_distance.append(distance)
        if hit[0]['type'] == '16S':
            ssu_distance = distance

    return (ssu_distance, hits_distance)

def plot_global_distance(list_of_species, sp_index):
    focus = list_of_species[sp_index]
    list_copy = list(list_of_species)
    list_copy.remove(focus)

    ssus = []
    median_ys = []
    avg_ys = []
    plt.close()

    for sp in list_copy:
        x_points = []
        y_points = []
        print sp
        data = global_distance(focus, sp)
        print data[0]
        if data[0]:
            for i in range(len(data[1])):
                x_points.append(data[0])
                y_points.append(data[1][i])
        plt.scatter(x_points, y_points, marker='o', color='#B0384D')
        ssus.append(data[0])
        median_ys.append(np.median(data[1]))
        avg_ys.append(np.average(data[1]))

    plt.title(focus)
    plt.plot([x for x in ssus], [y for y in ssus], marker='s', color='#599FA7')
    plt.plot([x for x in ssus], [y for y in median_ys], '^m')
    plt.xlabel("Species Distance")
    plt.ylabel("Gene Distance")
    plt.xlim(xmax=1.5, xmin=0)
    plt.ylim(ymin=0)

    plt.savefig('{}_distances.pdf'.format(focus))

def pairwise_distance(species_1, species_2):
    print "working on {} vs {}".format(species_1, species_2)
    hits = all_by_all(species_1, species_2)
    xy_points = []
    for hit in hits:
        x = len(hit[0]['dna_seq'])
        y = o.get_gene_distance(str(hit[0]['dna_seq']), str(hit[1]['dna_seq']))
        xy_points.append((x,y))
    


    x, y = zip(*xy_points)
    y_average = np.average(y)
    y_std = np.std(y)
    plt.close()
    plt.title("{} vs {}".format(species_1, species_2))
    plt.scatter(x, y, marker='o', color='#B0384D')
    plt.axhline(y=y_average, color='#599FA7')
    plt.axhline(y=y_average+y_std, color='#599FA7', ls=':')
    plt.axhline(y=y_average-y_std, color='#599FA7', ls=':')

    plt.xlabel("Gene Length (bp)")
    plt.ylabel("Gene Distance")
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    plt.savefig('{}_{}_genedistance.pdf'.format(species_1, species_2))

def global_hist(list_of_species, sp_index):
    focus = list_of_species[sp_index]
    list_copy = list(list_of_species)
    list_copy.remove(focus)

    plt.close()

    for sp in list_copy:
        print sp
        x_points = []
        data = global_distance(focus, sp)
        print data[0]
        if data[0] < 1.0:
            data_list = []
            for i in range(len(data[1])):
                data_list.append(data[1][i])
            x_points.append(data_list)
            
    plt.hist(x_points, histtype='bar', stacked=True)

    plt.title(focus)
    plt.xlabel("Species Distance")
    plt.ylabel("Gene Distance")
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
          ncol=4, fancybox=True)

    plt.savefig('{}_hist.pdf'.format(focus))

def pairwise_hist(species_1, species_2):
    print "working on {} vs {}".format(species_1, species_2)
    hits = all_by_all(species_1, species_2)
    xs = []
    for hit in hits:
        x = o.get_gene_distance(str(hit[0]['dna_seq']), str(hit[1]['dna_seq']))
        xs.append(x)
    
    plt.close()
    plt.title("{} vs {}".format(species_1, species_2))
    plt.hist(xs, color='#B0384D',)
    
    plt.xlabel("Gene Distance")
    plt.ylabel("frequency")
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    plt.savefig('{}_{}_distance_hist.pdf'.format(species_1, species_2))

def count_records(blast_xml):
    with open(blast_xml, 'r') as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        counter = 0
        for i in blast_records:
            try:
                i.alignments[0]
                counter +=1
            except IndexError:
                pass
        return counter

def plot_many(list_of_species_pairs):
    fig = plt.figure()
    ypos = 1
    for pair in list_of_species_pairs:    
        if pair[0] == pairs[0][0]:
            ax = fig.add_axes([1,ypos,1,1])
            output_loc_hist(pair[0], pair[1], ax)

    plot_url = py.plot_mpl(fig)
    print plot_url

    # plt.xlabel("Position")
    # plt.ylabel("percent identity")
    # plt.savefig('/Users/KBLaptop/Desktop/try.pdf')

if __name__ == '__main__':
    kv.mongo_init('more_genomes')
    os.chdir('/Users/KBLaptop/computation/kvasir/data/output/more_genomes/')
    ls = kv.get_species_collections()
    print ls
    ls.remove('Arthrobacter_arilaitensis_Re117')
    pairs = []
    for pair in combinations(ls, 2):
        pairs.append((pair[0], pair[1]))
    plot_many(pairs)


    #     if os.path.isfile('{}_{}.pdf'.format(pair[0], pair[1])):
    #         continue
    #     try:
    #         output_loc_hist(pair[0], pair[1])
    #     except RuntimeError:
    #         print "Couldn't compare {} and {}".format(pair[0], pair[1])

    # test()