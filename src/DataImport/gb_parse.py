from Bio import SeqIO
import os
from re import search


def parse_genbank(genbank_file):
    """
    Retrieves data from a genbank file to prepare for import to MongoDB
    :param genbank_file: a file ending with ".gb" or ".gbk" that contains genomic information
    :rtype list: list of gene records for insertion to MongoDB
    """
    with open(genbank_file, 'r') as in_handle:
        records = SeqIO.parse(in_handle, 'gb')
        record_list = []
        contig_ids = []

        # ToDo: Need a log file to record added locus_tags, species names etc, preferably with a way to reference back to original file
        contig_counter = 0
        locus_tag_counter = 0

        for contig in records:
            contig_counter += 1
            print contig_counter
            try:
                species = contig.annotations['source']
            except KeyError:
                # uses filename (without extension) as species name
                species = os.path.splitext(os.path.basename(genbank_file)[0])

            print species

            if contig.id in contig_ids:
                contig.id = "{}_{}".format(contig.id, contig_counter)
            else:
                contig_ids.append(contig.id)

            # todo: does it make sense to make a class for these record types? Ultimately we're just making a dict but...
            record_list.append({
                'type': 'contig',
                'dna_seq': str(contig.seq),
                'contig_id': contig.id,
                'species': species
            })

            for feature in contig.features:
                locus_tag_counter += 1

                try:
                    aa_seq = feature.qualifiers['translation'][0]
                except KeyError:
                    aa_seq = None

                try:
                    locus_tag = feature.qualifiers['locus_tag'][0]
                except KeyError:
                    locus_tag = 'kv_{}'.format(str(locus_tag_counter).zfill(5))

                try:
                    annotation = feature.qualifiers['product'][0]
                except KeyError:
                    annotation = None

                feature_type = None
                if feature.type == 'CDS':
                    feature_type = 'CDS'
                elif feature.type == 'rRNA':
                    if check_16S(feature):
                        feature_type = '16s'
                    else:
                        feature_type = 'rRNA'

                # grabs DNA sequence from record object
                dna_seq = str(feature.extract(contig).seq)
                record_list.append({
                    'type': feature_type,
                    'dna_seq': dna_seq,
                    'aa_seq': aa_seq,
                    'locus_tag': locus_tag,
                    'annotation': annotation,
                    'location': {
                        'start': int(feature.location.start),
                        'end': int(feature.location.end),
                        'strand': feature.location.strand,
                        'contig': contig.id},
                    'species':species
                })
    return record_list


def parse_genbank_and_insert(genbank_file, collection):
    from DataImport.mongo_import import mongo_import_record

    with open(genbank_file, 'r') as in_handle:
        records = SeqIO.parse(in_handle, 'gb')
        record_list = []
        contig_ids = []

        contig_counter = 0
        locus_tag_counter = 0

        for contig in records:
            contig_counter += 1
            print contig_counter
            try:
                species = contig.annotations['source']
            except KeyError:
                # uses filename (without extension) as species name
                species = os.path.splitext(os.path.basename(genbank_file)[0])

            print species

            if contig.id in contig_ids:
                contig.id = "{}_{}".format(contig.id, contig_counter)
            else:
                contig_ids.append(contig.id)

            mongo_import_record({
                'type': 'contig',
                'dna_seq': str(contig.seq),
                'contig_id': contig.id,
                'species': species
                },
                collection
            )

            for feature in contig.features:
                locus_tag_counter += 1

                try:
                    aa_seq = feature.qualifiers['annotation'][0]
                except KeyError:
                    aa_seq = None

                try:
                    locus_tag = feature.qualifiers['locus_tag'][0]
                except KeyError:
                    locus_tag = 'kv_{}'.format(str(locus_tag_counter).zfill(5))

                try:
                    annotation = feature.qualifiers['product'][0]
                except KeyError:
                    annotation = None

                feature_type = None
                if feature.type == 'CDS':
                    feature_type = 'CDS'
                elif feature.type == 'rRNA':
                    if check_16S(feature):
                        feature_type = '16s'
                    else:
                        feature_type = 'rRNA'

                # grabs DNA sequence from record object
                dna_seq = str(feature.extract(contig).seq)

                mongo_import_record({
                    'type': feature_type,
                    'dna_seq': dna_seq,
                    'aa_seq': aa_seq,
                    'locus_tag': locus_tag,
                    'annotation': annotation,
                    'location': {
                        'start': int(feature.location.start),
                        'end': int(feature.location.end),
                        'strand': feature.location.strand,
                        'contig': contig.id},
                    'species':species
                    },
                    collection
                )

def check_16S(feature):
    """
    Check if rRNA feature is a 16s gene
    :param feature: Biopython genbank feature
    :rtype Boolean: True if annotation indicates 16s and length >1000
    """
    if feature.type == 'rRNA':
        try:
            annotation = str(feature.qualifiers['product'])
        except KeyError:
            return False

        if search(r"16[sS]|ssu|SSU", annotation):
            if len(feature) > 1000:
                return True
            else:
                return False