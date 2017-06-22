from Bio import SeqIO
import os
from re import search
import logging

logger = logging.getLogger(__name__)
logger.propagate = True # passes up to parent logger

def check_16S(feature):
    """Check if rRNA feature is a 16s gene

    Args:
      feature: Biopython genbank feature

    Returns:
      Bool
    """
    if feature.type == 'rRNA':
        try:
            annotation = str(feature.qualifiers['product'])
        except KeyError:
            return False

        if search(r"16[sS]|ssu|SSU", annotation) and len(feature) > 700:
            return True
        else:
            return False
    else:
        return False


def parse_genbank(genbank_file):
    """Generator yielding contig and gene records

    Retrieves data from a genbank file to prepare for import to MongoDB

    Args:
      genbank_file: path to a file ending with ".gb" or ".gbk" that contains genomic information

    Returns:
      generator - each round yields dictionary with a contig record
    """
    with open(genbank_file, 'r') as in_handle:
        records = SeqIO.parse(in_handle, 'gb')

        # TODO: Need a log file to record added locus_tags, species names etc, preferably with a way to reference back to original file
        for record in add_contig_data(records, genbank_file):
            yield record


def add_contig_data(records, genbank_file):
    """Called by `parse_genbank()` - yields contig record, and then individual gene records
    for each gene within a contig.

    Args:
      records: generator - SeqIO parse object
      genbank_file: open file handle for genbank file

    Returns:
      generator - yields contig and gene records (as Dicts)
    """

    current_species = None
    for contig in records:
        if "source" in contig.annotations:
            species = contig.annotations["source"]
        else:
            # uses filename (without extension) as species name
            f = os.path.splitext(os.path.basename(genbank_file))
            logger.warning('{} in file {} does not have \"source\" attribute, using "{}" as species'.format(
                contig.id,
                genbank_file,
                f[0]
            ))
            species = f[0]

        if not species == current_species:
            logger.info("Importing {}".format(species))
            current_species = species

        # TODO: append list of gene records contained within contig?
        contig_record = {
            'type': 'contig',
            'dna_seq': str(contig.seq),
            'contig_id': contig.id,
            'species': species
        }
        yield contig_record

        for feature in add_features(contig, species):
            yield feature

def add_features(contig, species):
    """Called by `add_contig_data()` - yields individual gene records for each gene within a contig.

    Args:
      contig: SeqIO contig record
      species: name of species

    Returns:
      generator (Dict) - gene records for import to MongoDB
    """

    for feature in contig.features:

        if "translation" in feature.qualifiers:
            aa_seq = feature.qualifiers['translation'][0]
        else:
            aa_seq = None

        if "locus_tag" in feature.qualifiers:
            locus_tag = feature.qualifiers['locus_tag'][0]
        else:
            locus_tag = None
        if "product" in feature.qualifiers:
            annotation = feature.qualifiers['product'][0]
        else:
            annotation = None

        feature_type = None
        if feature.type == 'rRNA':
            if check_16S(feature):
                feature_type = '16s'
            else:
                feature_type = 'rRNA'
        else:
            feature_type = feature.type

        dna_seq = str(feature.extract(contig).seq) # grabs DNA sequence from record object
        feature_record = {
            'type': feature_type,
            'dna_seq': dna_seq,
            'aa_seq': aa_seq,
            'locus_tag': locus_tag,
            'annotation': annotation,
            'species': species,
            'location': {
                'start': int(feature.location.start),
                'end': int(feature.location.end),
                'strand': feature.location.strand,
                'contig': contig.id},
        }

        yield feature_record
