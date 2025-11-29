import bisect

from pathlib import Path
from typing import Dict

import gffutils
from Bio import SeqIO
from BCBio import GFF
import warnings

from .consts import YEAST_FASTA_PATH, YEAST_FIVE_PRIME_UTR, YEAST_THREE_PRIME_UTR, YEAST_GFF_DB_PATH, YEAST_GFF_PATH

from .file_utils import read_yeast_genome_fasta_dict

from .timer import Timer
from .util import get_longer_string
from .validate import validate_yeast_files
from .genome.LocusInfo import LocusInfo

def cond_print(text, verbose=False):
    if verbose:
        print(text)


def get_locus_to_data_dict_alternative(create_introns=False):
    db_path = Path(YEAST_GFF_DB_PATH)

    if not db_path.exists():
        warnings.warn("Database not found, creating database, this is slow")
        with Timer() as t:
            db = gffutils.create_db(str(YEAST_GFF_PATH), dbfn=str(db_path), force=True, keep_order=True,
                                    merge_strategy='merge', sort_attribute_values=True)
            if create_introns:
                db.update(list(db.create_introns()))  # Uncomment to create introns
        print(f"DB create took: {t.elapsed_time}s")
    else:
        print("Opening DB")
        db = gffutils.FeatureDB(str(db_path))

    fasta_dict = read_yeast_genome_fasta_dict()
    locus_to_data = dict()
    locus_to_strand = dict()

    for feature in db.features_of_type(('CDS', 'intron')):
        chrom = feature.seqid
        if chrom == 'NC_001224.1':  # skip mitochondrial DNA
            continue

        locus_tags = feature.attributes['locus_tag']
        if len(locus_tags) != 1:
            raise ValueError(f"Multiple locuses: {locus_tags}")
        locus_tag = locus_tags[0]

        seq = fasta_dict[chrom].seq[feature.start - 1: feature.end]

        if feature.strand == '-':
            seq = seq.reverse_complement()
        seq = seq.upper()

        if feature.featuretype == 'CDS':
            cds = feature

            if locus_tag not in locus_to_data:
                locus_info = LocusInfo()
                locus_info.exons = [(cds.start, seq)]
                locus_to_data[locus_tag] = locus_info

                locus_to_strand[locus_tag] = cds.strand
            else:
                bisect.insort(locus_to_data[locus_tag].exons, (cds.start, seq))
        elif feature.featuretype == 'intron':
            intron = feature

            if locus_tag not in locus_to_data:
                locus_info = LocusInfo()
                locus_info.introns = [(feature.start, seq)]
                locus_to_data[locus_tag] = locus_info

                locus_to_strand[locus_tag] = intron.strand
            else:
                bisect.insort(locus_to_data[locus_tag].introns, (intron.start, seq))

    for locus_tag in locus_to_data:
        locus_info = locus_to_data[locus_tag]
        if locus_to_strand[locus_tag] == '-':
            locus_info.exons.reverse()
            locus_info.introns.reverse()

        locus_info = locus_to_data[locus_tag]
        locus_info.exons = [element for _, element in locus_info.exons]
        locus_info.introns = [element for _, element in locus_info.introns]
        locus_info.exon_concat = "".join(str(seq) for seq in locus_info.exons)

    return locus_to_data


def get_locus_to_data_dict_yeast():
    locus_to_data = dict()

    with Timer() as t:
        with open(YEAST_FASTA_PATH, 'r') as fasta_handle, open(YEAST_GFF_PATH, 'r') as gff_handle:
            for record in GFF.parse(gff_handle, base_dict=SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))):
                if "NC_001224.1" == record.id:  # Mitochondrial DNA, skip
                    continue

                cond_print(f"> {record.id}")

                cond_print(record.features)
                for feature in record.features:
                    cond_print(f"> {feature.type}")
                    if feature.type == 'gene':
                        cond_print(f" > {feature}")
                        cond_print(f" > {feature.id}")
                        for sub_feature in feature.sub_features:
                            for sub_sub_feature in sub_feature.sub_features:
                                cond_print(f">> {sub_sub_feature}")
                                if sub_sub_feature.type == 'CDS':
                                    cond_print(f"  CDS ID: {sub_sub_feature.id}")
                                    cond_print(f"  Location: {sub_sub_feature.location}")
                                    cond_print(f"  Parent Gene: {sub_sub_feature.qualifiers.get('Parent')}")
                                    if len(sub_sub_feature.qualifiers['locus_tag']) != 1:
                                        raise ValueError(f"More than 1 locus tag on CDS ID: {sub_sub_feature.id}")
                                    locus_tag = sub_sub_feature.qualifiers['locus_tag'][0]
                                    seq = sub_sub_feature.extract(record.seq)

                                    locus_to_data.setdefault(locus_tag, LocusInfo()).exons.append(seq.upper())

    print(f"Took: {t.elapsed_time}s")

    return locus_to_data


def load_three_prime_utr(locus_to_data):
    for record in SeqIO.parse(YEAST_THREE_PRIME_UTR, "fasta"):
        locus = record.name.split('_')[4]
        if locus in locus_to_data:
            locus_info = locus_to_data[locus]

            locus_info.three_prime_utr = get_longer_string(record.seq.upper()[:-1],
                                                           locus_info.three_prime_utr)


def load_five_prime_utr(locus_to_data):
    for record in SeqIO.parse(YEAST_FIVE_PRIME_UTR, "fasta"):
        locus = record.name.split('_')[4]
        if locus in locus_to_data:
            locus_info = locus_to_data[locus]

            locus_info.five_prime_utr = get_longer_string(record.seq.upper()[:-1],
                                                          locus_info.five_prime_utr)


def get_full_locus_to_data() -> Dict[str, LocusInfo]:
    validate_yeast_files()
    locus_to_data = get_locus_to_data_dict_alternative()

    load_five_prime_utr(locus_to_data)
    load_three_prime_utr(locus_to_data)
    for locus_name, locus_info in locus_to_data.items():
        locus_info.full_mrna = f"{locus_info.five_prime_utr}{locus_info.exon_concat}{locus_info.three_prime_utr}"
    return locus_to_data
