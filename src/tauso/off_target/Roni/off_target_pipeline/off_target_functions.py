from src.tauso.hybridization.fast_hybridization import get_trigger_mfe_scores_by_risearch
from Bio import SeqIO
import pandas as pd
from io import StringIO
import os
import pickle

"""
This file contains helper functions both getting the mRNA sequences and for off-target calculations
"""

def dna_to_rna_reverse_complement(seq: str) -> str:
    seq = seq.upper()
    translation_table = str.maketrans("ATGC", "UACG")
    # Translate and reverse
    return seq.translate(translation_table)[::-1]


def normalize_chrom(chrom: str) -> str:
    # Lowercase safe removal of exact "chr" prefix only
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]

    # Mitochondrial chromosome normalization
    if chrom in {"M", "m"}:
        return "MT"

    return chrom

def get_min_max_coords(tup_list):
    flatt = (x for tup in tup_list for x in tup)
    abs_min = abs_max = next(flatt)

    for x in flatt:
        if x < abs_min:
            abs_min = x
        if x > abs_max:
            abs_max = x

    return (abs_min, abs_max)

def parse_risearch_output(output_str: str) -> pd.DataFrame:
    columns = ["trigger", "trigger_start", "trigger_end", "target", "target_start", "target_end", "score", "energy"]
    df = pd.read_csv(StringIO(output_str.strip()), sep="\t", header=None, names=columns)
    return df

def aggregate_off_targets(df: pd.DataFrame) -> pd.DataFrame:
    # Aggregate: sum score (if it's hybridization hits) and take minimum (strongest) energy
    grouped = df.groupby("target").agg({
        "score": "sum",
        "energy": "min"
    }).reset_index()
    return grouped

def parse_gtf(gtf_path, cache_path=None):
    import pickle, os

    if not cache_path:
        cache_path = os.path.splitext(os.path.basename(gtf_path))[0] + "_gtf_annotations.pkl"

    if os.path.exists(cache_path):
        with open(cache_path, 'rb') as f:
            return pickle.load(f)

    annotations = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) != 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            if feature not in ["exon", "CDS", "UTR"]:
                continue

            attr_dict = {
                key: val.strip('"') for key, val in
                [attr.strip().split(' ') for attr in attributes.split(';') if attr.strip()]
            }
            transcript_id = attr_dict.get("transcript_id")
            if not transcript_id:
                continue

            start, end = int(start), int(end)

            if transcript_id not in annotations:
                annotations[transcript_id] = {
                    "chrom": chrom,
                    "strand": strand,
                    "exons": [],
                    "cds": [],
                    "utrs": [],
                    "feature_coords": []  # will collect all coords for accurate span
                }

            # Add to appropriate list
            feature_map = {
                "exon": "exons",
                "CDS": "cds",
                "UTR": "utrs"
            }
            key = feature_map.get(feature)
            if key:
                annotations[transcript_id][key].append((start, end))

            # Collect for overall start/end calculation
            annotations[transcript_id]["feature_coords"].append((start, end))

    # Now compute correct start/end from all features
    for coords in annotations.values():
        coords['start'], coords['end'] = get_min_max_coords(coords['feature_coords'])

    with open(cache_path, 'wb') as f:
        pickle.dump(annotations, f)

    return annotations


def parse_fasta(fasta_path, output_pickle):
    fasta_dict = {}

    for record in SeqIO.parse(fasta_path, "fasta"):
        chrom_id = record.id.split()[0]  # e.g., "1", "KI270728.1", "GL000225.1"
        fasta_dict[chrom_id] = str(record.seq).upper()

    print(f"✅ Saved {len(fasta_dict)} main chromosomes to {output_pickle}")
    with open(output_pickle, "wb") as f:
        pickle.dump(fasta_dict, f)

    return fasta_dict

name2accession = {
    'chr1':  'CM000663.2',
    'chr2':  'CM000664.2',
    'chr3':  'CM000665.2',
    'chr4':  'CM000666.2',
    'chr5':  'CM000667.2',
    'chr6':  'CM000668.2',
    'chr7':  'CM000669.2',
    'chr8':  'CM000670.2',
    'chr9':  'CM000671.2',
    'chr10': 'CM000672.2',
    'chr11': 'CM000673.2',
    'chr12': 'CM000674.2',
    'chr13': 'CM000675.2',
    'chr14': 'CM000676.2',
    'chr15': 'CM000677.2',
    'chr16': 'CM000678.2',
    'chr17': 'CM000679.2',
    'chr18': 'CM000680.2',
    'chr19': 'CM000681.2',
    'chr20': 'CM000682.2',
    'chr21': 'CM000683.2',
    'chr22': 'CM000684.2',
    'chrX': 'CM000685.2',
    'chrY': 'CM000686.2',
    'chrM': 'J01415.2'
}
