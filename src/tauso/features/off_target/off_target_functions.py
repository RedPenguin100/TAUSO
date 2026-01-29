from Bio import SeqIO
import pandas as pd
from io import StringIO
import pickle

"""
This file contains helper functions both getting the mRNA sequences and for off-target calculations
"""

def dna_to_rna_reverse_complement(seq: str) -> str:
    seq = seq.upper()
    translation_table = str.maketrans("ATGC", "UACG")
    # Translate and reverse
    return seq.translate(translation_table)[::-1]


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
            gene_name = attr_dict.get("gene_name")

            if not transcript_id:
                continue

            start, end = int(start), int(end)

            if transcript_id not in annotations:
                annotations[transcript_id] = {
                    "chrom": chrom,
                    "strand": strand,
                    "gene_name": gene_name,
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
        all_features = coords['feature_coords']
        # Start is the smallest of all 'start' values
        coords['start'] = min(f[0] for f in all_features)

        # End is the largest of all 'end' values
        coords['end'] = max(f[1] for f in all_features)

    if cache_path:
        with open(cache_path, 'wb') as f:
            pickle.dump(annotations, f)

    print("Read GTF file")
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
