import json
from pathlib import Path

import pandas as pd

from ...consts import DATA_PATH
from ..consts_dataframe import SEQUENCE, CANONICAL_GENE
from ...features.cai import calc_CAI
from ...util import _to_str_seq, get_antisense, _norm_rna_to_dna


PARDIR = Path(__file__).parent

def load_mrna_by_gene_from_files(files: list[str | Path], seq_column: str = "Original Transcript Sequence") -> dict[
    str, str]:
    """
    Load {Gene -> <seq_column>} from a manual list of CSV paths.
    - Expects columns: 'Gene' and <seq_column>
    - Returns DNA alphabet (A/C/G/T) after U->T via _norm_rna_to_dna
    - If multiple rows per gene: keeps the *longest* sequence
    """
    files = [Path(f) for f in files]
    rows = []
    for f in files:
        df = pd.read_csv(f, usecols=['Gene', seq_column])
        df[seq_column] = df[seq_column].map(_norm_rna_to_dna)
        # Keep only clean sequences
        df = df[df[seq_column].str.fullmatch(r'[ACGT]+', na=False)]
        rows.append(df)

    if not rows:
        return {}

    big = pd.concat(rows, ignore_index=True)
    big['len'] = big[seq_column].str.len()
    # Pick longest per gene
    chosen = big.sort_values(['Gene', 'len'], ascending=[True, False]).drop_duplicates('Gene')
    return dict(zip(chosen['Gene'], chosen[seq_column]))


# ---- Choose which mRNA to use for mRNA-based features (tAI/windows on mRNA, etc.) ----
def choose_preferred_mrna(gene_name: str, mrna_built_from_exons: str, gene_to_mrna_real: dict[str, str]) -> str:
    """
    Prefer the real (external) mRNA when available; otherwise fall back to exon-joined.
    Does NOT touch your genome->(mRNA/CDS) mappings or pre-mRNA flanks.
    """
    ext = gene_to_mrna_real.get(gene_name)
    return ext if ext else mrna_built_from_exons


def _build_spliced_mrna_from_exons(pre_mrna: str, exon_indices):
    """
    Build exon-joined mRNA by concatenating exon slices out of pre_mrna.
    Keeps your original assumptions: pre_mrna corresponds to genomic strand and
    starts at exon_indices[0][0]; exon intervals are used directly.
    """
    if not exon_indices:
        return ""
    pre_genome_start = exon_indices[0][0]
    parts = []
    for exon_start, exon_end in exon_indices:
        pm_start = exon_start - pre_genome_start
        pm_end = exon_end - pre_genome_start
        parts.append(pre_mrna[pm_start:pm_end])
    return "".join(parts)


def populate_cai_for_aso_dataframe(aso_df, locus_info, cell_line='A431'):
    """
    At the moment calculation is on pre-mRNA, which is not precisely the Codon Adaptation Index method,
    but we saw this as a useful proxy method for expression and distinction between cell lines.
    """
    SUPPORTED_CELL_LINES = ['A431']
    if cell_line == 'A431':
        MRNA_FILENAME = DATA_PATH / 'human' / 'transcripts' / 'ACH-001328_transcriptome.csv'
    else:
        raise ValueError(f'Supporting only {SUPPORTED_CELL_LINES} at the moment.')

    with open(PARDIR / 'cai_cache' / 'weights_cache.json') as f:
        weights_flat_dict = json.load(f)

    weights_flat = weights_flat_dict[cell_line]

    # Column names
    SENSE_LENGTH = 'sense_length'  # Length of the ASO (nt)
    CDS_SEQUENCE = 'cds_sequence'  # CDS string (joined exons within CDS range)
    IN_CODING_REGION = 'in_coding_region'  # site is within CDS on a real exon

    # Flank sizes
    FLANK_SIZES_PREMRNA = [20, 30, 40, 50, 60, 70]
    FLANK_SIZES_CDS = [20, 30, 40, 50, 60, 70]

    aso_df[CDS_SEQUENCE] = ""
    aso_df[IN_CODING_REGION] = False

    for fs in FLANK_SIZES_PREMRNA:
        aso_df[f"flank_sequence_{fs}"] = ""
    for fs in FLANK_SIZES_CDS:
        aso_df[f"local_coding_region_around_ASO_{fs}"] = ""

    # Cache CDS per gene
    gene_to_cds_info = {}

    # ---- main loop ----
    for index, row in aso_df.iterrows():
        gene_name = row[CANONICAL_GENE]

        # Keep using your current pre-mRNA for flanks/exon-intron logic (coerced to clean string)
        pre_mrna = _to_str_seq(locus_info.full_mrna)
        antisense = _to_str_seq(row[SEQUENCE])
        sense = _to_str_seq(get_antisense(antisense))

        # Locate site on pre-mRNA
        idx = pre_mrna.find(sense)
        aso_df.at[index, SENSE_LENGTH] = len(antisense)

        if idx != -1:
            # Genomic correction (kept as-is)
            genome_corrected_index = idx + locus_info.exon_indices[0][0]

            # pre-mRNA flanks (now using .at and guaranteed string slices)
            for fs in FLANK_SIZES_PREMRNA:
                flank_start = max(0, idx - fs)
                flank_end = min(len(pre_mrna), idx + len(sense) + fs)
                flank_seq = pre_mrna[flank_start:flank_end]
                aso_df.at[index, f"flank_sequence_{fs}"] = flank_seq

            # Build CDS + genome->mRNA map (kept identical to your approach)
            if gene_name not in gene_to_cds_info:
                cds_seq = []  # build as list for speed, join at end
                genome_to_mrna_map = {}
                mrna_idx = 0
                for exon_start, exon_end in locus_info.exon_indices:
                    for gpos in range(exon_start, exon_end):
                        if mrna_idx >= len(pre_mrna):
                            break
                        if locus_info.cds_start <= gpos <= locus_info.cds_end:
                            cds_seq.append(pre_mrna[mrna_idx])
                            genome_to_mrna_map[gpos] = len(cds_seq) - 1
                        mrna_idx += 1
                cds_seq = ''.join(cds_seq)
                gene_to_cds_info[gene_name] = (cds_seq, genome_to_mrna_map)
            else:
                cds_seq, genome_to_mrna_map = gene_to_cds_info[gene_name]

            # Save CDS
            aso_df.at[index, CDS_SEQUENCE] = _to_str_seq(cds_seq)

            # If within CDS, extract local CDS context (unchanged logic; .at + str)
            if (
                    locus_info.cds_start <= genome_corrected_index <= locus_info.cds_end
                    and genome_corrected_index in genome_to_mrna_map
            ):
                aso_df.at[index, IN_CODING_REGION] = True
                cds_idx = genome_to_mrna_map[genome_corrected_index]
                for fs in FLANK_SIZES_CDS:
                    start = max(0, cds_idx - fs)
                    end = min(len(cds_seq), cds_idx + len(sense) + fs)
                    local_seq = cds_seq[start:end]
                    aso_df.at[index, f"local_coding_region_around_ASO_{fs}"] = _to_str_seq(local_seq)

    CDS_WINDOWS = FLANK_SIZES_CDS

    # Loop over each flank window size
    for flank in CDS_WINDOWS:
        local_col = f"local_coding_region_around_ASO_{flank}"
        is_local_flag_col = f"region_is_local_{flank}"

        # Create the binary flag: 1 if local exists, 0 otherwise
        aso_df[is_local_flag_col] = aso_df[local_col].apply(
            lambda x: isinstance(x, str) and x.strip() != ""
        ).astype(int)

    for flank in CDS_WINDOWS:
        local_col = f"local_coding_region_around_ASO_{flank}"
        CAI_col = f"CAI_score_{flank}_CDS"
        aso_df[CAI_col] = (
            aso_df[local_col].astype(str).apply(lambda s: calc_CAI(s, weights_flat))
        )
    aso_df["CAI_score_global_CDS"] = (
        aso_df["cds_sequence"].astype(str).apply(lambda s: calc_CAI(s, weights_flat))
    )
