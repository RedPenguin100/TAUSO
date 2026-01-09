import os
import pandas as pd
from off_target_functions import dna_to_rna_reverse_complement, normalize_chrom, name2accession
from gtf_functions import build_gene_to_transcript_index, select_best_transcript
import pickle
from pathlib import Path

"""
This code extracts the mRNA sequences 
"""
# Get the directory of the current script
#PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# Full path to cell line data folder
#DATA_DIR = os.path.join(PROJECT_ROOT, "scripts", "data_genertion", "cell_line_expression")


# cell_line2data = {
#     "ACH-001328": [os.path.join(DATA_DIR, "ACH-001328_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-001328_mutations.csv")],
#     "ACH-000463": [os.path.join(DATA_DIR, "ACH-000463_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-000463_mutations.csv")],
#     "ACH-001188": [os.path.join(DATA_DIR, "ACH-001188_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-001188_mutations.csv")],
#     "ACH-001086": [os.path.join(DATA_DIR, "ACH-001086_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-001086_mutations.csv")],
#     "ACH-000739": [os.path.join(DATA_DIR, "ACH-000739_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-000739_mutations.csv")],
#     "ACH-000232": [os.path.join(DATA_DIR, "ACH-000232_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-000232_mutations.csv")]
# }

#cell_line2data = {"ACH-000681": [os.path.join(DATA_DIR, "ACH-000681_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-000681_mutations.csv")]}


# Reference file paths
#fasta_path = os.path.join(DATA_DIR, "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa")
#gtf_path = os.path.join(DATA_DIR, "gencode.v48.chr_patch_hapl_scaff.annotation.gtf")

def load_pickle_if_needed(obj, name="object"):
    """
    If obj is a dict, return it.
    If obj is a path (str or Path), load pickle and return it.
    """
    if isinstance(obj, dict):
        return obj

    if isinstance(obj, (str, Path)):
        with open(obj, "rb") as f:
            loaded = pickle.load(f)
        print(f"Loaded {name} from pickle: {obj}")
        return loaded

    raise TypeError(
        f"{name} must be either a dict or a path to a pickle file, "
        f"got {type(obj)}"
    )

def get_relevant_exp_data(cell_line_dict, cell_line):
    exp_fl = cell_line_dict[cell_line][0]
    mut_fl = cell_line_dict[cell_line][1]
    print(exp_fl, mut_fl)
    if isinstance(cell_line_dict[cell_line][0], pd.DataFrame):
        exp_data = exp_fl
        print("hi")
    else:
        print("hmmm")
        exp_data = pd.read_csv(cell_line_dict[cell_line][0], usecols=[0, 1, 2, 3])

    if isinstance(cell_line_dict[cell_line][1], pd.DataFrame):
        mut_data = mut_fl
    else:
        mut_data = pd.read_csv(cell_line_dict[cell_line][1])

    exp_data = exp_data[exp_data["Transcript_ID"].str.contains("ENST", na=False)]
    ###
    #exp_data = exp_data.head(15)
    return exp_data, mut_data


def get_premrna_coords(annotation):
    coords = []
    for feature in ["5utrs", "exons", "introns", "CDS", "3utrs"]:
        coords.extend(annotation.get(feature, []))
    if coords:
        return min(s for s, _ in coords), max(e for _, e in coords)
    return None, None

def add_original_sequence(exp_data, fasta_dict):
    def fetch_sequence(row):
        tid = row["Transcript_ID"].split('.')[0]
        return fasta_dict.get(tid)
    exp_data["Original Transcript Sequence"] = exp_data.apply(fetch_sequence, axis=1)
    return exp_data



def final_func_premrna(cell_line_dict, genome, annotation):
    # Load dicts if needed
    gtf_dict = load_pickle_if_needed(annotation, name="GTF annotations")
    fasta_dict = load_pickle_if_needed(genome, name="Genome FASTA")

    print('GTF and FASTA ready.')
    print(f'fasta keys (sample): {list(fasta_dict.keys())[:5]}')

    all_results = {}

    for cell_line, _ in cell_line_dict.items():
        print(f"Processing {cell_line}...")

        # Load expression and mutation data
        exp_data, mut_data = get_relevant_exp_data(cell_line_dict, cell_line)

        original_seqs = []

        for _, row in exp_data.iterrows():
            transcript_id_full = row["Transcript_ID"]
            transcript_id = transcript_id_full.split('.')[0]

            transcript_info = (
                gtf_dict.get(transcript_id_full)
                or gtf_dict.get(transcript_id)
            )

            if not transcript_info:
                original_seqs.append(None)
                continue

            chrom = transcript_info["chrom"]
            strand = transcript_info["strand"]
            start = transcript_info["start"]
            end = transcript_info["end"]

            if chrom not in name2accession:
                original_seqs.append(None)
                print(
                    f"Skipping {transcript_id_full}: "
                    f"chromosome {chrom} not in accession map."
                )
                continue

            norm_chrom = name2accession[chrom]

            try:
                seq = fasta_dict[norm_chrom][start - 1:end]

                if strand == "-":
                    seq = dna_to_rna_reverse_complement(seq)
                else:
                    seq = seq.replace("T", "U")

                original_seqs.append(seq)

            except Exception as e:
                print(f"Error extracting {transcript_id_full}: {e}")
                original_seqs.append(None)

        exp_data["Original Transcript Sequence"] = original_seqs
        exp_data["Mutated Transcript Sequence"] = None

        output_path = f"{cell_line}_transcriptome_premRNA.csv"
        exp_data.to_csv(output_path, index=False)
        print(f"Saved {output_path}")

        all_results[cell_line] = exp_data

    return all_results


def enrich_expression_data_with_sequence(cell_line_dict, curr_cell_line, fasta_dict, gtf_annotations):
    """
    Modifies cell_line_dict in-place.
    Adds 'Transcript_ID' and 'Original Transcript Sequence' to the exp_data DataFrame.

    Args:
        cell_line_dict: { 'ACH-###': [exp_df, mut_df] }
        fasta_dict: { 'chr1': 'ATCG...' } (Keys must match GTF chrom names)
        gtf_annotations: Dictionary from parse_gtf (MUST include 'gene_name')
    """

    # 1. Build the index once to speed up lookups
    print("Building Gene -> Transcript index...")
    gene_to_tids = build_gene_to_transcript_index(gtf_annotations)

    exp_data = cell_line_dict[curr_cell_line][0]
    mut_data = cell_line_dict[curr_cell_line][1]
    # 2. Iterate through each cell line

    print(f"Processing {curr_cell_line}...")

    # Prepare lists to collect new column data
    transcript_ids = []
    sequences = []

    # 3. Iterate through genes in the expression dataframe
    # Assuming format: "MT-ATP8 (4509)"
    for full_gene_str in exp_data['Gene']:

        # Extract pure symbol "MT-ATP8" from "MT-ATP8 (4509)"
        gene_symbol = full_gene_str.split(' ')[0]

        # Find best transcript
        tid = select_best_transcript(gene_symbol, gene_to_tids, gtf_annotations)

        seq_str = None

        if tid:
            # Get coordinates
            annot = gtf_annotations[tid]
            chrom = annot['chrom']
            start = annot['start']  # 1-based from GTF
            end = annot['end']  # 1-based from GTF
            strand = annot['strand']

            # Check if chromosome exists in FASTA (handle 'chr' prefix mismatches if needed)
            if chrom in fasta_dict:
                # Python slicing is 0-based: [start-1 : end]
                # We grab the whole genomic span (pre-mRNA)
                raw_seq = fasta_dict[chrom][start - 1: end]

                # Handle Strand
                if strand == '-':
                    seq_str = dna_to_rna_reverse_complement(raw_seq)
                else:
                        seq_str = raw_seq
            else:
                print(f"Warning: Chromosome {chrom} not found in Fasta.")

        transcript_ids.append(tid)
        sequences.append(seq_str)

        # 4. Assign new columns
    exp_data['Transcript_ID'] = transcript_ids
    exp_data['Original Transcript Sequence'] = sequences

    # Update the dictionary (redundant since lists are mutable, but clear)
    cell_line_dict[curr_cell_line][0] = exp_data

    print("Enrichment complete.")
    return cell_line_dict

# genome_path = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/chr_1_4.pkl'
# annotation_path = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/_gtf_annotations.pkl'
# final_func_premrna(cell_line2data, genome_path, annotation_path)
#
# print('test')

def general_exp_data(EXP_path):
    # load expression matrix
    exp_data = pd.read_csv(EXP_path)

    # keep only gene columns (gene name + number in parentheses)
    gene_cols = exp_data.columns[exp_data.columns.str.contains(r"\(")]
    exp_data = exp_data[gene_cols]

    # average across cell lines → one value per gene
    mean_exp = exp_data.mean(axis=0)

    # convert Series → DataFrame with two columns
    mean_exp_data = mean_exp.reset_index()
    mean_exp_data.columns = ["Gene", "general_expression_norm"]
    mean_exp_data = mean_exp_data.sort_values("general_expression_norm", ascending=False)

    # back-transform to TPM (assuming log2(TPM))
    mean_exp_data["expression_TPM"] = 2 ** mean_exp_data["general_expression_norm"]

    return mean_exp_data