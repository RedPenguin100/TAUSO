import os
import pandas as pd
from off_target_functions import dna_to_rna_reverse_complement, normalize_chrom, name2accession
import pickle
from pathlib import Path

"""
This code extracts the mRNA sequences 
"""
# Get the directory of the current script
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# Full path to cell line data folder
DATA_DIR = os.path.join(PROJECT_ROOT, "scripts", "data_genertion", "cell_line_expression")


# cell_line2data = {
#     "ACH-001328": [os.path.join(DATA_DIR, "ACH-001328_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-001328_mutations.csv")],
#     "ACH-000463": [os.path.join(DATA_DIR, "ACH-000463_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-000463_mutations.csv")],
#     "ACH-001188": [os.path.join(DATA_DIR, "ACH-001188_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-001188_mutations.csv")],
#     "ACH-001086": [os.path.join(DATA_DIR, "ACH-001086_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-001086_mutations.csv")],
#     "ACH-000739": [os.path.join(DATA_DIR, "ACH-000739_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-000739_mutations.csv")],
#     "ACH-000232": [os.path.join(DATA_DIR, "ACH-000232_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-000232_mutations.csv")]
# }

cell_line2data = {"ACH-000681": [os.path.join(DATA_DIR, "ACH-000681_transcriptome.csv"), os.path.join(DATA_DIR, "ACH-000681_mutations.csv")]}


# Reference file paths
fasta_path = os.path.join(DATA_DIR, "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa")
gtf_path = os.path.join(DATA_DIR, "gencode.v48.chr_patch_hapl_scaff.annotation.gtf")

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

def get_relevant_exp_data(cell_line):
    exp_fl = cell_line2data[cell_line][0]
    mut_fl = cell_line2data[cell_line][1]
    if isinstance(cell_line2data[cell_line][0], pd.DataFrame):
        exp_data = exp_fl
    else:
        exp_data = pd.read_csv(cell_line2data[cell_line][0], usecols=[0, 1, 2, 3])

    if isinstance(cell_line2data[cell_line][1], pd.DataFrame):
        mut_data = mut_fl
    else:
        mut_data = pd.read_csv(cell_line2data[cell_line][1])

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

"""

def final_func_premrna(cell_line_dict, genome_pkl, annotation_pkl):
    # Load from pre-generated pickle files
    with open(annotation_pkl, 'rb') as f:
        gtf_dict = pickle.load(f)
    print('Read GTF.')
    with open(genome_pkl, 'rb') as f:
        print("About to read FASTA")
        fasta_dict = pickle.load(f)
        print("Successfully read FASTA")

    print(f'fasta keys: {list(fasta_dict.keys())[:5]}')  # limit to 5 keys
    print('Read FASTA.')

    print(fasta_dict.keys())

    all_results = {}

    for cell_line, _ in cell_line_dict.items():
        print(f"Processing {cell_line}...")

        # Load expression and mutation data
        exp_data, mut_data = get_relevant_exp_data(cell_line)

        original_seqs = []
        for idx, row in exp_data.iterrows():
            transcript_id_full = row["Transcript_ID"]
            transcript_id = transcript_id_full.split('.')[0]  # remove version

            transcript_info = gtf_dict.get(transcript_id_full) or gtf_dict.get(transcript_id)
            #print(transcript_info)

            if not transcript_info:
                original_seqs.append(None)
                #print(f"Transcript {transcript_id_full} not found in GTF.")
                continue

            chrom = transcript_info["chrom"]
            strand = transcript_info["strand"]
            start = transcript_info["start"]
            end = transcript_info["end"]

            #norm_chrom = normalize_chrom(chrom)
            if chrom not in name2accession:
                original_seqs.append(None)
                print(f"Skipping transcript {transcript_id_full}: chromosome {chrom} not in accession map.")
                continue
            norm_chrom = name2accession[chrom]
            try:
                seq = fasta_dict[norm_chrom][start - 1:end]  # extract using 0-based indexing
                if strand == "-":
                    seq = dna_to_rna_reverse_complement(seq)
                else:
                    seq = seq.replace("T", "U")
                original_seqs.append(seq)
            except Exception as e:
                print(norm_chrom in fasta_dict)  # Should return True or False
                print("Trying chrom:", norm_chrom)
                print(f"Error extracting for {transcript_id_full}: {e}")
                original_seqs.append(None)

        exp_data["Original Transcript Sequence"] = original_seqs
        exp_data["Mutated Transcript Sequence"] = None  # for future use

        # Save output
        output_path = f"{cell_line}_transcriptome_premRNA.csv"
        exp_data.to_csv(output_path, index=False)
        print(f"Saved {output_path}")

        all_results[cell_line] = exp_data

    return all_results
"""

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
        exp_data, mut_data = get_relevant_exp_data(cell_line)

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



# genome_path = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/chr_1_4.pkl'
# annotation_path = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/_gtf_annotations.pkl'
# final_func_premrna(cell_line2data, genome_path, annotation_path)
#
# print('test')