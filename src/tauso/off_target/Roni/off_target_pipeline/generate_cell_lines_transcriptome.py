import pandas as pd
import os
from pyprojroot import here
from src.tauso.off_target.Roni.off_target_pipeline.off_target_functions import parse_gtf, parse_fasta
from src.tauso.off_target.Roni.off_target_pipeline.get_premRNA_sequences import enrich_expression_data_with_sequence
from src.tauso.off_target.Roni.off_target_pipeline.mutate_cell_line_transcriptome import get_expression_of_cell_line, get_mutations_of_cell_line, mutate_transcriptome

PROJECT_ROOT = here()
DATA_DIR = os.path.join(PROJECT_ROOT, "src", "tauso", "off_target", "Roni", "data")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "src", "tauso", "off_target", "Roni", "outputs")

fasta_path = os.path.join(DATA_DIR, "GRCh38.p13.genome.fa")
gtf_path = os.path.join(DATA_DIR, "gencode.v34.chr_patch_hapl_scaff.annotation.gtf")
exp_path = os.path.join(DATA_DIR, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
mut_path = os.path.join(DATA_DIR, "OmicsSomaticMutations.csv")

cell_line_lst = ["ACH-000681"]
save_csv = True

def generate_cell_lines_transcriptomes(
    fasta_fl_path, gtf_fl_path, expression_fl_path, mutation_fl_path,
    output_dir, cell_line_list, save_csv=False):
    """
    Generate cell line–specific transcriptomes by integrating genomic sequences,
    annotations, expression profiles, and somatic mutations.

    Parameters
    ----------
    fasta_fl_path : str
        Path to the reference human genome FASTA file.
    gtf_fl_path : str
        Path to the GTF file containing transcript annotations.
    expression_fl_path : str
        Path to the DepMap CSV file with gene expression data.
    mutation_fl_path : str
        Path to the DepMap CSV file with somatic mutation data.
    output_dir : str
        Directory for cached files and optional outputs.
    cell_line_list : list of str
        DepMap cell line identifiers (e.g., ACH-######).
    save_csv : bool, optional
        If True, save each generated transcriptome as a CSV file.

    Returns
    -------
    dict
        A dictionary mapping each cell line to a list of pandas DataFrames:
        [expression_data, mutation_data, mutated_transcriptome].
    """
    transcriptomes = {}
    sequences = parse_fasta(fasta_fl_path, os.path.join(output_dir, "fasta_sequences.pkl"))
    annotations = parse_gtf(gtf_fl_path, os.path.join(output_dir, "gtf_annotations.pkl"))

    for cell_line in cell_line_list:

        exp_data = get_expression_of_cell_line(cell_line, expression_fl_path, output_dir)
        mut_data = get_mutations_of_cell_line(cell_line, mutation_fl_path, output_dir)

        transcriptomes[cell_line] = [exp_data, mut_data]

        transcriptomes = enrich_expression_data_with_sequence(transcriptomes, cell_line, sequences, annotations)

        exp_df = transcriptomes[cell_line][0]
        exp_df["Original Transcript Sequence"] = (
            exp_df["Original Transcript Sequence"].astype(str).str.replace("T", "U"))
        mut_df = transcriptomes[cell_line][1]

        mut_exp_data = mutate_transcriptome(exp_df, mut_df, annotations)
        transcriptomes[cell_line].append(mut_exp_data)

        if save_csv:
            df = pd.DataFrame(transcriptomes[cell_line][-1])
            df.to_csv(os.path.join(OUTPUT_DIR, cell_line + "_transcriptome.csv"), index=False)

    return transcriptomes


generate_cell_lines_transcriptomes(fasta_path, gtf_path, exp_path, mut_path, OUTPUT_DIR, cell_line_lst, save_csv=save_csv)