import os
import pandas as pd
from src.tauso.off_target.Roni.off_target_pipeline.off_target_functions import parse_gtf, parse_fasta
from src.tauso.off_target.Roni.off_target_pipeline.get_premRNA_sequences import enrich_expression_data_with_sequence, general_exp_data
from src.tauso.off_target.Roni.off_target_pipeline.mutate_cell_line_transcriptome import get_expression_of_cell_line, get_mutations_of_cell_line, mutate_transcriptome


def run_transcriptome_pipeline(
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

        if cell_line == "general":
            exp_data = general_exp_data(expression_fl_path)
            mut_data = exp_data
        else:
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
            df.to_csv(os.path.join(output_dir, cell_line + "_transcriptome.csv"), index=False)

    return transcriptomes