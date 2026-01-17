import pandas as pd
import os
from pyprojroot import here

from src.tauso.off_target.Roni.off_target_pipeline.off_target_feature import run_off_target_pipeline

PROJECT_ROOT = here()
DATA_DIR = os.path.join(PROJECT_ROOT, "src", "tauso", "off_target", "Roni", "data")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "src", "tauso", "off_target", "Roni", "outputs")

ASO_DATA_DIR = os.path.join(PROJECT_ROOT, "notebooks", "data")
aso_data_path = os.path.join(ASO_DATA_DIR, "data_asoptimizer_updated.csv")
aso_df = pd.read_csv(aso_data_path)

cellline_name_to_depmap = {
    'A431': 'ACH-001328',
    'A-431': 'ACH-001328',
    'NCI-H460': 'ACH-000463',
    'SH_SY5Y': 'ACH-001188',
    'SH-SY5Y': 'ACH-001188',
    'HeLa': 'ACH-001086',
    'Hela': 'ACH-001086',
    'HepG2': 'ACH-000739',
    'U-251MG': 'ACH-000232',
    'U251': 'ACH-000232'
}

def generate_off_target_vectors(transcriptomes, top_ns, cutoffs, method=0):
    """
    transcriptomes is a dictionary with cell line depmap ids as keys and a list of expression data, mutation data and transcriptome data as values
    top_ns is a list of top n transcripts to consider as potential off-targets
    cutoffs is a list of cutoff values to consider as meaningful off-targets
    """
    return run_off_target_pipeline(
        ASO_df=aso_df,
        cell_line2data=transcriptomes ,
        output_dir=OUTPUT_DIR, top_n_list=top_ns,
        cutoff_list=cutoffs, method=method,
        name_to_depmap=cellline_name_to_depmap
    )