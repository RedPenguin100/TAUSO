from src.tauso.off_target.Roni.off_target_pipeline.run_transcriptome_pipeline import run_transcriptome_pipeline
import os
from pyprojroot import here

PROJECT_ROOT = here()
DATA_DIR = os.path.join(PROJECT_ROOT, "src", "tauso", "off_target", "Roni", "data")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "src", "tauso", "off_target", "Roni", "outputs")

fasta_path = os.path.join(DATA_DIR, "GRCh38.p13.genome.fa")
gtf_path = os.path.join(DATA_DIR, "gencode.v34.chr_patch_hapl_scaff.annotation.gtf")
exp_path = os.path.join(DATA_DIR, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
mut_path = os.path.join(DATA_DIR, "OmicsSomaticMutations.csv")

#cell_line_lst = ["ACH-000681"]
cell_line_lst = ["general", "ACH-001328", "ACH-000463", "ACH-001188", "ACH-001086", "ACH-000739", "ACH-000232"]

def generate_cell_lines_transcriptomes(save_csv=False):
    return run_transcriptome_pipeline(fasta_path, gtf_path, exp_path, mut_path, OUTPUT_DIR, cell_line_lst, save_csv=save_csv)