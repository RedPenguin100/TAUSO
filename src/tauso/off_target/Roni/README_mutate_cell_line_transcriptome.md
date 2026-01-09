# Transcriptome Mutation Integration Pipeline

This script processes human cell line data to generate mutated transcriptome profiles by integrating:
- Normalized expression values (TPM),
- Somatic mutations (SNVs, deletions, insertions, etc.),
- Genomic annotations (GTF), and
- Reference transcript sequences (FASTA).

It outputs the top `N` most highly expressed transcripts (per cell line), with mutations applied directly to their sequences where possible.

---

## 📁 Inputs

| Input | Description |
|------|-------------|
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | Normalized gene expression matrix (`log2(TPM + 1)`) with genes as columns and cell lines as rows |
| `OmicsSomaticMutations.csv` | Somatic mutation file containing `ModelID`, `DNAChange`, and `VariantType` |
| `gencode.v48.chr_patch_hapl_scaff.annotation.gtf` | GTF annotation file containing `exon` and `CDS` features |
| `Homo_sapiens.GRCh38.cdna.all.fa` | FASTA file of cDNA transcript sequences from Ensembl |
| `ACH-*` | DepMap ID of the cell line to be analyzed. Represents a unique identifier assigned by the DepMap project. |
| `top_n` | Integer specifying the number of top-expressed genes (ranked by TPM) to include in the analysis for each cell line. |


## 🚀 How It Works

1. **Expression Parsing**  
   Extracts the top `N` most highly expressed protein-coding genes for a given cell line, calculates TPM, and sorts descending.

2. **Mutation Parsing**  
   Filters the mutation dataset for the given cell line, extracts standardized mutation fields (e.g. `SNV`, `deletion`, `insertion`, etc.).

3. **Transcript Annotation**  
   Maps gene symbols to transcript IDs and sequences using the FASTA file. Only the first match for each gene is used.

4. **CDS Localization**  
   Computes the offset from the transcript start to the CDS start using flattened exon structures in the GTF file.

5. **Sequence Mutation**  
   Applies mutations (where possible) directly to the transcript sequences based on position and type. Unsupported or mismatched cases are logged and skipped.

6. **Result Export**  
   Outputs a CSV for each cell line with the top `N` transcripts, their original and mutated sequences, expression levels, and transcript IDs.

---

## 📦 Output Format

Each resulting CSV file (e.g., `ACH-001328_transcriptome_top500.csv`) contains:

| Column | Description |
|--------|-------------|
| `Gene` | Gene name |
| `ACH-XXXX_expression_norm` | Original log-normalized expression |
| `expression_TPM` | Transformed TPM expression |
| `Transcript_ID` | Ensembl transcript ID |
| `Original Transcript sequence` | cDNA sequence from FASTA |
| `Mutated Transcript sequence` | Modified sequence with mutations applied, if applicable |

---

## 🧪 Example Cell Lines Processed

- A431 (`ACH-001328`)
- NCI-H460 (`ACH-000463`)
- SH-SY5Y (`ACH-001188`)
- HeLa (`ACH-001086`)
- HepG2 (`ACH-000739`)
- U-251MG (`ACH-000232`)

---

## 🛠️ Requirements

- Python 3.7+
- `pandas`
- `biopython`

Install dependencies via:

```bash
pip install pandas biopython

## 📜 License

This project is part of the ASO Design pipeline and is for research use only.
