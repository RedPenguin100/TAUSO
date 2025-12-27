import pandas as pd
from Bio import SeqIO
import re
import csv
import os

"""
This code might not work if RAM has less than 8 GB (relevant for working with pre-mRNA (long) sequences)
"""


def get_expression_of_cell_line(cell_line, expression_file, output_dir):
    """
    Reads a large CSV line-by-line to find a specific cell_line (ModelID),
    ignores metadata columns, and returns a clean Gene Expression DataFrame.
    """
    # 1. Normalize input to match file format (e.g., 'ach-001' -> 'ACH-001')
    target_id = cell_line.upper()

    target_row = None
    header = None

    # We use csv.reader which handles large files without loading everything to RAM
    with open(expression_file, 'r') as f:
        reader = csv.reader(f)

        # Read the header first
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("File appears to be empty.")

        # --- A. Find the Metadata Columns ---
        # Find which column index holds the 'ModelID'
        try:
            model_idx = header.index("ModelID")
        except ValueError:
            raise ValueError("Column 'ModelID' was not found in the CSV header.")

        # --- B. Find where Genes Start ---
        # Heuristic: Metadata cols are plain text, Genes usually look like "Gene (ID)"
        # We look for the first column that contains parentheses "(".
        # Based on your file, this should be index 6 (7th column).
        gene_start_idx = 0
        for i, col_name in enumerate(header):
            if "(" in col_name and ")" in col_name:
                gene_start_idx = i
                break

        # Fallback: if no parentheses found, assume standard DepMap layout (skip first 6 cols)
        if gene_start_idx == 0:
            gene_start_idx = 6

        # --- C. Stream the file to find the row ---
        for row in reader:
            # Check the specific ModelID column
            if len(row) > model_idx and row[model_idx] == target_id:
                target_row = row
                break

    if target_row is None:
        raise ValueError(f"Cell line '{cell_line}' not found in file.")

    # --- D. Extract only Gene Data ---
    # Slice the lists to keep only the gene columns
    gene_names = header[gene_start_idx:]
    expression_values = target_row[gene_start_idx:]

    # Create the Series (Genes as index)
    # coerce errors='coerce' turns non-numeric strings to NaN, just in case
    expression_series = pd.to_numeric(pd.Series(expression_values, index=gene_names), errors='coerce')
    expression_series.name = f'{cell_line}_expression_norm'

    # Convert to DataFrame
    expression_df = expression_series.reset_index()
    expression_df.columns = ['Gene', f'{cell_line}_expression_norm']

    # --- E. Calculate TPM (Reversing the log2(x+1)) ---
    # Note: Ensure your input data is indeed log2(TPM+1) before doing this
    expression_df['expression_TPM'] = (2 ** expression_df[f'{cell_line}_expression_norm']) - 1

    # --- F. Filter and Sort ---
    expression_df = expression_df[expression_df['expression_TPM'] > 0]
    expression_df = expression_df.sort_values(by='expression_TPM', ascending=False).reset_index(drop=True)

    # Save and return
    output_filename = f'{cell_line}_expression.csv'
    #expression_df.to_csv(output_filename, index=False)
    expression_df.to_csv(os.path.join(output_dir, output_filename), index=False)

    print(f"Successfully processed {target_id}. Found {len(expression_df)} expressed genes.")
    return expression_df


def get_mutations_of_cell_line(cell_line, mutation_file, output_dir):
    """
    Creates a DataFrame that contains all mutations of the cell line
    """
    necessary_cols = ['ModelID', 'Ref', 'Alt', 'VariantType', 'VariantInfo', 'DNAChange', 'HugoSymbol']

    # read data
    mutation_data = pd.read_csv(mutation_file, usecols=necessary_cols)

    # filter by cell line
    cell_line_mutation = mutation_data[mutation_data['ModelID'] == cell_line]
    #cell_line_mutation.to_csv(f'{cell_line}_mutations.csv', index=False)
    cell_line_mutation.to_csv(os.path.join(output_dir, f'{cell_line}_mutations.csv'), index=False)
    print(f"Obtained mutation data for {cell_line}")
    return cell_line_mutation


def get_record_by_hugo(indexed_fasta, hugo_symbol):
    for record in indexed_fasta:
        header = record.description
        if f"gene_symbol:{hugo_symbol}" in header:
            return record
    return None


def extract_cds_from_gtf(gtf_path):
    # Unused in final_func
    cds_rows = []

    with open(gtf_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # skip header
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            if feature != 'CDS':
                continue

            # extract transcript_id from attributes field
            attr_dict = {}
            for attr in attributes.strip().split(';'):
                if attr.strip():
                    key_value = attr.strip().split(' ', 1)
                    if len(key_value) == 2:
                        key, value = key_value
                        attr_dict[key] = value.strip('"')

            transcript_id = attr_dict.get('transcript_id')
            if not transcript_id:
                continue

            cds_rows.append({
                'transcript_id': transcript_id,
                'start': int(start),
                'end': int(end),
                'strand': strand
            })

    return pd.DataFrame(cds_rows)


def prepare_gtf_dataframe(gtf_file):
    """
    Load and preprocess GTF file, keeping only 'exon' and 'CDS' entries with parsed transcript IDs.
    """
    col_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    df = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=col_names)

    # Filter only 'exon' and 'CDS'
    df = df[df['feature'].isin(['exon', 'CDS'])]

    # Extract transcript_id from attributes
    df['transcript_id'] = df['attribute'].str.extract(r'transcript_id "([^"]+)"')

    # Drop entries without transcript_id (should be rare or malformed)
    df = df.dropna(subset=['transcript_id'])

    return df


# def find_shift(tid_dict, mut_idx):
#     mut_idx = int(mut_idx)
#     if mut_idx <= 0:
#         print(f"Warning: mut_idx must be 1-based (>= 1), got {mut_idx}")
#         return None
#
#     accu_len = 0
#     cds_list = tid_dict['cds']
#
#     for s, e in cds_list:
#         seg_len = e - s + 1
#         if accu_len + seg_len >= mut_idx:
#             offset = mut_idx - accu_len  # 0-based offset
#             genomic_pos = s + offset
#             shift = genomic_pos - tid_dict['start'] - mut_idx
#             return shift
#         accu_len += seg_len
#
#     print(f"Warning: mut_idx {mut_idx} is out of bounds for total CDS length.")
#     return None


# def find_shift(tid_dict, mut_idx):
#     mut_idx = int(mut_idx)
#     if mut_idx <= 0:
#         print(f"Warning: mut_idx must be 1-based (>= 1), got {mut_idx}")
#         return None
#
#     cds_list = tid_dict['cds']
#     strand = tid_dict['strand']
#     start = tid_dict['start']
#     end = tid_dict['end']
#
#     # For minus strand, reverse the CDS list and adjust the order
#     if strand == '-':
#         cds_list = sorted(cds_list, reverse=True)
#
#     accu_len = 0
#     for s, e in cds_list:
#         seg_len = e - s + 1
#         if accu_len + seg_len >= mut_idx:
#             offset = mut_idx - accu_len  # offset into this CDS segment
#
#             if strand == '+':
#                 genomic_pos = s + offset
#                 shift = genomic_pos - start - mut_idx
#
#             else:  # strand == '-'
#                 genomic_pos = e - offset
#                 shift = genomic_pos - start - mut_idx
#
#             return shift
#
#         accu_len += seg_len
#
#     print(f"Warning: mut_idx {mut_idx} is out of bounds for total CDS length.")
#     return None

def find_shift(tid_dict, mut_idx):
    """
    This function analyzes tid_dict, which has all the exons' coordinates.
    It maps the mutation index relative to the beginning of CDS to its location
    within the entire pre-mRNA sequence.
    """

    mut_idx = int(mut_idx)
    if mut_idx == 0:
        print(f"Warning: mut_idx must be 1-based (>= 1), got {mut_idx}")
        return None

    cds_list = tid_dict['cds']
    strand = tid_dict['strand']
    start = tid_dict['start']
    end = tid_dict['end']

    # For minus strand, reverse the CDS list and adjust the order
    if strand == '-':
        cds_list = sorted(cds_list, reverse=True)

    accu_len = 0
    for s, e in cds_list:

        curr_s = s
        curr_e = e

        seg_len = e - s + 1
        if accu_len + seg_len >= mut_idx:
            offset = mut_idx - accu_len  # offset into this CDS segment

            if strand == '+':
                genomic_pos = s + offset
                shift = genomic_pos - start - mut_idx

            else:  # strand == '-'
                genomic_pos = e - offset
                shift = end - genomic_pos - mut_idx

            return shift

        accu_len += seg_len

    if accu_len < mut_idx:
        offset = mut_idx - accu_len + 1

        if strand == '+':
            genomic_pos = curr_e + offset
            shift = genomic_pos - start - mut_idx

        else:  # strand == '-'
            genomic_pos = curr_s - offset
            shift = end - genomic_pos - mut_idx

        return shift

    print(f"Warning: mut_idx {mut_idx} is out of bounds for total CDS length.")
    return None


def parse_cdna_jumps(cdna_str): # Robust but was originally made for splicing mutations
    # Remove leading "c."
    if cdna_str.startswith("c."):
        cdna_str = cdna_str[2:]

    # Find positions (start and optional end)
    if "_" in cdna_str:
        start_part, end_part = cdna_str.split("_", 1)
    else:
        start_part = cdna_str
        end_part = start_part  # same if no range given

    # Helper to extract coordinate and jump
    def extract_parts(part):
        match = re.match(r"(\d+)([+-]\d+)?", part)
        if match:
            coord = int(match.group(1))
            jump = int(match.group(2)) if match.group(2) else 0
            return coord, jump
        return None, None  # if parsing fails

    start, start_jump = extract_parts(start_part)
    end, end_jump = extract_parts(end_part)

    return start, start_jump, end, end_jump


def safe_int(x):
    try:
        return int(x)
    except (ValueError, TypeError):
        return None


def get_distance_to_start_codon_from_df(df, transcript_id):
    """
    ========================== Unused in final_func, relevant for mature mRNA sequences only ==========================
    Compute CDS start position (in transcript/mRNA coordinates) by flattening exons and locating CDS start.
    """
    # get exons and CDS for this transcript
    exons = df[(df['transcript_id'] == transcript_id) & (df['feature'] == 'exon')]
    cds = df[(df['transcript_id'] == transcript_id) & (df['feature'] == 'CDS')]

    if exons.empty or cds.empty:
        raise ValueError(f"Transcript {transcript_id} not found or missing exons/CDS in GTF.")

    strand = exons.iloc[0]['strand']

    # sort exons and CDS by genomic position (according to strand)
    if strand == '+':
        exons = exons.sort_values(by='start')
        cds_start = cds['start'].min()
    else:
        exons = exons.sort_values(by='start', ascending=False)
        cds_start = cds['end'].max()

    # compute CDS start position in transcript (mRNA) coordinates
    transcript_pos = 0
    for _, exon in exons.iterrows():
        exon_start = exon['start']
        exon_end = exon['end']
        exon_len = exon_end - exon_start + 1

        if strand == '+':
            if exon_start <= cds_start <= exon_end:
                return transcript_pos + (cds_start - exon_start)
        else:  # reverse strand
            if exon_start <= cds_start <= exon_end:
                return transcript_pos + (exon_end - cds_start)

        transcript_pos += exon_len

    raise ValueError(f"Could not find CDS start position in exons for transcript {transcript_id}")#


# returns a dict that contains the orders for mutation

# def mutation_dict(mut_row):
#     # Skip if DNAChange is missing or empty
#     if pd.isna(mut_row.get('DNAChange')) or mut_row['DNAChange'].strip() == '':
#         return None
#
#     try:
#         # Split into ID and mutation string
#         transcript_id, mutation = mut_row['DNAChange'].split(':')
#     except ValueError:
#         return None  # if it can't unpack into two parts
#
#     if mut_row['VariantType'] == 'SNV':
#         match = re.search(r'\d+', mutation)
#         if match:
#             start, end = match.span()
#             location = mutation[start:end]
#             from_nt, to_nt = mutation[end:].split('>')
#             result = {'id': transcript_id, 'start': location, 'end': location, 'ref': from_nt, 'alt': to_nt, 'type': mut_row['VariantType']}
#         else:
#             result = {'id': transcript_id, 'start': None, 'end': None, 'ref': None, 'alt': None, 'type': None}
#
#     elif mut_row['VariantType'] == 'deletion':
#         numbers = re.findall(r'\d+', mutation)
#         if len(numbers) == 1:
#             result = {'id': transcript_id, 'start': numbers[0], 'end': numbers[0], 'type': mut_row['VariantType']}
#         elif len(numbers) >= 2:
#             result = {'id': transcript_id, 'start': numbers[0], 'end': numbers[1], 'type': mut_row['VariantType']}
#         else:
#             result = {'id': transcript_id, 'start': None, 'end': None, 'type': None}
#
#     elif mut_row['VariantType'] == 'substitution':
#         numbers = re.findall(r'\d+', mutation)
#         if len(numbers) == 1:
#             result = {'id': transcript_id, 'start': numbers[0], 'end': numbers[0]}
#         elif len(numbers) >= 2:
#             result = {'id': transcript_id, 'start': numbers[0], 'end': numbers[1]}
#         else:
#             result = {'id': transcript_id, 'start': None, 'end': None, 'type': None}
#
#         if 'delins' in mutation:
#             result['type'] = 'delins'
#             result['subs'] = mutation.split('delins')[-1]
#         elif 'inv' in mutation:
#             result['type'] = 'inversion'
#
#     elif mut_row['VariantType'] == 'insertion':
#         numbers = re.findall(r'\d+', mutation)
#         if len(numbers) == 1:
#             result = {'id': transcript_id, 'start': numbers[0], 'end': numbers[0]}
#         elif len(numbers) >= 2:
#             result = {'id': transcript_id, 'start': numbers[0], 'end': numbers[1]}
#         else:
#             result = {'id': transcript_id, 'start': None, 'end': None, 'type': None}
#
#         if 'ins' in mutation:
#             result['type'] = 'insertion'
#             result['ins'] = mutation.split('ins')[-1]
#         elif 'dup' in mutation:
#             result['type'] = 'duplication'
#     else:
#         return None
#
#     result['info'] = mut_row.get('VariantInfo')
#
#     if isinstance(result['info'], str) and (
#             "splice_donor_variant" in result['info'] or "splice_acceptor_variant" in result['info']
#     ):
#         start, start_jump, end, end_jump = parse_cdna_jumps(mutation)
#         result['start'] = start
#         result['end'] = end
#         result['start_jump'] = start_jump
#         result['end_jump'] = end_jump
#
#     return result

def mutation_dict(mut_row):
    # Skip if DNAChange is missing or empty
    if pd.isna(mut_row.get('DNAChange')) or mut_row['DNAChange'].strip() == '':
        return None

    try:
        # Split into ID and mutation string
        transcript_id, mutation = mut_row['DNAChange'].split(':')
    except ValueError:
        return None  # if it can't unpack into two parts

    result = {'id': transcript_id, 'start_jump': 0, 'end_jump': 0}  # default jumps = 0

    if mut_row['VariantType'] == 'SNV':
        match = re.search(r'\d+', mutation)
        if match:
            start, end = match.span()
            location = mutation[start:end]
            from_nt, to_nt = mutation[end:].split('>')
            result.update({'start': int(location), 'end': int(location), 'ref': from_nt, 'alt': to_nt, 'type': mut_row['VariantType']})
        else:
            result.update({'start': None, 'end': None, 'ref': None, 'alt': None, 'type': None})

    elif mut_row['VariantType'] == 'deletion':
        numbers = re.findall(r'\d+', mutation)
        if len(numbers) == 1:
            result.update({'start': int(numbers[0]), 'end': int(numbers[0]), 'type': mut_row['VariantType']})
        elif len(numbers) >= 2:
            result.update({'start': int(numbers[0]), 'end': int(numbers[1]), 'type': mut_row['VariantType']})
        else:
            result.update({'start': None, 'end': None, 'type': None})

    elif mut_row['VariantType'] == 'substitution':
        numbers = re.findall(r'\d+', mutation)
        if len(numbers) == 1:
            result.update({'start': int(numbers[0]), 'end': int(numbers[0])})
        elif len(numbers) >= 2:
            result.update({'start': int(numbers[0]), 'end': int(numbers[1])})
        else:
            result.update({'start': None, 'end': None, 'type': None})

        if 'delins' in mutation:
            result['type'] = 'delins'
            result['subs'] = mutation.split('delins')[-1]
        elif 'inv' in mutation:
            result['type'] = 'inversion'

    elif mut_row['VariantType'] == 'insertion':
        numbers = re.findall(r'\d+', mutation)
        if len(numbers) == 1:
            result.update({'start': int(numbers[0]), 'end': int(numbers[0])})
        elif len(numbers) >= 2:
            result.update({'start': int(numbers[0]), 'end': int(numbers[1])})
        else:
            result.update({'start': None, 'end': None, 'type': None})

        if 'ins' in mutation:
            result['type'] = 'insertion'
            result['ins'] = mutation.split('ins')[-1]
        elif 'dup' in mutation:
            result['type'] = 'duplication'
    else:
        return None

    result['info'] = mut_row.get('VariantInfo')

    # --- Handle splice donor/acceptor ---
    if isinstance(result['info'], str) and (
        "splice_donor_variant" in result['info'] or "splice_acceptor_variant" in result['info']
    ):
        start, start_jump, end, end_jump = parse_cdna_jumps(mutation)
        result['start'] = start
        result['end'] = end
        result['start_jump'] = start_jump
        result['end_jump'] = end_jump

    return result

def comp_strand(seq):
    seq = seq[::-1]
    new_seq = ''
    for n in range(0,len(seq)):
        if seq[n] == 'A' or seq[n] =='a':
            new_seq += 'T'
        elif seq[n] == 'T' or seq[n] =='t' :
            new_seq += 'A'
        elif seq[n] == 'C' or seq[n] == 'c':
            new_seq += 'G'
        elif seq[n] == 'G' or seq[n] == 'g':
            new_seq += 'C'
        else:
            new_seq += seq[n]
    return new_seq

# ================================================= OLD MUTATE FUNCTION ================================================
# def mutate(mut_dict, shift, sequence):
#     # Prepare numeric values once
#     start = int(mut_dict.get('start', 0)) + int(mut_dict.get('start_jump', 0))
#     end = int(mut_dict.get('end', 0)) + int(mut_dict.get('end_jump', 0))
#     shift = int(shift)  # ensure int
#
#     mut_seq = list(sequence)
#
#     if mut_dict['type'] == 'SNV':
#
#         loc = int(mut_dict['start']) + shift - 1
#         ref = mut_dict['ref']
#         alt = mut_dict['alt']
#
#         mut_seq[loc] = alt
#
#         # if sequence[loc] == ref:
#         #     mut_seq[loc] = alt
#         #
#         # else:
#         #     mut_seq = 'mismatch'
#         #     print(f'mismatch at {loc}, expected {ref}, got {mut_seq[loc]}')
#
#     elif mut_dict['type'] == 'deletion':
#
#         start = int(mut_dict['start']) + shift - 1
#         end = int(mut_dict['end']) + shift
#         del mut_seq[start:end]
#
#     elif mut_dict['type'] == 'insertion':
#
#         start = int(mut_dict['start']) + shift - 1
#         mut_seq[start] += mut_dict['ins']
#
#     elif mut_dict['type'] == 'duplication':
#
#         start = int(mut_dict['start']) + shift - 1
#         end = int(mut_dict['end']) + shift
#         dup = mut_seq[start:end]
#         mut_seq[end - 1] += ''.join(dup)
#
#     elif mut_dict['type'] == 'delins':
#
#         start = int(mut_dict['start']) + shift - 1
#         end = int(mut_dict['end']) + shift
#         del mut_seq[start:end]
#         mut_seq[start - 1] += mut_dict['subs']
#
#     elif mut_dict['type'] == 'inversion':
#
#         start = int(mut_dict['start']) + shift - 1
#         end = int(mut_dict['end']) + shift
#         sub = ''.join(mut_seq[start:end])
#         new = comp_strand(sub)
#         del mut_seq[start:end]
#         mut_seq[start - 1] += new
#
#
#     return ''.join(mut_seq)

def mutate(mut_dict, shift, sequence):
    """
    Apply a text-based sequence modification from a mutation dictionary.
    This is purely string manipulation for data analysis.
    """

    # Prepare numeric values once
    start = int(mut_dict.get('start', 0)) + int(mut_dict.get('start_jump', 0))
    end = int(mut_dict.get('end', 0)) + int(mut_dict.get('end_jump', 0))
    shift = int(shift)  # ensure int

    mut_seq = list(sequence)  # convert to mutable list of chars

    if mut_dict['type'] == 'SNV':
        loc = start + shift - 1
        mut_seq[loc] = mut_dict['alt']

    elif mut_dict['type'] == 'deletion':
        start_i = start + shift - 1
        end_i = end + shift
        del mut_seq[start_i:end_i]

    elif mut_dict['type'] == 'insertion':
        start_i = start + shift - 1
        mut_seq[start_i] += mut_dict.get('ins', '')

    elif mut_dict['type'] == 'duplication':
        start_i = start + shift - 1
        end_i = end + shift
        dup = mut_seq[start_i:end_i]
        mut_seq[end_i - 1] += ''.join(dup)

    elif mut_dict['type'] == 'delins':
        start_i = start + shift - 1
        end_i = end + shift
        del mut_seq[start_i:end_i]
        mut_seq[start_i - 1] += mut_dict.get('subs', '')

    elif mut_dict['type'] == 'inversion':
        start_i = start + shift - 1
        end_i = end + shift
        sub = ''.join(mut_seq[start_i:end_i])
        new = comp_strand(sub)  # assuming you have this for text complement
        del mut_seq[start_i:end_i]
        mut_seq[start_i - 1] += new

    return ''.join(mut_seq)


def final_func(seq_fasta, expression_file, mutation_file, cds_gtf_file, cell_line, top_n):
    # read the mutation and expression data
    mut_data = get_mutations_of_cell_line(cell_line, mutation_file)
    print('read mutation data')
    exp_data = get_expression_of_cell_line(cell_line, expression_file)
    exp_data = exp_data.head(top_n)
    print('read expression data')

    # add new columns
    exp_data['Transcript_ID'] = None
    exp_data['Mutated Transcript Sequence'] = None
    exp_data['Original Transcript Sequence'] = None

    records = list(SeqIO.parse(seq_fasta, "fasta"))
    record_by_hugo = {}
    for record in records:
        match = re.search(r'gene_symbol:([^\s]+)', record.description)
        if match:
            hugo = match.group(1)
            if hugo not in record_by_hugo:  # only take the first occurrence
                record_by_hugo[hugo] = record
    # add to exp_data the id
    for idx, row in exp_data.iterrows():
        hugo_symbol = row['Gene'].split()[0]

        record = record_by_hugo.get(hugo_symbol)
        if record:
            exp_data.at[idx, 'Transcript_ID'] = record.id
            exp_data.at[idx, 'Original Transcript Sequence'] = str(record.seq)
        else:
            exp_data.at[idx, 'Transcript_ID'] = None
            exp_data.at[idx, 'Original Transcript Sequence'] = None
        print('got sequence')
    print('finished getting sequences')

    ############################################################################################################

    # Parse GTF file
    cds_data = prepare_gtf_dataframe(cds_gtf_file)
    print('parsed GTF')

    # Convert exp_data to dict indexed by transcript ID for O(1) lookup
    exp_data_indexed = exp_data.set_index('Transcript_ID')

    for idx, row in mut_data.iterrows():
        mut_dict = mutation_dict(row)

        if mut_dict is not None:
            tid = mut_dict['id']

            if tid in exp_data_indexed.index:
                try:
                    shift = get_distance_to_start_codon_from_df(cds_data, tid)
                    print(f"{tid} shift: {shift}")
                    if pd.isna(exp_data_indexed.at[tid, 'Mutated Transcript Sequence']):
                        seq = exp_data_indexed.at[tid, 'Original Transcript Sequence']
                    else:
                        seq = exp_data_indexed.at[tid, 'Mutated Transcript Sequence']

                    mutated_seq = mutate(mut_dict, shift, seq)
                    exp_data_indexed.at[tid, 'Mutated Transcript Sequence'] = mutated_seq
                    print(f'mutated {tid}')
                except Exception as e:
                    print(f"Skipping {tid} due to error: {e}")

    ##############################################################################################################

    exp_data = exp_data_indexed.reset_index()
    exp_data.to_csv(f'{cell_line}_transcriptome_top{top_n}.csv', index=False)
    print('finished!')
    return exp_data


# celline_list = ['ACH-001328',
#                 'ACH-000463',
#                 'ACH-001188',
#                 'ACH-001086',
#                 'ACH-000739',
#                 'ACH-000232',
#                 ]

celline_list = ['ACH-000681']

# =================================== Create the mutated cell line transcriptome files =================================

# for line in celline_list:
#     final_func(
#         'Homo_sapiens.GRCh38.cdna.all.fa',
#         'OmicsExpressionProteinCodingGenesTPMLogp1.csv',
#         'OmicsSomaticMutations.csv',
#         'gencode.v48.chr_patch_hapl_scaff.annotation.gtf',
#         line, 500
#     )

A549 = "ACH-000681"
path_1 = "/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/"
exp_file = path_1 + "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
mut_file = path_1 + "OmicsSomaticMutations.csv"
seq_file = path_1 + "Homo_sapiens.GRCh38.cdna.all.fa"
gtf_file = path_1 + "gencode.v48.chr_patch_hapl_scaff.annotation.gtf"


# final_func(seq_file, exp_file, mut_file, gtf_file, A549, 30000)

def mutate_transcriptome(expression, mutations, annotation_dict):
    print('Mutating transcriptomes...')
    exp_data_indexed = expression.set_index('Transcript_ID')
    exp_data_indexed["Mutated Transcript Sequence"] = ""
    for idx, row in mutations.iterrows():
        mut_dict = mutation_dict(row)

        if mut_dict is not None:
            tid = mut_dict['id']
            mut_idx = mut_dict['start']

            if tid in exp_data_indexed.index:
                try:
                    shift = find_shift(annotation_dict[tid], mut_idx)
                    print(f"{tid} shift: {shift}")
                    if pd.isna(exp_data_indexed.at[tid, 'Mutated Transcript Sequence']):
                        seq = exp_data_indexed.at[tid, 'Original Transcript Sequence']
                    else:
                        seq = exp_data_indexed.at[tid, 'Mutated Transcript Sequence']

                    mutated_seq = mutate(mut_dict, shift, seq)
                    exp_data_indexed.at[tid, 'Mutated Transcript Sequence'] = mutated_seq
                    if seq != mutated_seq:
                        print(f'really mutated {tid}')
                    else:
                        print(f'didnt really do anything... check {tid}')
                except Exception as e:
                    print(f"Skipping {tid} due to error: {e}")
    return exp_data_indexed