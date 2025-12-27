def select_best_transcript(full_gene_str, gene_map, gtf_annotations):
    """
    Args:
        full_gene_str: "TSPAN6 (7105)"
        gene_map: Dictionary built from GTF mapping gene_names to transcript_ids
    """
    # 1. Get the clean symbol: "TSPAN6"
    gene_symbol = full_gene_str.split(' ')[0].strip()

    # 2. Try direct match
    possible_tids = gene_map.get(gene_symbol)

    # 3. Fallback: If not found, try a case-insensitive match
    if not possible_tids:
        # Some symbols in files are lowercase, GTF is usually uppercase
        possible_tids = gene_map.get(gene_symbol.upper())

    if not possible_tids:
        return None

    # Find the longest transcript (pre-mRNA span)
    best_tid = None
    max_len = -1
    for tid in possible_tids:
        data = gtf_annotations[tid]
        length = data['end'] - data['start']
        if length > max_len:
            max_len = length
            best_tid = tid

    return best_tid


def build_gene_to_transcript_index(gtf_annotations):
    """
    Creates a helper map: Gene Symbol -> List of Transcript IDs.
    Returns: {'MT-ATP8': ['ENST000...', 'ENST000...'], ...}
    """
    gene_map = {}
    for tid, data in gtf_annotations.items():
        gname = data.get('gene_name')
        if gname:
            if gname not in gene_map:
                gene_map[gname] = []
            gene_map[gname].append(tid)
    return gene_map


def get_reverse_complement(seq):
    """Reverses and complements a DNA sequence."""
    complement_map = str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn")
    return seq.translate(complement_map)[::-1]