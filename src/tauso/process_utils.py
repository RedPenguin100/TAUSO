from concurrent.futures import ProcessPoolExecutor, as_completed

from tqdm import tqdm

from .timer import Timer


def get_simplified_fasta_dict(fasta_dict):
    simplified_fasta_dict = dict()
    for locus_tag, locus_info in fasta_dict.items():
        simplified_fasta_dict[locus_tag] = str(locus_info.upper().seq)
    return simplified_fasta_dict


def validated_get_simplified_fasta_dict(fasta_dict, simplified_fasta_dict):
    if simplified_fasta_dict is None and fasta_dict is None:
        raise ValueError('Either simplified_fasta_dict or fasta_dict must be specified')

    if simplified_fasta_dict is None:
        return get_simplified_fasta_dict(fasta_dict)
    return simplified_fasta_dict



def validate_organism(organism: str):
    organisms = ['human', 'yeast']
    if organism not in organisms:
        raise ValueError(f'Organism={organism} must be in {organisms}')

