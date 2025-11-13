from .consts import YEAST_GFF_PATH, YEAST_FASTA_PATH, YEAST_THREE_PRIME_UTR, YEAST_FIVE_PRIME_UTR, \
    YEAST_README, YEAST_GFF_DB_PATH


def validate_yeast_files():
    if not YEAST_GFF_PATH.is_file() and not YEAST_GFF_DB_PATH.is_file():
        raise FileNotFoundError(
            f'Neither {YEAST_GFF_PATH.name} nor {YEAST_GFF_DB_PATH.name} are present, please read the {YEAST_README} file')
    if not YEAST_FASTA_PATH.is_file():
        raise FileNotFoundError(f'{YEAST_FASTA_PATH} is not present, please read the {YEAST_README} file')

    if not YEAST_THREE_PRIME_UTR.is_file():
        raise FileNotFoundError(
            f'{YEAST_THREE_PRIME_UTR} is not present, please read the {YEAST_README} file'
        )
    if not YEAST_FIVE_PRIME_UTR.is_file():
        raise FileNotFoundError(f'{YEAST_FIVE_PRIME_UTR} is not present, please read the {YEAST_README} file')
