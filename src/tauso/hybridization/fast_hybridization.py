import random
import os
import subprocess

from enum import Enum
from typing import Dict

from Bio.Seq import Seq

from ..consts import RISEARCH1_BINARY_PATH, TMP_PATH, RISEARCH1_PATH


class Interaction(Enum):
    RNA_RNA = "RNA_RNA"
    DNA_RNA_NO_WOBBLE = "DNA_RNA_NO_WOBBLE"
    DNA_DNA = "DNA_DNA"
    MODIFIED = "MODIFIED"


def dump_target_file(target_filename: str, name_to_sequence: Dict[str, str]):
    tmp_path = TMP_PATH / target_filename
    TMP_PATH.mkdir(exist_ok=True)
    with open(tmp_path, "w") as f:
        for name, sequence in name_to_sequence.items():
            f.write(f">{name}\n{sequence}\n")
    return tmp_path


def get_trigger_mfe_scores_by_risearch(trigger: str, name_to_sequence: Dict[str, str],
                                       interaction_type: Interaction = Interaction.DNA_RNA_NO_WOBBLE,
                                       minimum_score: int = 900, neighborhood: int = 0, parsing_type=None,
                                       target_file_cache=None, e_modified=(249, 100)) -> str:
    if not RISEARCH1_BINARY_PATH.is_file():
        raise FileNotFoundError(
            f'RIsearch binary is not found at {RISEARCH1_BINARY_PATH}. Please read the "RIsearch" part of the project readme.')
    # used to dump cached files
    TMP_PATH.mkdir(exist_ok=True)

    hash = random.getrandbits(64)

    if target_file_cache is None:
        target_filename = f"target-{hash}.fa"
        target_path = dump_target_file(target_filename, name_to_sequence)
    else:
        target_filename = target_file_cache
        target_path = target_filename

    query_filename = f"query-{hash}.fa"
    query_path = TMP_PATH / query_filename

    with open(query_path, "w") as f:
        f.write(f">trigger\n{Seq(trigger).reverse_complement_rna()}\n")

    if interaction_type == Interaction.DNA_RNA_NO_WOBBLE:
        m = 'su95_noGU'
    elif interaction_type == Interaction.RNA_RNA:
        m = 't04'
    elif interaction_type == Interaction.MODIFIED:
        e1, e2 = e_modified
        m = f'modified {e1} {e2}'
    else:
        raise ValueError(f"Unsupported interaction type: {interaction_type}")

    try:
        args = [RISEARCH1_BINARY_PATH, "-q", str(query_path), "-t", str(target_path), "-s", f"{minimum_score}",
                "-d", "30", "-m", m, '-n', f"{neighborhood}"]
        if parsing_type is not None:
            args.append(f'-p{parsing_type}')

        result = subprocess.check_output(args, universal_newlines=True, text=True, cwd=str(RISEARCH1_PATH)
                                         )
    finally:
        if target_file_cache is None:
            os.remove(target_path)
        query_path.unlink()
    # logger.info(f"finished evaluating a single trigger: {trigger}")
    return result
