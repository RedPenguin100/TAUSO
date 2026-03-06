import platform
import random
import os
import re
import subprocess
import tempfile
import uuid

from pathlib import Path
from typing import Dict, List

from Bio.Seq import Seq

from .Interaction import Interaction
from ...common.consts import OUT_FOLDER

if platform.system() == "Linux" and os.path.exists("/dev/shm"):
    TMP_PATH = Path("/dev/shm/tauso_risearch_tmp")
else:
    TMP_PATH = Path(tempfile.gettempdir()) / "tauso_risearch_tmp"



def dump_target_file(target_filename: str, name_to_sequence: Dict[str, str]):
    tmp_path = TMP_PATH / target_filename
    TMP_PATH.mkdir(exist_ok=True)
    with open(tmp_path, "w+") as f:
        for name, sequence in name_to_sequence.items():
            f.write(f">{name}\n{sequence}\n")
    return tmp_path


def get_risearch_path():
    binary_path = str(OUT_FOLDER / 'risearch_executable')

    if not os.path.exists(binary_path):
        raise FileNotFoundError(f"Binary missing at {binary_path}")
    return binary_path


def get_trigger_mfe_scores_by_risearch(trigger: str, name_to_sequence: Dict[str, str],
                                       interaction_type: Interaction = Interaction.RNA_DNA_NO_WOBBLE,
                                       minimum_score: int = 900, neighborhood: int = 0, parsing_type=None,
                                       target_file_cache=None, transpose=False, unique_id=None) -> str:
    if not name_to_sequence:
        raise ValueError("name_to_sequence is empty!")
    TMP_PATH.mkdir(parents=True, exist_ok=True)

    risearch_path = get_risearch_path()

    if unique_id is None:
        unique_id = uuid.uuid4().hex

    if target_file_cache is None:
        target_filename = f"target-{unique_id}.fa"
        target_path = Path(dump_target_file(target_filename, name_to_sequence)).resolve()
    else:
        target_path = Path(target_file_cache).resolve()

    query_path = (TMP_PATH / f"query-{unique_id}.fa").resolve()

    # Use 'w' and flush to ensure data is written before subprocess starts
    with open(query_path, "w") as f:
        f.write(f">trigger\n{Seq(trigger).reverse_complement_rna()}\n")
        f.flush()
        os.fsync(f.fileno())

    # --- MINIMAL ADDITION: Validation Exceptions ---
    for p, name in [(target_path, "Target"), (query_path, "Query")]:
        if not p.exists():
            raise FileNotFoundError(f"{name} file was not created at {p}")
        if p.stat().st_size == 0:
            raise ValueError(f"{name} file is empty (0 bytes) at {p}. Disk might be full or write failed.")

    if interaction_type == Interaction.RNA_DNA_NO_WOBBLE:
        m = 'su95_noGU'
    elif interaction_type == Interaction.RNA_RNA:
        m = 't04'
    else:
        raise ValueError(f"Unsupported interaction type: {interaction_type}={str({interaction_type})}")

    args = [str(risearch_path), "-q", str(query_path), "-t", str(target_path), "-s", f"{minimum_score}",
            "-d", "30", "-m", m, '-n', f"{neighborhood}"]
    if transpose:
        args.append("-R")

    try:
        if parsing_type is not None:
            args.append(f'-p{parsing_type}')
        # Capture stderr to see the actual RIsearch error if it returns 255
        result = subprocess.check_output(args, universal_newlines=True, text=True, cwd=str(OUT_FOLDER), stderr=subprocess.STDOUT)
    finally:
        if target_file_cache is None and target_path.exists():
            os.remove(target_path)
        if query_path.exists():
            query_path.unlink()
    return result

'''
note: I used hashing here to be able to run multiple sequences (triggers in my case) in parallel (for multiple triggers) without conflicts with the files,
you can remove it if it's not an issue for you 
note: in my case I first preformed reverse_complement_rna to the trigger because my goal was to find sequences that might undesirably bind to the
trigger binding site of the toehold, adjust it for your particular case
note: play with the d and s parameters I passed here, and with other parameters that might be relevant for your usecase
'''


def _parse_mfe_scores_2(result):
    if not result:
        return [[]]

    lines = result.split('\n')
    target_to_energies = dict()
    for line in lines:
        if line == '':
            continue
        line_parts = line.split('\t')
        target_name = line_parts[3]
        target_energy = line_parts[-1]
        if line_parts[3] in target_to_energies:
            target_to_energies[target_name].append(float(target_energy))
        else:
            target_to_energies[target_name] = [float(target_energy)]

    # ignore the keys at this point
    return list(target_to_energies.values())


# in my case I was only interested in the energy scores and didn't care about the actual sequences, this is the parsing I used in case it helps you
def get_mfe_scores(result: str, parsing_type=None) -> List[List[float]]:
    mfe_results = []

    if parsing_type is None:
        for gene_result in result.split("\n\nquery trigger")[1:]:
            stripped_result = gene_result.strip()
            regex_results = re.findall("Free energy \[kcal/mol\]: [0-9-.]+ ", stripped_result)
            mfe_results.append(
                [float(regex_result.replace('Free energy [kcal/mol]: ', '').strip()) for regex_result in regex_results])
    elif parsing_type == '2':
        return _parse_mfe_scores_2(result)

    return mfe_results
