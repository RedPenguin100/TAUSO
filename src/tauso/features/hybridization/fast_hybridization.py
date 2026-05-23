import os
import platform
import re
import subprocess
import tempfile
import uuid
from importlib.resources import files
from pathlib import Path
from typing import Dict, List, Tuple

from ...util import get_antisense_rna
from .Interaction import Interaction

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


def get_risearch_path() -> str:
    binary_path = files("tauso") / "out" / "risearch_executable"
    if not binary_path.is_file():
        raise FileNotFoundError(f"Binary missing at {binary_path}")
    return str(binary_path)


def _interaction_mode(interaction_type: Interaction) -> str:
    if interaction_type == Interaction.RNA_DNA_NO_WOBBLE:
        return "su95_noGU"
    if interaction_type == Interaction.RNA_RNA:
        return "t04"
    raise ValueError(f"Unsupported interaction type: {interaction_type}")


def _build_risearch_args(
    query_path: Path,
    target_path: Path,
    min_score: int,
    mode: str,
    neighborhood: int,
    transpose: bool,
    parsing_type,
) -> List[str]:
    args = [
        get_risearch_path(),
        "-q",
        str(query_path),
        "-t",
        str(target_path),
        "-s",
        str(min_score),
        "-d",
        "30",
        "-m",
        mode,
        "-n",
        str(neighborhood),
    ]
    if transpose:
        args.append("-R")
    if parsing_type is not None:
        args.append(f"-p{parsing_type}")
    return args


def _run_risearch(args: List[str]) -> str:
    TMP_PATH.mkdir(parents=True, exist_ok=True)
    return subprocess.check_output(
        args,
        universal_newlines=True,
        text=True,
        cwd=str(TMP_PATH),
        stderr=subprocess.STDOUT,
        stdin=subprocess.DEVNULL,
    )


def get_trigger_mfe_scores_by_risearch(
    trigger: str,
    name_to_sequence: Dict[str, str],
    interaction_type: Interaction = Interaction.RNA_DNA_NO_WOBBLE,
    minimum_score: int = 900,
    neighborhood: int = 0,
    parsing_type=None,
    target_file_cache=None,
    transpose=False,
    unique_id=None,
) -> str:
    if not name_to_sequence:
        raise ValueError("name_to_sequence is empty!")
    TMP_PATH.mkdir(parents=True, exist_ok=True)

    if unique_id is None:
        unique_id = uuid.uuid4().hex

    if target_file_cache is None:
        target_path = Path(dump_target_file(f"target-{unique_id}.fa", name_to_sequence)).resolve()
    else:
        target_path = Path(target_file_cache).resolve()

    query_path = (TMP_PATH / f"query-{unique_id}.fa").resolve()
    with open(query_path, "w") as f:
        f.write(f">trigger\n{get_antisense_rna(trigger)}\n")

    for p, name in [(target_path, "Target"), (query_path, "Query")]:
        if not p.exists():
            raise FileNotFoundError(f"{name} file was not created at {p}")
        if p.stat().st_size == 0:
            raise ValueError(f"{name} file is empty at {p}. Disk might be full.")

    mode = _interaction_mode(interaction_type)
    args = _build_risearch_args(query_path, target_path, minimum_score, mode, neighborhood, transpose, parsing_type)

    try:
        return _run_risearch(args)
    finally:
        if target_file_cache is None and target_path.exists():
            os.remove(target_path)
        if query_path.exists():
            query_path.unlink()


def get_triggers_mfe_scores_batch(
    trigger_id_seq_pairs: List[Tuple[str, str]],
    target_file_path,
    interaction_type: Interaction = Interaction.RNA_DNA_NO_WOBBLE,
    minimum_score: int = 900,
    neighborhood: int = 0,
    parsing_type=None,
    transpose=False,
    batch_id=None,
) -> str:
    """Run RIsearch once for all (id, trigger) pairs against a pre-built target file."""
    if not trigger_id_seq_pairs:
        return ""

    TMP_PATH.mkdir(parents=True, exist_ok=True)

    if batch_id is None:
        batch_id = uuid.uuid4().hex

    query_path = (TMP_PATH / f"query-batch-{batch_id}.fa").resolve()
    lines: List[str] = []
    for query_id, trigger in trigger_id_seq_pairs:
        lines.append(f">{query_id}")
        lines.append(get_antisense_rna(trigger))
    query_path.write_text("\n".join(lines) + "\n")

    mode = _interaction_mode(interaction_type)
    args = _build_risearch_args(
        query_path, target_file_path, minimum_score, mode, neighborhood, transpose, parsing_type
    )

    try:
        return _run_risearch(args)
    finally:
        if query_path.exists():
            query_path.unlink()


def stream_triggers_mfe_hits(
    trigger_id_seq_pairs: List[Tuple[str, str]],
    target_file_path,
    interaction_type: Interaction = Interaction.RNA_DNA_NO_WOBBLE,
    minimum_score: int = 900,
    neighborhood: int = 0,
    parsing_type="2",
    transpose=False,
    batch_id=None,
    parse_chunk_rows: int = 100_000,
    usecols=("trigger", "energy"),
):
    """Streaming counterpart to get_triggers_mfe_scores_batch.

    Yields pandas DataFrames (chunks of `parse_chunk_rows` rows each) with
    columns ['trigger', 'energy'] parsed directly from the RIsearch subprocess
    stdout via pd.read_csv(chunksize=...). This combines:

      - Bounded memory (one chunk in flight, vs the full TSV materialized
        by get_triggers_mfe_scores_batch + parse_risearch_output)
      - C-level parsing speed (pd.read_csv is much faster than pure-Python
        line iteration on multi-million-line outputs)

    Assumes parsing_type="2" — the standard 8-column TSV emitted by the
    scoring callers: trigger, t_start, t_end, target, ta_start, ta_end, score, energy.
    Pass `usecols` to load only the columns the caller needs (defaults to
    trigger+energy; the off_target_specific path needs trigger+target+energy).
    """
    if not trigger_id_seq_pairs:
        return

    TMP_PATH.mkdir(parents=True, exist_ok=True)
    if batch_id is None:
        batch_id = uuid.uuid4().hex

    query_path = (TMP_PATH / f"query-stream-{batch_id}.fa").resolve()
    lines: List[str] = []
    for query_id, trigger in trigger_id_seq_pairs:
        lines.append(f">{query_id}")
        lines.append(get_antisense_rna(trigger))
    query_path.write_text("\n".join(lines) + "\n")

    mode = _interaction_mode(interaction_type)
    args = _build_risearch_args(
        query_path, target_file_path, minimum_score, mode, neighborhood, transpose, parsing_type
    )

    columns = [
        "trigger", "trigger_start", "trigger_end",
        "target", "target_start", "target_end",
        "score", "energy",
    ]
    usecols_list = list(usecols)
    dtype_map = {"trigger": "string", "target": "string", "energy": "float64"}
    dtypes = {k: v for k, v in dtype_map.items() if k in usecols_list}

    try:
        import pandas as pd

        with subprocess.Popen(
            args,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            stdin=subprocess.DEVNULL,
            cwd=str(TMP_PATH),
            bufsize=-1,
        ) as proc:
            try:
                yield from pd.read_csv(
                    proc.stdout,
                    sep="\t",
                    header=None,
                    names=columns,
                    usecols=usecols_list,
                    dtype=dtypes,
                    chunksize=parse_chunk_rows,
                )
            except pd.errors.EmptyDataError:
                # No hits at all — RIsearch produced empty stdout.
                pass
            finally:
                # Drain any remaining bytes pandas didn't read so RIsearch can exit
                if proc.stdout is not None and not proc.stdout.closed:
                    try:
                        proc.stdout.read()
                    except Exception:
                        pass
            rc = proc.wait()
            if rc != 0:
                raise subprocess.CalledProcessError(rc, args)
    finally:
        if query_path.exists():
            query_path.unlink()


def _parse_mfe_scores_2(result):
    if not result:
        return [[]]

    lines = result.split("\n")
    target_to_energies = dict()
    for line in lines:
        if line == "":
            continue
        line_parts = line.split("\t")
        target_name = line_parts[3]
        target_energy = line_parts[-1]
        if line_parts[3] in target_to_energies:
            target_to_energies[target_name].append(float(target_energy))
        else:
            target_to_energies[target_name] = [float(target_energy)]

    return list(target_to_energies.values())


def get_mfe_scores(result: str, parsing_type=None) -> List[List[float]]:
    mfe_results = []

    if parsing_type is None:
        for gene_result in result.split("\n\nquery trigger")[1:]:
            stripped_result = gene_result.strip()
            regex_results = re.findall("Free energy \\[kcal/mol\\]: [0-9-.]+ ", stripped_result)
            mfe_results.append([float(r.replace("Free energy [kcal/mol]: ", "").strip()) for r in regex_results])
    elif parsing_type == "2":
        return _parse_mfe_scores_2(result)

    return mfe_results
