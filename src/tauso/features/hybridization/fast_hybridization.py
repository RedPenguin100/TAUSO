import os
import platform
import re
import subprocess
import tempfile
import uuid
from contextlib import contextmanager
from pathlib import Path
from typing import Callable, Dict, List, NamedTuple, Tuple

import risearch_tauso

from ...util import get_antisense_rna
from .interaction import Interaction

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
    return risearch_tauso.executable_path()


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


# The fixed 8-column parsing_type="2" TSV that RIsearch emits.
RISEARCH_COLUMNS = (
    "trigger",
    "trigger_start",
    "trigger_end",
    "target",
    "target_start",
    "target_end",
    "score",
    "energy",
)


@contextmanager
def _risearch_stdout(
    trigger_id_seq_pairs: List[Tuple[str, str]],
    target_file_path,
    *,
    interaction_type: Interaction = Interaction.RNA_DNA_NO_WOBBLE,
    minimum_score: int = 900,
    neighborhood: int = 0,
    parsing_type="2",
    transpose=False,
    batch_id=None,
):
    """The single RIsearch invocation. Writes the batched query FASTA, builds the args,
    spawns the subprocess, and yields its live stdout for the caller to consume. On exit
    it drains stdout, checks the exit code, and removes the query file. Every RIsearch
    consumer (streaming aggregation, raw-hits collection) goes through here so the
    process/FASTA/cleanup logic lives in exactly one place.
    """
    TMP_PATH.mkdir(parents=True, exist_ok=True)
    if batch_id is None:
        batch_id = uuid.uuid4().hex

    query_path = (TMP_PATH / f"query-{batch_id}.fa").resolve()
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
        with subprocess.Popen(
            args,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            stdin=subprocess.DEVNULL,
            cwd=str(TMP_PATH),
            bufsize=-1,
        ) as proc:
            try:
                yield proc.stdout
            finally:
                # Drain any unread stdout so the child can exit cleanly.
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


def _pa_column_types():
    """pyarrow column types for the RISEARCH_COLUMNS TSV (built lazily — pyarrow is a
    heavy import kept out of module import time)."""
    import pyarrow as pa

    return {
        "trigger": pa.string(),
        "trigger_start": pa.int64(),
        "trigger_end": pa.int64(),
        "target": pa.string(),
        "target_start": pa.int64(),
        "target_end": pa.int64(),
        "score": pa.int64(),
        "energy": pa.float64(),
    }


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


class RisearchAggregation(NamedTuple):
    """How to reduce streamed RIsearch hits into a {key: value} dict.

    columns  — which of the 8 RIsearch TSV columns to parse (a subset of
               trigger/target/energy).
    combine  — (pyarrow.RecordBatch) -> partial pyarrow.Table, applied per parsed
               block so millions of hits never materialize at once.
    finalize — (concatenated partial Table) -> the result dict.

    See min_energy_by_pair_multi_cutoff and sum_exp_by_trigger_multi_cutoff.
    """

    columns: Tuple[str, ...]
    combine: Callable
    finalize: Callable


def parse_risearch_hits_pyarrow(
    trigger_id_seq_pairs: List[Tuple[str, str]],
    target_file_path,
    *,
    aggregation: RisearchAggregation,
    interaction_type: Interaction = Interaction.RNA_DNA_NO_WOBBLE,
    minimum_score: int = 900,
    neighborhood: int = 0,
    parsing_type="2",
    transpose=False,
    batch_id=None,
    block_size: int = 64 << 20,
):
    """Run RIsearch and stream its stdout through pyarrow, block by block.

    Pure mechanism — no energy biology lives here. It owns the subprocess, the
    pyarrow CSV reader, the per-block loop, the concatenation of partials and the
    temp-file cleanup. pyarrow's CSV reader runs in C++ and releases the GIL, so
    callers scale under a ThreadPoolExecutor (the old pandas parse held the GIL
    and plateaued at ~1.7x). Memory is bounded by `block_size` plus the small
    partials `aggregation` returns.

    What to compute is supplied by `aggregation` (a RisearchAggregation): it
    declares which columns to parse, reduces each parsed block to a partial table
    (`combine`), and reduces the concatenated partials to the result dict
    (`finalize`). The RIsearch output is always the 8-column parsing_type="2" TSV:
    trigger, t_start, t_end, target, ta_start, ta_end, score, energy.
    """
    if not trigger_id_seq_pairs:
        return {}

    import pyarrow as pa
    import pyarrow.csv as pacsv

    type_map = _pa_column_types()
    read_opts = pacsv.ReadOptions(column_names=list(RISEARCH_COLUMNS), use_threads=False, block_size=block_size)
    parse_opts = pacsv.ParseOptions(delimiter="\t")
    convert_opts = pacsv.ConvertOptions(
        include_columns=list(aggregation.columns),
        column_types={col: type_map[col] for col in aggregation.columns},
    )

    parts = []
    with _risearch_stdout(
        trigger_id_seq_pairs,
        target_file_path,
        interaction_type=interaction_type,
        minimum_score=minimum_score,
        neighborhood=neighborhood,
        parsing_type=parsing_type,
        transpose=transpose,
        batch_id=batch_id,
    ) as stdout:
        try:
            reader = pacsv.open_csv(
                stdout, read_options=read_opts, parse_options=parse_opts, convert_options=convert_opts
            )
            for batch in reader:
                if batch.num_rows == 0:
                    continue
                parts.append(aggregation.combine(batch))
        except pa.ArrowInvalid:
            # No hits — RIsearch produced empty stdout.
            pass

    if not parts:
        return {}

    return aggregation.finalize(pa.concat_tables(parts))


def risearch_hits_dataframe(
    trigger_id_seq_pairs: List[Tuple[str, str]],
    target_file_path,
    *,
    interaction_type: Interaction = Interaction.RNA_DNA_NO_WOBBLE,
    minimum_score: int = 900,
    neighborhood: int = 0,
    parsing_type="2",
    transpose=False,
    batch_id=None,
    block_size: int = 64 << 20,
):
    """Run one batched RIsearch and return ALL hits as a DataFrame (the 8 RISEARCH_COLUMNS).

    Same single batched invocation as parse_risearch_hits_pyarrow and streamed through
    pyarrow, but materialises the full table instead of reducing it — intended for tests
    and debugging, not the hot path (production uses parse_risearch_hits_pyarrow with an
    aggregation so memory stays bounded). Replaces the old get_triggers_mfe_scores_batch
    + parse_risearch_output pair.
    """
    import pandas as pd
    import pyarrow as pa
    import pyarrow.csv as pacsv

    empty = pd.DataFrame(columns=list(RISEARCH_COLUMNS))
    if not trigger_id_seq_pairs:
        return empty

    read_opts = pacsv.ReadOptions(column_names=list(RISEARCH_COLUMNS), use_threads=False, block_size=block_size)
    parse_opts = pacsv.ParseOptions(delimiter="\t")
    convert_opts = pacsv.ConvertOptions(column_types=_pa_column_types())

    tables = []
    with _risearch_stdout(
        trigger_id_seq_pairs,
        target_file_path,
        interaction_type=interaction_type,
        minimum_score=minimum_score,
        neighborhood=neighborhood,
        parsing_type=parsing_type,
        transpose=transpose,
        batch_id=batch_id,
    ) as stdout:
        try:
            reader = pacsv.open_csv(
                stdout, read_options=read_opts, parse_options=parse_opts, convert_options=convert_opts
            )
            for batch in reader:
                if batch.num_rows:
                    tables.append(pa.Table.from_batches([batch]))
        except pa.ArrowInvalid:
            pass

    return pa.concat_tables(tables).to_pandas() if tables else empty


def min_energy_by_pair_multi_cutoff(cutoffs) -> RisearchAggregation:
    """Off-target score for several cutoffs from ONE loose RIsearch pass.

    Run RIsearch at the loosest cutoff and, in the same streaming pass, keep the
    min-energy hit per (trigger, target) *separately for each cutoff*, counting a
    hit toward cutoff ``c`` only when ``score > c``. RIsearch's ``-s`` is exclusive
    and its cutoffs are nested, so this reproduces N independent per-cutoff runs
    bit-for-bit (see tests/complete/test_off_target_derivation_equivalence.py).

    Returns ``{cutoff: {(trigger, target): min_energy}}``. Memory is bounded by the
    number of distinct (cutoff, trigger, target) pairs, not the raw hit count.
    """
    sorted_cutoffs = sorted(set(int(c) for c in cutoffs))

    def _empty():
        import pyarrow as pa

        return pa.table(
            {
                "cutoff": pa.array([], pa.int64()),
                "trigger": pa.array([], pa.string()),
                "target": pa.array([], pa.string()),
                "energy": pa.array([], pa.float64()),
            }
        )

    def combine(batch):
        import pyarrow as pa
        import pyarrow.compute as pc

        table = pa.Table.from_batches([batch])
        parts = []
        for c in sorted_cutoffs:
            sub = table.filter(pc.greater(table.column("score"), c))
            if sub.num_rows == 0:
                continue
            agg = (
                sub.group_by(["trigger", "target"])
                .aggregate([("energy", "min")])
                .rename_columns(["trigger", "target", "energy"])
            )
            agg = agg.append_column("cutoff", pa.array([c] * agg.num_rows, pa.int64()))
            parts.append(agg.select(["cutoff", "trigger", "target", "energy"]))
        return pa.concat_tables(parts) if parts else _empty()

    def finalize(table):
        final = (
            table.group_by(["cutoff", "trigger", "target"])
            .aggregate([("energy", "min")])
            .rename_columns(["cutoff", "trigger", "target", "energy"])
        )
        result: dict = {c: {} for c in sorted_cutoffs}
        for c, trig, tgt, e in zip(
            final.column("cutoff").to_pylist(),
            final.column("trigger").to_pylist(),
            final.column("target").to_pylist(),
            final.column("energy").to_pylist(),
        ):
            result[c][(trig, tgt)] = float(e)
        return result

    return RisearchAggregation(columns=("trigger", "target", "score", "energy"), combine=combine, finalize=finalize)


def sum_exp_by_trigger_multi_cutoff(cutoffs, rt: float = 0.616) -> RisearchAggregation:
    """Single-target-gene score for several cutoffs from ONE loose RIsearch pass.

    A Boltzmann-style sum of exp(-rt*energy) over every hit per trigger (not just the
    strongest), kept separately for each cutoff — a hit counts toward cutoff ``c`` only
    when ``score > c``. RIsearch's ``-s`` is exclusive and its cutoffs are nested, so this
    reproduces N independent per-cutoff runs (within FP rounding, as the sum is
    float-order-dependent).

    Returns ``{cutoff: {trigger: sum_exp}}``.
    """
    sorted_cutoffs = sorted(set(int(c) for c in cutoffs))

    def _empty():
        import pyarrow as pa

        return pa.table(
            {
                "cutoff": pa.array([], pa.int64()),
                "trigger": pa.array([], pa.string()),
                "_exp": pa.array([], pa.float64()),
            }
        )

    def combine(batch):
        import pyarrow as pa
        import pyarrow.compute as pc

        exp_col = pc.exp(pc.multiply(batch.column("energy"), -rt))
        base = pa.table({"trigger": batch.column("trigger"), "score": batch.column("score"), "_exp": exp_col})
        parts = []
        for c in sorted_cutoffs:
            sub = base.filter(pc.greater(base.column("score"), c))
            if sub.num_rows == 0:
                continue
            agg = sub.group_by("trigger").aggregate([("_exp", "sum")]).rename_columns(["trigger", "_exp"])
            agg = agg.append_column("cutoff", pa.array([c] * agg.num_rows, pa.int64()))
            parts.append(agg.select(["cutoff", "trigger", "_exp"]))
        return pa.concat_tables(parts) if parts else _empty()

    def finalize(table):
        final = (
            table.group_by(["cutoff", "trigger"])
            .aggregate([("_exp", "sum")])
            .rename_columns(["cutoff", "trigger", "_exp"])
        )
        result: dict = {c: {} for c in sorted_cutoffs}
        for c, trig, s in zip(
            final.column("cutoff").to_pylist(),
            final.column("trigger").to_pylist(),
            final.column("_exp").to_pylist(),
        ):
            result[c][trig] = float(s)
        return result

    return RisearchAggregation(columns=("trigger", "score", "energy"), combine=combine, finalize=finalize)


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
