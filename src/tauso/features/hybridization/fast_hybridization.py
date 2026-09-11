import os
import platform
import re
import tempfile
import uuid
from pathlib import Path
from typing import Callable, Dict, List, NamedTuple, Tuple

import pyrisearch_tauso

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


def _interaction_mode(interaction_type: Interaction) -> str:
    if interaction_type == Interaction.RNA_DNA_NO_WOBBLE:
        return "su95_noGU"
    if interaction_type == Interaction.RNA_RNA:
        return "t04"
    raise ValueError(f"Unsupported interaction type: {interaction_type}")


# Per-nucleotide extension penalty, in dacal/mol. RIsearch leaves this at 0;
# hybridization here is calibrated with it set.
EXTENSION_PENALTY = 30

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


# TAUSO names the hit columns for what they mean here; pyrisearch_tauso names
# them after the query and target sides of the search.
_LIBRARY_COLUMN = {
    "trigger": "qname",
    "trigger_start": "qbeg",
    "trigger_end": "qend",
    "target": "tname",
    "target_start": "tbeg",
    "target_end": "tend",
    "score": "score",
    "energy": "energy",
}


def _query_sequences(trigger_id_seq_pairs: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """The triggers as RIsearch has to see them: a trigger binds its target, so
    what is searched for is the antisense of it."""
    return [(query_id, get_antisense_rna(trigger)) for query_id, trigger in trigger_id_seq_pairs]


def _under_tauso_names(batch, columns):
    """A batch of hits under the names the aggregations read.

    pyarrow hands back the columns in the order they were asked for, so the
    names line up by position.
    """
    import pyarrow as pa

    return pa.RecordBatch.from_arrays(list(batch.columns), names=list(columns))


def _search_options(
    trigger_id_seq_pairs,
    target_file_path,
    *,
    interaction_type,
    minimum_score,
    neighborhood,
    transpose,
    block_size,
):
    """What a RIsearch search takes, under the options TAUSO runs it with."""
    return dict(
        queries=_query_sequences(trigger_id_seq_pairs),
        targets=target_file_path,
        min_score=minimum_score,
        matrix=_interaction_mode(interaction_type),
        extension_penalty=EXTENSION_PENALTY,
        neighborhood=neighborhood,
        transpose=transpose,
        block_size=block_size,
    )


def _library_reduction(aggregation: "RisearchAggregation"):
    """The aggregation as pyrisearch_tauso reads it.

    The library names the hit columns after the query and target sides of the
    search, so each batch is renamed before an aggregation written against
    TAUSO's names sees it.
    """
    return pyrisearch_tauso.Reduction(
        columns=[_LIBRARY_COLUMN[column] for column in aggregation.columns],
        combine=lambda batch: aggregation.combine(_under_tauso_names(batch, aggregation.columns)),
        finalize=aggregation.finalize,
        empty={},
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

    args = [
        "-q",
        str(query_path),
        "-t",
        str(target_path),
        "-s",
        str(minimum_score),
        "-d",
        str(EXTENSION_PENALTY),
        "-m",
        _interaction_mode(interaction_type),
        "-n",
        str(neighborhood),
    ]
    if transpose:
        args.append("-R")
    if parsing_type is not None:
        args.append(f"-p{parsing_type}")

    try:
        return pyrisearch_tauso.run(args, cwd=TMP_PATH).stdout
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

    See aggregate_by_pair_multi_cutoff and stats_by_trigger_multi_cutoff.
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

    return pyrisearch_tauso.search_reduced(
        reduction=_library_reduction(aggregation),
        **_search_options(
            trigger_id_seq_pairs,
            target_file_path,
            interaction_type=interaction_type,
            minimum_score=minimum_score,
            neighborhood=neighborhood,
            transpose=transpose,
            block_size=block_size,
        ),
    )


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
    """Run a batched RIsearch and return ALL hits as a DataFrame (the 8 RISEARCH_COLUMNS).

    Same batched invocation as parse_risearch_hits_pyarrow and streamed through
    pyarrow, but materialises the full table instead of reducing it — intended for tests
    and debugging, not the hot path (production uses parse_risearch_hits_pyarrow with an
    aggregation so memory stays bounded). Replaces the old get_triggers_mfe_scores_batch
    + parse_risearch_output pair.
    """
    import pandas as pd

    if not trigger_id_seq_pairs:
        return pd.DataFrame(columns=list(RISEARCH_COLUMNS))

    table = pyrisearch_tauso.hits_table(
        columns=[_LIBRARY_COLUMN[column] for column in RISEARCH_COLUMNS],
        **_search_options(
            trigger_id_seq_pairs,
            target_file_path,
            interaction_type=interaction_type,
            minimum_score=minimum_score,
            neighborhood=neighborhood,
            transpose=transpose,
            block_size=block_size,
        ),
    )
    return table.rename_columns(list(RISEARCH_COLUMNS)).to_pandas()


# ---------------------------------------------------------------------------
# Multi-cutoff aggregations
#
# Each factory returns a RisearchAggregation that reduces a loose RIsearch pass
# into ``{cutoff: {key: value}}``. They share one streaming skeleton,
# ``_by_key_multi_cutoff``: group the hits with ``score > cutoff`` by ``keys``,
# aggregate, and pack each group into a stored value. Deriving every cutoff from a
# loose run is equivalent to independent per-cutoff runs, because RIsearch's
# ``-s`` is exclusive and its cutoffs are nested (see
# tests/complete/test_off_target_derivation_equivalence.py); Boltzmann sums match
# only within FP rounding, being float-order dependent.
# ---------------------------------------------------------------------------

# RT at body temperature (310 K), in kcal/mol -- the kT of the Boltzmann factor
# exp(-energy/RT). Single source of truth for hybridization; passed to no function.
RT_KCAL_MOL = 0.616

# pyarrow aggregate func -> the func that merges per-block partials in finalize.
_MERGE_FUNC = {"sum": "sum", "min": "min", "count": "sum"}


def _by_key_multi_cutoff(
    cutoffs,
    *,
    keys: Tuple[str, ...],
    value_aggs: Tuple[Tuple[str, str, str], ...],
    pack: Callable,
) -> RisearchAggregation:
    """Streaming per-cutoff, per-key RIsearch reducer shared by every aggregation factory.

    keys       group-by columns; the result key is that single value or a tuple of them.
    value_aggs ``(name, source, func)`` triples aggregating ``source`` with pyarrow ``func``
               into a partial column ``name``. ``source`` is ``energy`` or ``_exp`` (the
               Boltzmann weight ``exp(-energy/RT_KCAL_MOL)``, materialized only when some
               triple sources ``_exp``).
    pack       maps a finalized group's ``{name: value}`` to the stored value.

    Returns a RisearchAggregation producing ``{cutoff: {key: pack(...)}}``.
    """
    sorted_cutoffs = sorted({int(c) for c in cutoffs})
    keys = tuple(keys)
    names = [name for name, _src, _func in value_aggs]
    needs_energy = any(src == "energy" for _n, src, _f in value_aggs)
    needs_exp = any(src == "_exp" for _n, src, _f in value_aggs)

    def _empty_partial():
        import pyarrow as pa

        cols = {"cutoff": pa.array([], pa.int64())}
        for k in keys:
            cols[k] = pa.array([], pa.string())
        for name, _src, func in value_aggs:
            cols[name] = pa.array([], pa.int64() if func == "count" else pa.float64())
        return pa.table(cols)

    def combine(batch):
        import pyarrow as pa
        import pyarrow.compute as pc

        cols = {k: batch.column(k) for k in keys}
        cols["score"] = batch.column("score")
        if needs_exp:
            cols["_exp"] = pc.exp(pc.divide(batch.column("energy"), -RT_KCAL_MOL))
        if needs_energy:
            cols["energy"] = batch.column("energy")
        base = pa.table(cols)

        parts = []
        for c in sorted_cutoffs:
            sub = base.filter(pc.greater(base.column("score"), c))
            if sub.num_rows == 0:
                continue
            agg = sub.group_by(list(keys)).aggregate([(src, func) for _n, src, func in value_aggs])
            out = {"cutoff": pa.array([c] * agg.num_rows, pa.int64())}
            for k in keys:
                out[k] = agg.column(k)
            for name, src, func in value_aggs:
                out[name] = agg.column(f"{src}_{func}")
            parts.append(pa.table(out))
        return pa.concat_tables(parts) if parts else _empty_partial()

    def finalize(table):
        merged = table.group_by(["cutoff", *keys]).aggregate(
            [(name, _MERGE_FUNC[func]) for name, _src, func in value_aggs]
        )
        cutoff_col = merged.column("cutoff").to_pylist()
        key_cols = [merged.column(k).to_pylist() for k in keys]
        value_cols = {name: merged.column(f"{name}_{_MERGE_FUNC[func]}").to_pylist() for name, _src, func in value_aggs}
        result: dict = {c: {} for c in sorted_cutoffs}
        for i, c in enumerate(cutoff_col):
            key = key_cols[0][i] if len(keys) == 1 else tuple(col[i] for col in key_cols)
            result[c][key] = pack({name: value_cols[name][i] for name in names})
        return result

    return RisearchAggregation(columns=(*keys, "score", "energy"), combine=combine, finalize=finalize)


def aggregate_by_pair_multi_cutoff(cutoffs) -> RisearchAggregation:
    """Per-(trigger, target) Boltzmann reducer for several cutoffs from a loose RIsearch pass.

    Keeps ``sum(exp(-energy/RT_KCAL_MOL))`` per pair (always >= 0).

    Returns ``{cutoff: {(trigger, target): sum_exp}}``.
    """
    return _by_key_multi_cutoff(
        cutoffs,
        keys=("trigger", "target"),
        value_aggs=(("_exp", "_exp", "sum"),),
        pack=lambda r: float(r["_exp"]),
    )


def stats_by_trigger_multi_cutoff(cutoffs) -> RisearchAggregation:
    """Site-resolved stats per trigger, per cutoff: the Boltzmann sum, the best-site energy
    and the hit count -- a strict superset of the per-trigger Boltzmann sum alone.

    Callers derive the effective number of sites (target multiplicity)
    ``log_eff = log(sum_exp) - (-min_energy / RT_KCAL_MOL)`` -- 0 for a single dominant site,
    growing when several comparable sites share the binding.

    Returns ``{cutoff: {trigger: (sum_exp, min_energy, n_sites)}}``.
    """
    return _by_key_multi_cutoff(
        cutoffs,
        keys=("trigger",),
        value_aggs=(("s", "_exp", "sum"), ("e", "energy", "min"), ("n", "energy", "count")),
        pack=lambda r: (float(r["s"]), float(r["e"]), int(r["n"])),
    )


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
