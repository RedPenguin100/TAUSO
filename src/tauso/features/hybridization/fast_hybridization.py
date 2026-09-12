import os
import platform
import tempfile
from pathlib import Path
from typing import Callable, Dict, Tuple

import pyrisearch_tauso

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


# What an ASO forms with its target: an RNA:DNA duplex, wobble pairs left out.
ASO_TARGET_MATRIX = pyrisearch_tauso.Matrix.SU95_NO_GU

# Per-nucleotide extension penalty, in dacal/mol. RIsearch leaves this at 0;
# hybridization here is calibrated with it set.
EXTENSION_PENALTY = 30

# The columns RIsearch's hits come back in.
RISEARCH_COLUMNS = pyrisearch_tauso.HIT_COLUMNS


# A reduction over RIsearch hits: which columns to read, how to reduce a batch
# to a partial, and how to reduce the partials. `empty` is what no hits means.
RisearchAggregation = pyrisearch_tauso.Reduction


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

    return RisearchAggregation(columns=(*keys, "score", "energy"), combine=combine, finalize=finalize, empty={})


def aggregate_by_pair_multi_cutoff(cutoffs) -> RisearchAggregation:
    """Per-(query, target) Boltzmann reducer for several cutoffs from a loose RIsearch pass.

    Keeps ``sum(exp(-energy/RT_KCAL_MOL))`` per pair (always >= 0).

    Returns ``{cutoff: {(query, target): sum_exp}}``.
    """
    return _by_key_multi_cutoff(
        cutoffs,
        keys=("query", "target"),
        value_aggs=(("_exp", "_exp", "sum"),),
        pack=lambda r: float(r["_exp"]),
    )


def stats_by_query_multi_cutoff(cutoffs) -> RisearchAggregation:
    """Site-resolved stats per query, per cutoff: the Boltzmann sum, the best-site energy
    and the hit count -- a strict superset of the per-query Boltzmann sum alone.

    Callers derive the effective number of sites (target multiplicity)
    ``log_eff = log(sum_exp) - (-min_energy / RT_KCAL_MOL)`` -- 0 for a single dominant site,
    growing when several comparable sites share the binding.

    Returns ``{cutoff: {query: (sum_exp, min_energy, n_sites)}}``.
    """
    return _by_key_multi_cutoff(
        cutoffs,
        keys=("query",),
        value_aggs=(("s", "_exp", "sum"), ("e", "energy", "min"), ("n", "energy", "count")),
        pack=lambda r: (float(r["s"]), float(r["e"]), int(r["n"])),
    )
