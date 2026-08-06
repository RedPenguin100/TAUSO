"""Per-gene RIsearch scan engine.

Scores ASOs against their target gene(s) with one loose RIsearch pass per gene, in parallel: a
target FASTA per gene, ASOs split into worker chunks over a thread pool, hits reduced to
site-resolved ``SiteStats`` (Boltzmann occupancy, best-site energy, hit count) per ASO and cutoff.
That superset feeds every per-gene hybridization feature, which derive their columns via
``emit_site_columns`` (see off_target_specific_gene.py).
"""

import logging
import os
import uuid
from collections import defaultdict
from typing import NamedTuple

import pandas as pd

from ....data.consts import ASO_SEQUENCE
from ....util import get_antisense
from ..fast_hybridization import (
    TMP_PATH,
    Interaction,
    dump_target_file,
    parse_risearch_hits_pyarrow,
    stats_by_trigger_multi_cutoff,
)
from .parallel import run_tasks_parallel

logger = logging.getLogger(__name__)


class SiteStats(NamedTuple):
    """Site-resolved reduction of one ASO's hits against one gene, at one cutoff."""

    sum_exp: float  # Sum of exp(-energy/RT) over sites (Boltzmann occupancy, >= 0)
    min_energy: float  # best (most negative) single-site energy
    n_sites: int  # number of hits above the cutoff


def _validate_genes_found(target_genes, gene_to_data):
    not_found = [g for g in target_genes if g not in gene_to_data]
    if not_found:
        raise ValueError(f"The following genes are not found in gene_to_data: {not_found}")


def _scan_one_gene_chunk(row_triggers, target_path, cutoffs, chunk_size=100):
    """Site-resolved stats for one gene's ASOs against its target, in query chunks, per cutoff.

    row_triggers: [(aso_index, trigger_seq)] for ASOs targeting this gene.
    Returns ``{cutoff: {trigger_id: (sum_exp, min_energy, n_sites)}}``. chunk_size caps the
    per-call query batch (and thus peak RIsearch output held while parsing).
    """
    cutoffs = [int(c) for c in cutoffs]
    minimum_score = min(cutoffs)
    combined: dict = {c: {} for c in cutoffs}
    for i in range(0, len(row_triggers), chunk_size):
        chunk = row_triggers[i : i + chunk_size]
        partial = parse_risearch_hits_pyarrow(
            trigger_id_seq_pairs=[(str(idx), trig) for idx, trig in chunk],
            target_file_path=target_path,
            aggregation=stats_by_trigger_multi_cutoff(cutoffs),
            minimum_score=minimum_score,
            parsing_type="2",
            interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            transpose=True,
            batch_id=f"{os.getpid()}-{uuid.uuid4().hex}",
        )
        for c in cutoffs:
            dst = combined[c]
            for trigger, (s, e, n) in partial.get(c, {}).items():
                if trigger in dst:
                    ps, pe, pn = dst[trigger]
                    dst[trigger] = (ps + s, min(pe, e), pn + n)
                else:
                    dst[trigger] = (s, e, n)
    return combined


def build_gene_chunk_tasks(gene_to_row_triggers, gene_to_target_path, n_jobs):
    """Split each gene's ASOs into sub-batches so all n_jobs workers are used.

    When there are fewer genes than workers, each gene's rows are divided into
    ceil(n_jobs / n_genes) sub-batches. Returns [(gene, sub_triggers, target_path)].
    """
    n_genes = max(1, len(gene_to_row_triggers))
    tasks_per_gene = max(1, (n_jobs + n_genes - 1) // n_genes)
    tasks = []
    for gene, row_triggers in gene_to_row_triggers.items():
        sub_size = max(1, (len(row_triggers) + tasks_per_gene - 1) // tasks_per_gene)
        target_path = gene_to_target_path[gene]
        for i in range(0, len(row_triggers), sub_size):
            tasks.append((gene, row_triggers[i : i + sub_size], target_path))
    return tasks


def scan_gene_sites(aso_df, gene_to_data, target_genes, get_gene_fn, cutoffs, n_jobs):
    """Run the per-gene RIsearch scan and return site-resolved stats per ASO and cutoff.

    get_gene_fn(row) returns the target gene for a row (its own canonical gene for on-target, or
    a fixed gene for single off-target). Returns ``{cutoff: {row_index: SiteStats}}`` covering
    every ASO that had at least one hit above the cutoff against its target.
    """
    cutoffs = [int(c) for c in cutoffs]
    _validate_genes_found(target_genes, gene_to_data)
    TMP_PATH.mkdir(exist_ok=True)

    gene_to_target_path = {}
    try:
        for gene in target_genes:
            gene_to_target_path[gene] = dump_target_file(
                f"target-{gene}-{uuid.uuid4().hex}.fa", {gene: gene_to_data[gene].full_mrna}
            )

        gene_to_row_triggers = defaultdict(list)
        for idx, seq, gene in zip(aso_df.index, aso_df[ASO_SEQUENCE], aso_df.apply(get_gene_fn, axis=1)):
            if pd.notna(gene) and gene in gene_to_target_path:
                gene_to_row_triggers[gene].append((idx, get_antisense(seq)))

        chunk_tasks = build_gene_chunk_tasks(gene_to_row_triggers, gene_to_target_path, n_jobs)
        logger.info(
            "RIsearch scan: %d genes, %d ASOs, %d tasks, %d cutoffs",
            len(gene_to_row_triggers),
            sum(len(v) for v in gene_to_row_triggers.values()),
            len(chunk_tasks),
            len(cutoffs),
        )

        tasks = [((gene, i), triggers, path, cutoffs) for i, (gene, triggers, path) in enumerate(chunk_tasks)]
        chunk_stats = run_tasks_parallel(tasks, _scan_one_gene_chunk, n_jobs)

        gene_stats: dict = defaultdict(lambda: {c: {} for c in cutoffs})
        for (gene, _), per_cutoff in chunk_stats.items():
            for c in cutoffs:
                gene_stats[gene][c].update(per_cutoff.get(c, {}))

        scan: dict = {c: {} for c in cutoffs}
        for gene, row_triggers in gene_to_row_triggers.items():
            for c in cutoffs:
                gene_cutoff_stats = gene_stats[gene][c]
                for idx, _ in row_triggers:
                    stat = gene_cutoff_stats.get(str(idx))
                    if stat is not None:
                        scan[c][idx] = SiteStats(*stat)
    finally:
        for path in gene_to_target_path.values():
            if os.path.exists(path):
                os.remove(path)

    return scan


def emit_site_columns(aso_df, scan, cutoffs, derivations):
    """Write one column per (derivation, cutoff) from a scan_gene_sites result.

    derivations: [(name_fn, derive_fn)]. For each cutoff, name_fn(cutoff) names the column and
    derive_fn(SiteStats) -> value; rows with no scored hit (or a None derivation) default to 0.0.
    Returns (aso_df, [feature_names]).
    """
    cutoffs = [int(c) for c in cutoffs]
    feature_names = []
    for cutoff in cutoffs:
        per_row = scan.get(cutoff, {})
        for name_fn, derive_fn in derivations:
            values = pd.Series(0.0, index=aso_df.index)
            for idx, stats in per_row.items():
                value = derive_fn(stats)
                if value is not None:
                    values[idx] = value
            name = name_fn(cutoff)
            aso_df[name] = values
            feature_names.append(name)
    return aso_df, feature_names
