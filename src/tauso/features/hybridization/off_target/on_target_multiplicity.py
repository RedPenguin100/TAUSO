"""Site-resolved on-target target-multiplicity feature.

Where ``on_target_total_hybridization`` keeps only the Boltzmann sum over an ASO's on-target
sites (dominated by the single strongest site), this scores the same RIsearch pass site-resolved:
per ASO and cutoff it keeps ``sum(exp(-E/RT))``, the best-site energy ``min(E)`` and the hit count,
then derives the log effective number of sites ``log_eff = log(sum_exp) - (-min_E / RT)`` -- 0 for a
single dominant site, growing when several comparable sites share the binding.

Reuses the per-gene chunking / parallel dispatch of ``gene_chunk_scoring``; only the aggregation
(three stats instead of one) and the final derivation differ.
"""

import logging
import os
import uuid
from collections import defaultdict

import numpy as np
import pandas as pd

from ....data.consts import ASO_SEQUENCE, CANONICAL_GENE_NAME
from ....util import get_antisense
from ..fast_hybridization import (
    TMP_PATH,
    Interaction,
    dump_target_file,
    parse_risearch_hits_pyarrow,
    stats_by_trigger_multi_cutoff,
)
from .gene_chunk_scoring import _RT, _validate_genes_found, build_gene_chunk_tasks
from .parallel import run_tasks_parallel

logger = logging.getLogger(__name__)


def _chunkify_stats_one_gene(row_triggers, target_path, cutoffs, chunk_size=100):
    """Site-resolved stats for one gene's ASOs against its target, in query chunks, per cutoff.

    Returns ``{cutoff: {trigger_id: (sum_exp, min_energy, n_sites)}}``.
    """
    cutoffs = [int(c) for c in cutoffs]
    minimum_score = min(cutoffs)
    combined: dict = {c: {} for c in cutoffs}
    for i in range(0, len(row_triggers), chunk_size):
        chunk = row_triggers[i : i + chunk_size]
        partial = parse_risearch_hits_pyarrow(
            trigger_id_seq_pairs=[(str(idx), trig) for idx, trig in chunk],
            target_file_path=target_path,
            aggregation=stats_by_trigger_multi_cutoff(cutoffs, RT=_RT),
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


def on_target_log_number_of_sites(aso_df, gene_to_data, cutoffs, n_jobs=1):
    """Log effective number of on-target sites per ASO, one column per cutoff.

    Each oligo is scored against its own canonical gene. Returns (aso_df, [feature_names]) with
    columns ``on_target_log_number_of_sites_{cutoff}``. ASOs with no gene or no qualifying hit get 0.
    """
    cutoffs = [int(c) for c in cutoffs]
    target_genes = aso_df[CANONICAL_GENE_NAME].dropna().unique()
    _validate_genes_found(target_genes, gene_to_data)
    TMP_PATH.mkdir(exist_ok=True)

    gene_to_target_path = {}
    try:
        for gene in target_genes:
            gene_to_target_path[gene] = dump_target_file(
                f"ontgt-{gene}-{uuid.uuid4().hex}.fa", {gene: gene_to_data[gene].full_mrna}
            )

        gene_to_row_triggers = defaultdict(list)
        for idx, seq, gene in zip(aso_df.index, aso_df[ASO_SEQUENCE], aso_df[CANONICAL_GENE_NAME]):
            if pd.notna(gene) and gene in gene_to_target_path:
                gene_to_row_triggers[gene].append((idx, get_antisense(seq)))

        chunk_tasks = build_gene_chunk_tasks(gene_to_row_triggers, gene_to_target_path, n_jobs)
        logger.info(
            "On-target multiplicity: %d genes, %d ASOs, %d tasks, %d cutoffs",
            len(gene_to_row_triggers),
            sum(len(v) for v in gene_to_row_triggers.values()),
            len(chunk_tasks),
            len(cutoffs),
        )

        tasks = [((gene, i), triggers, path, cutoffs) for i, (gene, triggers, path) in enumerate(chunk_tasks)]
        chunk_stats = run_tasks_parallel(tasks, _chunkify_stats_one_gene, n_jobs)

        gene_stats: dict = defaultdict(lambda: {c: {} for c in cutoffs})
        for (gene, _), per_cutoff in chunk_stats.items():
            for c in cutoffs:
                gene_stats[gene][c].update(per_cutoff.get(c, {}))

        feature_names = []
        for cutoff in cutoffs:
            values = pd.Series(0.0, index=aso_df.index)
            for gene, row_triggers in gene_to_row_triggers.items():
                gene_cutoff_stats = gene_stats[gene][cutoff]
                for idx, _ in row_triggers:
                    stat = gene_cutoff_stats.get(str(idx))
                    if stat is None:
                        continue
                    s, e, _n = stat
                    if s > 0:
                        log_eff = np.log(max(s, 1e-300)) - (-e / _RT)
                        if np.isfinite(log_eff):
                            values[idx] = log_eff
            fname = f"on_target_log_number_of_sites_{cutoff}"
            aso_df[fname] = values
            feature_names.append(fname)
    finally:
        for path in gene_to_target_path.values():
            if os.path.exists(path):
                os.remove(path)

    return aso_df, feature_names
