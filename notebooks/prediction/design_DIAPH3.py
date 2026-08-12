"""Generate v15 features for DIAPH3 ASOs and rank them for gymnotic delivery, cell line unknown.

DIAPH3 is looked up from the genome cache, so the tiled target is the full genomic span
(chr13:59,665,583-60,163,928 on the minus strand, 498,346 nt of pre-mRNA) and a 20-mer tiling
yields 498,327 candidates. That is ~57x the MALAT1 design, too large for one process, so this
script also runs partitioned: each shard tiles the whole gene, keeps a stride of the candidates,
and pushes them through the same pipeline `design_asos` runs (get_initial_data ->
_apply_standard_metadata -> generate_aso_features -> _fill_out_of_range_one_hots -> predict).
Every feature is computed per-candidate against the gene rather than against the scored batch, so
a shard reproduces the full-gene scores bit-for-bit and the merge is a concatenate-and-sort.

The screen context is deliberately unspecified: no cell line, no dose, no plating density. Those
leave `expr_*` and `off_target_score_specific_*` NaN, none of which the model guards as finite, so
the ranking stands on the sequence, structure and general off-target features.

Distributed (SLURM array of N tasks, then one merge job):

  python notebooks/prediction/design_DIAPH3.py --shard $SLURM_ARRAY_TASK_ID 16 --n-jobs 64
  python notebooks/prediction/design_DIAPH3.py --merge 16

Single process (small partial runs, e.g. a preflight canary):

  python notebooks/prediction/design_DIAPH3.py --first-n 200 --n-jobs 16

Shards write one parquet per block of candidates, so an interrupted run resumes by re-running the
same command: blocks already on disk are skipped unless --overwrite is passed.
"""
import argparse
import logging
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

from tauso.aso_generation import (
    AssetCache,
    Transfection,
    _apply_standard_metadata,
    _fill_out_of_range_one_hots,
    default_config,
    generate_aso_features,
    get_initial_data,
    summarize_design,
    tox_details,
)
from tauso.data.data import get_paths
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.inference import load_model, predict, score_column
from tauso.off_target.search import count_offtarget_matches_bulk

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE / "output" / "DIAPH3_nocell_gymnosis"

GENE = "DIAPH3"
GENOME = "GRCh38"
CELL_LINE = np.nan       # no cell line: the ranking is not conditioned on a screen context
DELIVERY = Transfection.GYMNOSIS
DOSE_NM = np.nan         # dose, nM
CELL_DENSITY = np.nan    # cells per well
CHEMISTRY = "2'MOE 5-10-5 gapmer, full PS"
ASO_SIZES = (20,)

# Bowtie match counts are reported for the best candidates only: the count pass aligns every query
# genome-wide at up to 2 mismatches reporting all hits, which does not scale to a 498k-row tiling.
OFFTARGET_TOP_N = 1000
OFFTARGET_MAX_MM = 2  # count genome matches up to 2 mismatches (Bowtie -v caps at 3)
OFFTARGET_COLS = ["perfect_matches", "off_targets_1mm", "off_targets_2mm"]
LIABILITY_COLS = ["offtarget_transcriptome", "offtarget_genomewide", "offtarget_rrna", "liabilities"]

SHARD_GLOB = "shard*_block*.parquet"


def build_config():
    cfg = default_config()                          # 2'-MOE 5-10-5 gapmer, full-PS 20-mer
    cfg.cell_line = CELL_LINE
    cfg.transfection_method = DELIVERY
    cfg.volume = DOSE_NM
    cfg.cell_per_well = CELL_DENSITY
    cfg.organism_name = "human"
    return cfg


def build_candidates(cache, config, first_n=None):
    """Tile the whole gene into candidates and stamp the standard metadata columns.

    Metadata is applied to the complete tiling before any sharding so `index_generated` numbers the
    candidates once, globally, and identifies a row across shards at merge time."""
    gene_data = cache.get_full_gene_data()
    if GENE not in gene_data:
        raise ValueError(f"gene {GENE!r} not found in the genome cache for {GENOME}.")
    gene_sequence = gene_data[GENE].full_mrna

    candidates = get_initial_data(gene_sequence, aso_sizes=list(ASO_SIZES), canonical_name=GENE)
    if first_n is not None:
        candidates = candidates.head(first_n).copy()
    _apply_standard_metadata(candidates, config)
    log(f"target {GENE}: {len(gene_sequence):,} nt -> {len(candidates):,} candidate ASOs {tuple(ASO_SIZES)}")
    return candidates


def score_block(block, cache, n_jobs, strict=True):
    """Featurize and score one block of candidates -- the tail of `design_asos`, on a subset."""
    featured, _ = generate_aso_features(block, cache, n_jobs=n_jobs)
    _, model_features = load_model()
    _fill_out_of_range_one_hots(featured, model_features)
    featured[score_column()] = predict(featured, strict=strict)
    return featured


def log(message):
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def _fmt(seconds):
    return f"{int(seconds) // 3600}h{int(seconds) % 3600 // 60:02d}m{int(seconds) % 60:02d}s"


def run_shard(args, out_dir):
    """Featurize + score this shard's stride of the tiling, one parquet per block."""
    cache = AssetCache(genome=GENOME)
    candidates = build_candidates(cache, build_config(), first_n=args.first_n)
    mine = candidates.iloc[args.shard :: args.shards]
    blocks = [mine.iloc[i : i + args.block_size].copy() for i in range(0, len(mine), args.block_size)]
    tag = f"shard {args.shard}/{args.shards}"
    log(f"{tag}: {len(mine):,} of {len(candidates):,} candidates in {len(blocks)} blocks "
        f"of <= {args.block_size:,}, n_jobs={args.n_jobs}")

    shard_dir = out_dir / "shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    started, done = time.time(), 0
    for i, block in enumerate(blocks, start=1):
        path = shard_dir / f"shard{args.shard:02d}_block{i:03d}.parquet"
        if path.exists() and not args.overwrite:
            done += len(block)
            log(f"{tag}: block {i}/{len(blocks)} already on disk, skipping ({done:,}/{len(mine):,})")
            continue
        block_started = time.time()
        log(f"{tag}: block {i}/{len(blocks)} starting, {len(block):,} candidates")
        scored = score_block(block, cache, n_jobs=args.n_jobs, strict=args.strict)
        scored.to_parquet(path, index=False)
        done += len(block)
        elapsed = time.time() - started
        eta = elapsed / done * (len(mine) - done)
        log(f"{tag}: block {i}/{len(blocks)} done in {_fmt(time.time() - block_started)} "
            f"-- {done:,}/{len(mine):,} candidates, elapsed {_fmt(elapsed)}, eta {_fmt(eta)}")

    log(f"{tag}: COMPLETE, {done:,} candidates in {_fmt(time.time() - started)} -> {shard_dir}")


def add_offtarget_counts(summary, top_n):
    """Add per-ASO off-target counts for the top `top_n` candidates: `perfect_matches` (0 mismatches),
    `off_targets_1mm`, and `off_targets_2mm` -- genome-wide Bowtie match counts at 0/1/2 mismatches
    (both strands, one pass), EXCLUDING any hit inside the on-target DIAPH3 locus so the intended site
    and intragenic near-matches are not counted. Candidates below the cut keep <NA>. If the Bowtie index
    is missing the columns are filled with <NA> and a note printed, rather than triggering a multi-GB
    index build (run `tauso setup-bowtie --genome GRCh38`)."""
    for col in OFFTARGET_COLS:
        summary[col] = pd.NA

    sentinel = Path(get_paths(GENOME)["fasta"]).parent / f"{GENOME}_bowtie_index" / "SUCCESS"
    if not sentinel.exists():
        log(f"Bowtie index for {GENOME} not found; skipping off-target counts "
            f"(run `tauso setup-bowtie --genome {GENOME}`).")
        return summary

    head = summary.head(top_n)
    log(f"counting genome-wide off-target matches for the top {len(head):,} candidates "
        f"(<= {OFFTARGET_MAX_MM} mismatches)")
    g = get_locus_to_data_dict(include_introns=True, gene_subset=[GENE], genome=GENOME)[GENE]
    counts = count_offtarget_matches_bulk(
        head["aso_sequence"].tolist(), genome=GENOME,
        max_mismatches=OFFTARGET_MAX_MM, exclude_regions=[(g.chrom, g.gene_start, g.gene_end)],
    )
    for position, col in enumerate(OFFTARGET_COLS):
        summary.loc[head.index, col] = head["aso_sequence"].map(lambda s: counts[s][position])
    return summary


def run_merge(args, out_dir):
    """Concatenate every shard block, rank globally, and write the features + consumer summary."""
    shard_dir = out_dir / "shards"
    paths = sorted(shard_dir.glob(SHARD_GLOB))
    if not paths:
        raise FileNotFoundError(f"no shard blocks found under {shard_dir}")
    shards_seen = sorted({int(p.name.split("_")[0][len("shard"):]) for p in paths})
    if args.merge is not None and len(shards_seen) != args.merge:
        missing = sorted(set(range(args.merge)) - set(shards_seen))
        raise ValueError(f"expected {args.merge} shards, found {len(shards_seen)}; missing shards {missing}")
    log(f"merging {len(paths)} shard blocks from {shard_dir}")

    frames = []
    for i, path in enumerate(paths, start=1):
        frames.append(pd.read_parquet(path))
        if i % 10 == 0 or i == len(paths):
            log(f"read {i}/{len(paths)} blocks, {sum(len(f) for f in frames):,} rows so far")
    merged = pd.concat(frames, ignore_index=True)
    del frames

    before = len(merged)
    merged = merged.drop_duplicates(subset=["index_generated"])
    if len(merged) != before:
        log(f"dropped {before - len(merged):,} duplicate candidates (overlapping shard blocks)")

    col = score_column()
    ranked = merged.sort_values(col, ascending=False, kind="stable").reset_index(drop=True)
    log(f"ranked {len(ranked):,} candidates by {col}")

    summary = summarize_design(ranked)
    summary = add_offtarget_counts(summary, args.offtarget_top_n)
    details = tox_details(ranked)                       # row-aligned to ranked -> summarize_design
    for col_name in LIABILITY_COLS:
        summary[col_name] = details[col_name].to_numpy()
    summary["chemistry"] = CHEMISTRY
    summary["transfection_method"] = DELIVERY
    summary["dosage_nm"] = DOSE_NM
    summary["cell_line"] = CELL_LINE
    summary["density_cells_per_well"] = CELL_DENSITY

    out_dir.mkdir(parents=True, exist_ok=True)
    ranked.to_parquet(out_dir / "DIAPH3_nocell_gymnosis_features.parquet", index=False)
    summary.to_csv(out_dir / "DIAPH3_nocell_gymnosis_ranked.csv", index=False)
    log(f"{len(ranked):,} {GENE} ASOs featurized + ranked "
        f"(off-target counts for the top {args.offtarget_top_n:,} + liabilities) -> {out_dir}")


def run_single(args, out_dir):
    """Featurize, score and rank in one process -- the design_MALAT1 flow, for small runs."""
    cache = AssetCache(genome=GENOME)
    candidates = build_candidates(cache, build_config(), first_n=args.first_n)
    started = time.time()
    ranked = score_block(candidates, cache, n_jobs=args.n_jobs, strict=args.strict)
    ranked = ranked.sort_values(score_column(), ascending=False, kind="stable").reset_index(drop=True)
    log(f"scored {len(ranked):,} candidates in {_fmt(time.time() - started)}")

    summary = summarize_design(ranked)
    summary = add_offtarget_counts(summary, args.offtarget_top_n)
    details = tox_details(ranked)
    for col_name in LIABILITY_COLS:
        summary[col_name] = details[col_name].to_numpy()
    summary["chemistry"] = CHEMISTRY
    summary["transfection_method"] = DELIVERY
    summary["dosage_nm"] = DOSE_NM
    summary["cell_line"] = CELL_LINE
    summary["density_cells_per_well"] = CELL_DENSITY

    out_dir.mkdir(parents=True, exist_ok=True)
    ranked.to_parquet(out_dir / "DIAPH3_nocell_gymnosis_features.parquet", index=False)
    summary.to_csv(out_dir / "DIAPH3_nocell_gymnosis_ranked.csv", index=False)
    log(f"{len(ranked):,} {GENE} ASOs featurized + ranked -> {out_dir}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shard", nargs=2, type=int, metavar=("K", "N"),
                    help="featurize shard K of N (candidates[K::N]) and write its block parquets")
    ap.add_argument("--merge", type=int, metavar="N",
                    help="merge the blocks of all N shards, rank globally, and write the outputs")
    ap.add_argument("--first-n", type=int, default=None, help="tile only the first N candidates")
    ap.add_argument("--block-size", type=int, default=8000,
                    help="candidates per block within a shard; one parquet + one progress line each")
    ap.add_argument("--n-jobs", type=int, default=min(os.cpu_count() or 1, 24), help="worker processes")
    ap.add_argument("--offtarget-top-n", type=int, default=OFFTARGET_TOP_N,
                    help="how many top-ranked candidates get genome-wide Bowtie match counts")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output directory")
    ap.add_argument("--overwrite", action="store_true", help="recompute blocks that are already on disk")
    ap.add_argument("--no-strict", dest="strict", action="store_false",
                    help="score candidates whose guarded features failed instead of raising")
    args = ap.parse_args()

    if args.shard is not None and args.merge is not None:
        ap.error("--shard and --merge are separate steps; pass one at a time")
    if args.shard is not None:
        args.shard, args.shards = args.shard
        if not 0 <= args.shard < args.shards:
            ap.error(f"shard index {args.shard} out of range for {args.shards} shards")

    # The pipeline reports its progress through the logging module: one line per calculator step,
    # which is the only visibility into a block that runs for minutes.
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s",
                        datefmt="%H:%M:%S", force=True)

    log(f"{GENE} design: cell_line={CELL_LINE} delivery={DELIVERY} dose={DOSE_NM} "
        f"density={CELL_DENSITY} chemistry={CHEMISTRY!r} out={args.out}")
    if args.shard is not None:
        run_shard(args, args.out)
    elif args.merge is not None:
        run_merge(args, args.out)
    else:
        run_single(args, args.out)


if __name__ == "__main__":
    main()
