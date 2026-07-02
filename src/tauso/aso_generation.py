import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

from tauso.data.consts import (
    ASO_SEQUENCE,
    CANONICAL_GENE_NAME,
    CELL_LINE,
    CELL_LINE_DEPMAP,
    CELL_LINE_DEPMAP_PROXY,
    CELL_LINE_ORGANISM,
    CHEMICAL_PATTERN,
    DENSITY_CELLS_PER_WELL,
    MODIFICATION_STRING,
    PS_PATTERN,
    STRUCTURE_SENSE_LENGTH,
    TRANSFECTION_RAW,
    VOLUME_NM,
    resolve_depmap_id,
    resolve_depmap_proxy,
)
from tauso.populate.calculators.cache import AssetCache
from tauso.populate.calculators.calculator import Calculator
from tauso.util import get_antisense


def get_initial_data(target_mrna, aso_sizes, canonical_name):
    candidates = []
    sense_lengths = []

    for aso_size in aso_sizes:
        for i in range(0, len(target_mrna) - (aso_size - 1)):
            target = target_mrna[i : i + aso_size]
            candidates.append(get_antisense(str(target)))
            sense_lengths.append(aso_size)

    df = pd.DataFrame(
        {
            ASO_SEQUENCE: candidates,
            STRUCTURE_SENSE_LENGTH: sense_lengths,
            CANONICAL_GENE_NAME: canonical_name,
        }
    )
    return df


def get_organism_name(genome):
    if genome == "GRCh38":
        return "human"
    elif genome == "GRCm39":
        return "mouse"
    else:
        raise ValueError(f"Unknown genome type: {genome}")


class Transfection:
    ELECTROPORATION = "Electroporation"
    LIPOFECTION = "Lipofection"
    GYMNOSIS = "Gymnosis"
    OTHER = "Other"


class Config:
    cell_per_well: int
    volume: int

    cell_line: str
    transfection_method: str

    standard_chemical_pattern: str
    standard_ps_pattern: str
    standard_modification: str

    organism_name: str


def default_config():
    data_config = Config()
    data_config.cell_per_well = 10000
    data_config.volume = 100
    data_config.transfection_method = Transfection.GYMNOSIS
    data_config.standard_modification = "MOE/5-methylcytosines/deoxy"
    data_config.standard_ps_pattern = 19 * "*"
    data_config.standard_chemical_pattern = "MMMMMddddddddddMMMMM"

    data_config.cell_line = "T24"

    data_config.organism_name = "human"

    return data_config


def generate_stub_data(
    target_gene: str,
    gene_sequence: str,
    first_n: int = None,
    version: str = None,
    data_config: Config = None,
):
    data = get_initial_data(gene_sequence, aso_sizes=[20], canonical_name=target_gene)

    if data_config is None:
        data_config = default_config()

    if first_n is not None:
        data = data[300 : 300 + first_n]  # TODO: change at some point

    data[MODIFICATION_STRING] = data_config.standard_modification
    data[CHEMICAL_PATTERN] = data_config.standard_chemical_pattern

    data[CELL_LINE] = data_config.cell_line  # TODO: handle empty case
    data[CELL_LINE_DEPMAP_PROXY] = data[CELL_LINE].map(resolve_depmap_proxy)
    data[CELL_LINE_DEPMAP] = data[CELL_LINE].map(resolve_depmap_id)

    data[PS_PATTERN] = data_config.standard_ps_pattern

    data[CELL_LINE_ORGANISM] = data_config.organism_name
    data[TRANSFECTION_RAW] = data_config.transfection_method

    if version is None:
        version = "generated"
    data.insert(0, f"index_{version}", range(1, len(data) + 1))
    return data


def generate_aso_features(data, cache: AssetCache, n_jobs=1, get_feature_dir_func=None):
    original_columns = set(data.columns)

    # TODO: consider utilizing save features
    logger.info("version is None, not saving features to disk")
    calculator = Calculator(
        data=data, data_version=None, overwrite=True, cpus=n_jobs, cache=cache, get_feature_dir=get_feature_dir_func
    )
    calculator.calculate_all()

    final_data = calculator.data

    # Dynamically determine all new features added by the pipeline
    new_features = list(set(final_data.columns) - original_columns)

    return final_data, new_features


def _apply_standard_metadata(data, config):
    """Set the standard chemistry / dose / cell-line / delivery columns the feature pipeline reads."""
    data[MODIFICATION_STRING] = config.standard_modification
    data[CHEMICAL_PATTERN] = config.standard_chemical_pattern
    data[PS_PATTERN] = config.standard_ps_pattern
    data[CELL_LINE] = config.cell_line
    data[CELL_LINE_DEPMAP_PROXY] = data[CELL_LINE].map(resolve_depmap_proxy)
    data[CELL_LINE_DEPMAP] = data[CELL_LINE].map(resolve_depmap_id)
    data[CELL_LINE_ORGANISM] = config.organism_name
    data[TRANSFECTION_RAW] = config.transfection_method
    data[VOLUME_NM] = config.volume
    data[DENSITY_CELLS_PER_WELL] = config.cell_per_well
    data.insert(0, "index_generated", range(1, len(data) + 1))


def _fill_out_of_range_one_hots(featured, model_features):
    """Add the model's high-position one-hot columns that short candidates lack, as all-zero (an
    out-of-range one-hot position is zero by definition). Errors on any other missing model feature."""
    missing = [f for f in model_features if f not in featured.columns]
    unexpected = [f for f in missing if not f.startswith("ohe_pos")]
    if unexpected:
        raise ValueError(f"generated candidates are missing non-positional model features: {unexpected[:5]}")
    for f in missing:
        featured[f] = 0.0
    if missing:
        logger.info(
            f"Filled {len(missing)} out-of-range one-hot position features with 0 (ASOs shorter than the grid)."
        )


def design_asos(
    target_gene,
    gene_sequence=None,
    *,
    cell_line=None,
    aso_sizes=(20,),
    first_n=None,
    top_n=None,
    model_version=None,
    cds_start=None,
    cds_end=None,
    config=None,
    cache=None,
    n_jobs=1,
    off_targets=False,
    off_target_max_distance=2,
    off_target_exclude_genes=None,
):
    """Design ASOs for a target end-to-end: tile candidate ASOs across the gene, compute their
    features, score them with the bundled efficacy model, and return them ranked best-first.

    target_gene:   canonical gene name (e.g. "DIAPH3").
    gene_sequence: target (pre-)mRNA to tile; if None it is looked up from the genome cache for
                   `target_gene`. If given, it is registered as a custom gene so the feature
                   pipeline can resolve it.
    cell_line:     overrides config.cell_line -- the screen context the ranking is conditioned on.
    aso_sizes:     ASO lengths to tile (default 20-mers).
    first_n:       featurize only the first N candidates (a full tiling can be thousands of ASOs).
    top_n:         return only the best N after scoring.
    model_version: which bundled model to score with (default the package default).
    cds_start/cds_end: optional 0-based, half-open coding span for a custom `gene_sequence`, so
                   5'UTR/CDS/3'UTR and start/stop features are annotated; ignored for a genome lookup.
    off_targets:   if True, also align every ranked candidate to the genome with Bowtie and return a
                   per-hit sequence off-target table (see below). Informational only -- it does not
                   affect the features or the ranking. Needs the Bowtie index (`tauso setup-bowtie`).
    off_target_max_distance: max mismatch distance to report off-target hits for (Bowtie -v).
    off_target_exclude_genes: extra gene symbols to treat as on-target and drop from the off-target
                   table. The candidates' own target gene is always excluded; for a custom
                   `gene_sequence` whose real symbol isn't `target_gene`, pass it here.

    Returns the candidates with their features and the `tauso_score_<version>` column, sorted
    best-first. If `off_targets` is True, returns a tuple `(ranked, offtargets)` where `offtargets`
    is a tidy table with one row per off-target hit -- columns rank, aso_sequence, off_target_gene,
    distance (exact mismatch count), region (exon/intron/...), chrom, start, strand -- so the hits
    for one ASO are `offtargets[offtargets.aso_sequence == seq]`.
    """
    from tauso.inference import DEFAULT_VERSION, load_model, predict, score_column

    version = model_version or DEFAULT_VERSION
    config = config or default_config()
    if cell_line is not None:
        config.cell_line = cell_line
    aso_sizes = _validate_aso_sizes(aso_sizes)
    if resolve_depmap_proxy(config.cell_line) is None:
        logger.warning(
            f"cell line {config.cell_line!r} does not resolve to a DepMap line; expression features will be uninformative."
        )

    if cache is None:
        cache = AssetCache(genome="GRCm39" if config.organism_name == "mouse" else "GRCh38")
    if gene_sequence is not None:
        gene_sequence = _normalize_target_sequence(gene_sequence)
        cache.set_custom_gene(name=target_gene, sequence=gene_sequence, cds_start=cds_start, cds_end=cds_end)
    else:
        gene_data = cache.get_full_gene_data()
        if target_gene not in gene_data:
            raise ValueError(
                f"gene {target_gene!r} not found in the genome cache; pass gene_sequence= for a custom target."
            )
        gene_sequence = gene_data[target_gene].full_mrna
    if max(aso_sizes) > len(gene_sequence):
        raise ValueError(f"max aso_size ({max(aso_sizes)}) exceeds target length ({len(gene_sequence)}).")

    candidates = get_initial_data(gene_sequence, aso_sizes=aso_sizes, canonical_name=target_gene)
    if first_n is not None:
        candidates = candidates.head(first_n).copy()
    _apply_standard_metadata(candidates, config)
    logger.info(
        f"design_asos: target={target_gene!r} cell_line={config.cell_line!r} model={version}; "
        f"scoring {len(candidates)} candidate ASOs through the feature pipeline."
    )

    featured, _ = generate_aso_features(candidates, cache, n_jobs=n_jobs)
    _, model_features = load_model(version)
    _fill_out_of_range_one_hots(featured, model_features)

    col = score_column(version)
    featured[col] = predict(featured, version)
    ranked = featured.sort_values(col, ascending=False, kind="stable").reset_index(drop=True)
    if top_n is not None:
        ranked = ranked.head(top_n)
    logger.info(f"design_asos: ranked {len(ranked)} candidate ASOs by {col}.")

    if off_targets:
        genome = "GRCm39" if config.organism_name == "mouse" else "GRCh38"
        offtargets = _sequence_offtarget_table(
            ranked, genome=genome, max_distance=off_target_max_distance, exclude_genes=off_target_exclude_genes
        )
        logger.info(
            f"design_asos: {len(offtargets)} sequence off-target hits (<= {off_target_max_distance} mismatches) "
            f"across {len(ranked)} candidates."
        )
        return ranked, offtargets
    return ranked


_DNA_ALPHABET = set("ACGT")

# ASO-length bounds = the range observed in the training data (shortest 12-mer, longest 28-mer);
# the model is not calibrated outside it, so design_asos rejects out-of-range lengths.
MIN_ASO_LENGTH = 12
MAX_ASO_LENGTH = 28


def _validate_aso_sizes(aso_sizes):
    sizes = list(aso_sizes)
    if not sizes or any(not isinstance(s, int) for s in sizes):
        raise ValueError(f"aso_sizes must be a non-empty list of integers, got {aso_sizes!r}")
    out_of_range = [s for s in sizes if not (MIN_ASO_LENGTH <= s <= MAX_ASO_LENGTH)]
    if out_of_range:
        raise ValueError(
            f"aso_sizes must be within [{MIN_ASO_LENGTH}, {MAX_ASO_LENGTH}] nt (the model's training range); got {out_of_range}"
        )
    return sizes


def _normalize_target_sequence(seq):
    """Uppercase and map U->T; validate the DNA alphabet. Rejects whitespace rather than stripping it."""
    s = str(seq)
    if any(c.isspace() for c in s):
        raise ValueError("gene_sequence must not contain whitespace")
    s = s.upper().replace("U", "T")
    if not s:
        raise ValueError("gene_sequence is empty")
    bad = sorted(set(s) - _DNA_ALPHABET)
    if bad:
        raise ValueError(f"gene_sequence has non-ACGU characters: {bad[:5]}")
    return s


# --- Result views for consumers -------------------------------------------------------------------

_REGION_COLS = [
    ("5'UTR", "struct_sense_in_5utr"),
    ("CDS", "struct_sense_in_cds"),
    ("3'UTR", "struct_sense_in_3utr"),
    ("intron", "struct_sense_in_intron"),
]

# Representative top-N / cutoff for the off-target summary (the store carries several; show one).
_OFFTARGET_SPECIFIC = "off_target_score_specific_BOLTZ_n100_c1000"
_OFFTARGET_GENERAL = "off_target_score_general_BOLTZ_n100_c1000"
_RRNA_TOTAL = "off_target_single_rRNA_total_c1000"
_RNASEH1_FIT = "rnase_krel_score_R4a_krel_dynamic"
G4HUNTER_LIABILITY = 1.5  # |G4Hunter| at/above ~1.5 indicates a likely G-quadruplex


def _col_or_nan(df, name):
    return df[name] if name in df.columns else pd.Series(np.nan, index=df.index)


def _target_region(ranked):
    """One region label per candidate from the (mutually exclusive) struct_sense_in_* one-hots."""
    region = pd.Series("unknown", index=ranked.index)
    for name, col in _REGION_COLS:
        if col in ranked.columns:
            region = region.mask((region == "unknown") & ranked[col].astype(bool), name)
    return region


def summarize_design(ranked, model_version=None):
    """Trim a `design_asos` result to the consumer-facing columns: rank, gene, ASO sequence, length,
    transcript start + region, and the efficacy ranking score (`tauso_score_<version>`). The score is
    a within-experiment RANK (higher = better predicted knockdown), not a percent-inhibition value."""
    from tauso.inference import DEFAULT_VERSION, score_column

    col = score_column(model_version or DEFAULT_VERSION)
    return pd.DataFrame(
        {
            "rank": range(1, len(ranked) + 1),
            "gene": ranked[CANONICAL_GENE_NAME].to_numpy(),
            "aso_sequence": ranked[ASO_SEQUENCE].to_numpy(),
            "length": ranked[STRUCTURE_SENSE_LENGTH].to_numpy(),
            "target_start": _col_or_nan(ranked, "structure_sense_start").to_numpy(),
            "target_region": _target_region(ranked).to_numpy(),
            col: ranked[col].to_numpy(),
        }
    )


def design_details(ranked):
    """Per-candidate safety / liability detail, computed from the full feature frame so it is surfaced
    whether or not these features were selected into the model: transcriptome & rRNA off-target burden,
    RNase H1 cleavage fit, and toxicity motifs (CpG/immune, G-quadruplex / G-run). Returns raw values
    plus liability flags and an implications note -- flags mark candidates to scrutinise, not confirmed
    toxicity."""

    def g(name):
        return _col_or_nan(ranked, name)

    out = pd.DataFrame(
        {
            "aso_sequence": ranked[ASO_SEQUENCE].to_numpy(),
            "offtarget_transcriptome": g(_OFFTARGET_SPECIFIC).to_numpy(),
            "offtarget_genomewide": g(_OFFTARGET_GENERAL).to_numpy(),
            "offtarget_rrna": g(_RRNA_TOTAL).to_numpy(),
            "rnaseh1_cleavage_fit": g(_RNASEH1_FIT).to_numpy(),
            "tox_cpg_count": g("tox_cpg_count").to_numpy(),
            "tox_tlr9_motif": g("tox_tlr9_gtcgtt").to_numpy(),
            "tox_g4hunter_max": g("tox_g4hunter_max").to_numpy(),
            "tox_grun_count": g("tox_gggg_motif_count").to_numpy(),
        }
    )
    out["flag_immune_cpg"] = out["tox_cpg_count"].fillna(0) > 0
    out["flag_hepatotox_g4_grun"] = (out["tox_g4hunter_max"].abs() >= G4HUNTER_LIABILITY) | (
        out["tox_grun_count"].fillna(0) > 0
    )
    out["flag_binds_rrna"] = out["offtarget_rrna"].fillna(0) > 0
    out["liabilities"] = [
        _liability_note(cpg, g4, rrna)
        for cpg, g4, rrna in zip(out["flag_immune_cpg"], out["flag_hepatotox_g4_grun"], out["flag_binds_rrna"])
    ]
    return out


def _liability_note(immune_cpg, hepatotox, binds_rrna):
    notes = []
    if immune_cpg:
        notes.append("CpG motif (TLR9/immunostimulation)")
    if binds_rrna:
        notes.append("rRNA off-target binding")
    if hepatotox:
        notes.append("G-quadruplex / G-run (hepatotoxicity/aggregation)")
    return "; ".join(notes) if notes else "none flagged"


_OFFTARGET_COLS = ["rank", "aso_sequence", "off_target_gene", "distance", "region", "chrom", "start", "strand"]


def _sequence_offtarget_table(ranked, *, genome, max_distance, exclude_genes):
    """Build the tidy sequence off-target table for a set of ranked candidates (see `design_asos`).

    Aligns each candidate ASO to `genome` with Bowtie up to `max_distance` mismatches and returns one
    row per hit to a gene OTHER than the intended target: rank, aso_sequence, off_target_gene, the
    exact mismatch `distance`, the genomic `region` (exon/intron/...), and the locus (chrom/start/
    strand). The candidates' own target gene plus any `exclude_genes` are dropped.

    Requires the Bowtie index; raises with a `tauso setup-bowtie` instruction if it is missing, rather
    than triggering a multi-GB index build.
    """
    import os

    from tauso.data.data import get_paths
    from tauso.off_target.search import annotate_hits, run_bowtie_search

    paths = get_paths(genome)
    sentinel = os.path.join(os.path.dirname(paths["fasta"]), f"{genome}_bowtie_index", "SUCCESS")
    if not os.path.exists(sentinel):
        raise FileNotFoundError(
            f"Bowtie index for {genome} not found; run `tauso setup-bowtie --genome {genome}` to use "
            "off_targets=True (it needs the sequence off-target index)."
        )

    exclude = set(exclude_genes or []) | set(pd.Series(ranked[CANONICAL_GENE_NAME]).dropna().unique())
    rank_of = {seq: i + 1 for i, seq in enumerate(ranked[ASO_SEQUENCE].tolist())}

    unique_seqs = list(dict.fromkeys(ranked[ASO_SEQUENCE].tolist()))
    all_hits = []
    for seq in unique_seqs:
        hits, _ = run_bowtie_search(seq, genome=genome, max_mismatches=max_distance)
        all_hits.extend(hits)
    annotated = annotate_hits(all_hits, genome=genome)

    if annotated.empty:
        return pd.DataFrame(columns=_OFFTARGET_COLS)
    keep = annotated["gene_name"].notna() & ~annotated["gene_name"].isin(exclude)
    tbl = annotated[keep]
    if tbl.empty:
        return pd.DataFrame(columns=_OFFTARGET_COLS)

    out = pd.DataFrame(
        {
            "rank": tbl["sequence"].map(rank_of).to_numpy(),
            "aso_sequence": tbl["sequence"].to_numpy(),
            "off_target_gene": tbl["gene_name"].to_numpy(),
            "distance": tbl["mismatches"].to_numpy(),
            "region": tbl["region_type"].to_numpy(),
            "chrom": tbl["chrom"].to_numpy(),
            "start": tbl["start"].to_numpy(),
            "strand": tbl["strand"].to_numpy(),
        }
    )
    return out.sort_values(["rank", "distance", "off_target_gene", "start"]).reset_index(drop=True)
