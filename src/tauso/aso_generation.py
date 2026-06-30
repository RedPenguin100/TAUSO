import logging

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
    config=None,
    cache=None,
    n_jobs=1,
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

    Returns the candidates with their features and the `tauso_score_<version>` column, sorted best-first.
    """
    from tauso.inference import DEFAULT_VERSION, load_model, predict, score_column

    version = model_version or DEFAULT_VERSION
    config = config or default_config()
    if cell_line is not None:
        config.cell_line = cell_line
    if cache is None:
        cache = AssetCache(genome="GRCm39" if config.organism_name == "mouse" else "GRCh38")
    if gene_sequence is not None:
        cache.set_custom_gene(name=target_gene, sequence=gene_sequence)
    else:
        gene_sequence = cache.get_full_gene_data()[target_gene].full_mrna

    candidates = get_initial_data(gene_sequence, aso_sizes=list(aso_sizes), canonical_name=target_gene)
    if first_n is not None:
        candidates = candidates.head(first_n).copy()
    _apply_standard_metadata(candidates, config)

    featured, _ = generate_aso_features(candidates, cache, n_jobs=n_jobs)
    _, model_features = load_model(version)
    _fill_out_of_range_one_hots(featured, model_features)

    col = score_column(version)
    featured[col] = predict(featured, version)
    ranked = featured.sort_values(col, ascending=False, kind="stable").reset_index(drop=True)
    return ranked.head(top_n) if top_n is not None else ranked
