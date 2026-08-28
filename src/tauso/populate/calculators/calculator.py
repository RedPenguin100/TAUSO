import logging
import os

import pandas as pd
import pyarrow.parquet as pq

from ...data.consts import (
    CANONICAL_GENE_NAME,
    CELL_LINE_DEPMAP,
    CHEMICAL_PATTERN,
    DENSITY_CELLS_PER_WELL,
    PS_PATTERN,
    STRUCT_SENSE_IN_3UTR,
    STRUCT_SENSE_IN_5UTR,
    STRUCT_SENSE_IN_CDS,
    STRUCT_SENSE_IN_CDS_NON_EXCLUSIVE,
    STRUCT_SENSE_IN_EXON,
    STRUCT_SENSE_IN_EXON_NON_EXCLUSIVE,
    STRUCT_SENSE_IN_INTRON,
    STRUCT_SENSE_IN_UTR,
    STRUCTURE_SENSE_DIST_TO_CANONICAL_START,
    STRUCTURE_SENSE_DIST_TO_CANONICAL_STOP,
    STRUCTURE_SENSE_DIST_TO_CLOSEST_SPLICE_JUNCTION,
    STRUCTURE_SENSE_DIST_TO_CLOSEST_START,
    STRUCTURE_SENSE_DIST_TO_CLOSEST_STOP,
    STRUCTURE_SENSE_DIST_TO_SPLICE_JUNCTION_EXONIC,
    STRUCTURE_SENSE_DIST_TO_SPLICE_JUNCTION_INTRONIC,
    STRUCTURE_SENSE_HOST_EXON_LOG_LENGTH,
    STRUCTURE_SENSE_HOST_INTRON_LOG_LENGTH,
    STRUCTURE_SENSE_JUNCTION_LOGDIST_CLOSEST,
    STRUCTURE_SENSE_JUNCTION_LOGDIST_EXONIC,
    STRUCTURE_SENSE_JUNCTION_LOGDIST_INTRONIC,
    STRUCTURE_SENSE_JUNCTION_SIGNED_LOGDIST_EXONIC,
    STRUCTURE_SENSE_JUNCTION_SIGNED_LOGDIST_INTRONIC,
    STRUCTURE_SENSE_LENGTH,
    STRUCTURE_SENSE_MRNA_DIST_TO_CANONICAL_STOP,
    STRUCTURE_SENSE_MRNA_DIST_TO_CLOSEST_STOP,
    STRUCTURE_SENSE_SIGNED_DIST_TO_CANONICAL_START,
    STRUCTURE_SENSE_SIGNED_DIST_TO_CANONICAL_STOP,
    STRUCTURE_SENSE_SIGNED_DIST_TO_CLOSEST_START,
    STRUCTURE_SENSE_SIGNED_DIST_TO_CLOSEST_STOP,
    STRUCTURE_SENSE_START,
    STRUCTURE_SENSE_START_FROM_END,
    STRUCTURE_SENSE_START_FROM_END_NORM,
    STRUCTURE_SENSE_START_NORM,
    STRUCTURE_SENSE_TYPE,
    TRANSFECTION_RAW,
    VOLUME_NM,
)
from ...debug import log_dataframe_memory
from ...features.hybridization.off_target import OFF_TARGET_TOP_NS, RISEARCH_SCORE_CUTOFFS
from ...features.interaction_features import internal_fold_gymnosis
from ...timer import Timer
from ...util import dna_to_rna
from ..feature_cache import cache_path_if_present, loose_shard_dir, save_feature_internal
from ..populate_context import (
    EXPRESSION_FEATURE_NAMES,
    populate_special_gene_expression,
    populate_target_expression,
    populate_transfection,
)
from ..populate_fold import (
    DEFAULT_ACCESS_SETTINGS,
    DEFAULT_SETTINGS,
    access_feature_name,
    mfe_feature_name,
    populate_access_features,
    populate_mfe_features,
)
from ..populate_sequence import FEATURE_SPECS, populate_sequence_features
from ..populate_structure import get_populated_df_with_structure_features
from .cache import AssetCache

logger = logging.getLogger(__name__)


class Calculator:
    def __init__(
        self,
        data: pd.DataFrame,
        data_version: str | None = None,
        overwrite=False,
        cpus: int | None = None,
        cache: AssetCache | None = None,
        get_feature_dir=None,
    ):
        self.data = data
        self.cpus = cpus or 32
        self.data_version = data_version
        self.index = f"index_{data_version}" if data_version else None
        self.overwrite = overwrite
        self.get_feature_dir_func = get_feature_dir
        self.cache = cache or AssetCache(genome="GRCh38")  # TODO: generalize for mice as well

        self._genes_u = None
        self._context_added = False

        if self.get_feature_dir_func is None:
            logger.warning(
                "[Calculator] get_feature_dir_func is None. To save features, please pass a function to the calculator."
            )
        logger.info("[Calculator] Initialized successfully.")

    def _save_calculated_feature(self, feature_name):
        # Will silently fail, but the warning will be displayed in the constructor
        if self.get_feature_dir_func is not None:
            save_feature_internal(
                self.data,
                feature_name=feature_name,
                overwrite=self.overwrite,
                version=self.data_version,
                saved_dir_func=self.get_feature_dir_func,
            )

    def _check_dependencies(self, required_columns: list):
        """Helper method to ensure upstream prep steps populated the necessary columns."""
        missing = [col for col in required_columns if col not in self.data.columns]
        if missing:
            raise ValueError(f"Missing required dependencies in dataframe: {missing}")

    def _require_str_columns(self, columns: list) -> None:
        """Validates the given columns exist and hold only strings; logs and raises otherwise."""
        for col in columns:
            if col not in self.data.columns:
                message = f"Missing required dependencies in dataframe: {[col]}"
                logger.error(message)
                raise ValueError(message)
            is_str = self.data[col].map(lambda x: isinstance(x, str))
            if not is_str.all():
                message = f"{col!r} must contain only strings; found {int((~is_str).sum())} non-string value(s)"
                logger.error(message)
                raise TypeError(message)

    def _cache_columns(self):
        """Column names in the locally-present wide cache for this run, or empty set."""
        cache = cache_path_if_present(self.data_version) if self.data_version else None
        if cache is None:
            return set(), None
        return set(pq.read_schema(cache).names), cache

    def _resolve_feature_source(self, feature_dir, feature, cache_cols, cache_path):
        """Return ('loose', path) | ('cache', cache_path) | None for `feature`. Loose wins.

        Loose shards are read from `<feature_dir>/_patches/` first, then from `feature_dir`
        itself (legacy location, migration window).
        """
        for parent in (loose_shard_dir(feature_dir), feature_dir):
            for ext in (".parquet", ".csv"):
                path = os.path.join(parent, f"{feature}{ext}")
                if os.path.exists(path):
                    return ("loose", path)
        if feature in cache_cols:
            return ("cache", cache_path)
        return None

    def _load_features_into_data(self, feature_names: list):
        """Reads saved feature shards (or the wide cache) back into self.data for in-memory dependencies."""
        to_load = [f for f in feature_names if f not in self.data.columns]
        if not to_load:  # dependencies already in memory (e.g. the in-memory generation path)
            return
        if self.index is None:
            raise ValueError(f"cannot load {to_load} from disk without a data_version")
        feature_dir = self.get_feature_dir_func(self.data_version)
        cache_cols, cache = self._cache_columns()
        for feature in to_load:
            src = self._resolve_feature_source(feature_dir, feature, cache_cols, cache)
            if src is None:
                raise FileNotFoundError(
                    f"No shard or cache column for '{feature}' in {feature_dir} (.parquet, .csv, or wide cache)."
                )
            kind, path = src
            if kind == "cache":
                feat_df = pd.read_parquet(path, columns=[self.index, feature])
            elif path.endswith(".parquet"):
                feat_df = pd.read_parquet(path)
            else:
                feat_df = pd.read_csv(path)
            self.data = self.data.merge(feat_df[[self.index, feature]], on=self.index, how="left")
            logger.debug("Loaded '%s' from %s into DataFrame.", feature, kind)

    def _get_missing_features(self, expected_features: list) -> list:
        if self.overwrite:
            return expected_features
        if self.get_feature_dir_func is None:
            # Features are never written without a feature dir, so none can be found on disk.
            return expected_features

        feature_dir = self.get_feature_dir_func(self.data_version)
        cache_cols, cache = self._cache_columns()
        missing = []

        for feature in expected_features:
            if self._resolve_feature_source(feature_dir, feature, cache_cols, cache) is not None:
                logger.debug("Skipping '%s': already present.", feature)
            else:
                missing.append(feature)

        return missing

    def _step(self, label, expected, compute, *, save_only_missing=False, load_if_present=False):
        """Run one feature step: skip when everything is already stored, else compute and save.

        `compute` receives the missing feature names and returns the ``(data, generated_names)``
        pair the ``populate_*`` functions produce. `save_only_missing` restricts saving to names
        that were absent, for steps whose populate function returns its whole family regardless
        of what was asked for. `load_if_present` reads the features back into ``self.data`` when
        the step is skipped, for the columns later steps depend on.
        """
        missing = self._get_missing_features(expected)
        if not missing:
            if load_if_present:
                logger.info("All %s features exist. Loading from disk...", label)
                self._load_features_into_data(expected)
            else:
                logger.info("All %s features exist. Skipping.", label)
            return []

        logger.info("Computing %d %s features...", len(missing), label)
        self.data, generated = compute(missing)
        for feature in generated:
            if not save_only_missing or feature in missing:
                self._save_calculated_feature(feature_name=feature)
        return generated

    def _get_unique_genes(self):
        """Lazy getter for the unique genes list."""
        if self._genes_u is None:
            self._genes_u = list(set(self.data[CANONICAL_GENE_NAME]))

            if "RNASEH1" not in self._genes_u:
                self._genes_u.append("RNASEH1")
        return self._genes_u

    def _ensure_genomic_context(self, cds_windows: list):
        """Ensures context columns exist in self.data before running codon algorithms."""
        if not self._context_added:
            logger.info("Adding external mRNA and genomic context columns...")
            from tauso.algorithms.genomic_context_windows import add_external_mrna_and_context_columns

            registry = self.cache.get_gene_registry(self._get_unique_genes())

            flank_sizes_premrna = [5, 20, 30, 40, 50, 60, 70]  # 5 = RBP footprint window

            self.data = add_external_mrna_and_context_columns(
                df=self.data,
                gene_registry=registry,
                flank_sizes_premrna=flank_sizes_premrna,
                flank_sizes_cds=cds_windows,
            )
            self._context_added = True

    def calculate_sequence(self):
        self._step(
            "sequence",
            [name for name, _ in FEATURE_SPECS],
            lambda missing: populate_sequence_features(self.data, features=missing, cpus=self.cpus),
        )

    def calculate_basic(self):
        """Calculates fast boolean/categorical and transfection features with dependency checking."""

        # chem_1st_gen completes the {1st, 2nd, 3rd}-generation triple: 1st-gen is a pure
        # DNA/PS oligo (no sugar mods); 2nd-gen has 2'-MOE; 3rd-gen has LNA or cEt. They're
        # not mutually exclusive in general but are in this dataset's gapmer chemistries.
        def compute_chemistry(missing):
            from tauso.data.consts import MODIFICATION_STRING

            self._check_dependencies([MODIFICATION_STRING])
            has_moe = self.data[MODIFICATION_STRING].str.contains("MOE", na=False)
            has_high_aff = self.data[MODIFICATION_STRING].str.contains("LNA|cEt", na=False)

            if "chem_1st_gen" in missing:
                self.data["chem_1st_gen"] = (~(has_moe | has_high_aff)).astype(int)
            if "chem_2nd_gen" in missing:
                self.data["chem_2nd_gen"] = has_moe.astype(int)
            if "chem_3rd_gen" in missing:
                self.data["chem_3rd_gen"] = has_high_aff.astype(int)
            return self.data, missing

        def compute_transfection(missing):
            self._check_dependencies([TRANSFECTION_RAW])
            data, _ = populate_transfection(self.data)
            return data, missing

        self._step("basic chemistry", ["chem_1st_gen", "chem_2nd_gen", "chem_3rd_gen"], compute_chemistry)
        self._step(
            "transfection",
            ["transfection_electroporation", "transfection_gymnosis", "transfection_lipofection"],
            compute_transfection,
        )

    def calculate_structure(self):
        expected_features = [
            STRUCTURE_SENSE_START,
            STRUCTURE_SENSE_START_FROM_END,
            STRUCTURE_SENSE_LENGTH,
            STRUCT_SENSE_IN_EXON,
            STRUCT_SENSE_IN_EXON_NON_EXCLUSIVE,
            STRUCT_SENSE_IN_INTRON,
            STRUCT_SENSE_IN_UTR,
            STRUCT_SENSE_IN_3UTR,
            STRUCT_SENSE_IN_5UTR,
            STRUCT_SENSE_IN_CDS,
            STRUCT_SENSE_IN_CDS_NON_EXCLUSIVE,
            STRUCTURE_SENSE_TYPE,
            STRUCTURE_SENSE_START_NORM,
            STRUCTURE_SENSE_START_FROM_END_NORM,
            STRUCTURE_SENSE_DIST_TO_CANONICAL_STOP,
            STRUCTURE_SENSE_DIST_TO_CLOSEST_STOP,
            STRUCTURE_SENSE_DIST_TO_CANONICAL_START,
            STRUCTURE_SENSE_DIST_TO_CLOSEST_START,
            STRUCTURE_SENSE_SIGNED_DIST_TO_CANONICAL_STOP,
            STRUCTURE_SENSE_SIGNED_DIST_TO_CLOSEST_STOP,
            STRUCTURE_SENSE_SIGNED_DIST_TO_CANONICAL_START,
            STRUCTURE_SENSE_SIGNED_DIST_TO_CLOSEST_START,
            STRUCTURE_SENSE_MRNA_DIST_TO_CANONICAL_STOP,
            STRUCTURE_SENSE_MRNA_DIST_TO_CLOSEST_STOP,
            STRUCTURE_SENSE_DIST_TO_SPLICE_JUNCTION_EXONIC,
            STRUCTURE_SENSE_DIST_TO_SPLICE_JUNCTION_INTRONIC,
            STRUCTURE_SENSE_DIST_TO_CLOSEST_SPLICE_JUNCTION,
            STRUCTURE_SENSE_JUNCTION_LOGDIST_EXONIC,
            STRUCTURE_SENSE_JUNCTION_LOGDIST_INTRONIC,
            STRUCTURE_SENSE_JUNCTION_LOGDIST_CLOSEST,
            STRUCTURE_SENSE_JUNCTION_SIGNED_LOGDIST_EXONIC,
            STRUCTURE_SENSE_JUNCTION_SIGNED_LOGDIST_INTRONIC,
            STRUCTURE_SENSE_HOST_EXON_LOG_LENGTH,
            STRUCTURE_SENSE_HOST_INTRON_LOG_LENGTH,
        ]

        def compute(missing):
            genes_u = self._get_unique_genes()
            gene_to_data = self.cache.get_lean_gene(genes_u=genes_u)
            data = get_populated_df_with_structure_features(self.data, genes_u, gene_to_data)
            # The populate function fills the whole family in one pass.
            return data, expected_features

        # Later steps read these columns, so they are reloaded when the step is skipped.
        self._step("structure", expected_features, compute, load_if_present=True)

    def calculate_expression(self):
        """Calculates mRNA expression features."""

        def compute(missing):
            self._check_dependencies([CELL_LINE_DEPMAP])
            cell_lines_depmap = self.data[CELL_LINE_DEPMAP].dropna().unique().tolist()
            transcriptomes = self.cache.get_transcriptomes(cell_lines_depmap=cell_lines_depmap)
            data, target_feats = populate_target_expression(self.data, transcriptomes)
            data, special_feats = populate_special_gene_expression(data, transcriptomes)
            return data, target_feats + special_feats

        self._step("expression", EXPRESSION_FEATURE_NAMES, compute)

    def calculate_rnase(self):
        """Calculates RNase H features."""
        expected_features = [
            "rnase_krel_dinucleotide_score_R4a_krel_dinuc_dynamic",
            "rnase_krel_dinucleotide_score_R4b_krel_dinuc_dynamic",
            "rnase_krel_dinucleotide_score_R7_krel_dinuc_dynamic",
            "rnase_score_dinucleotide_R4a_dinuc_dynamic",
            "rnase_score_dinucleotide_R4b_dinuc_dynamic",
            "rnase_score_dinucleotide_R7_dinuc_dynamic",
            "rnase_krel_score_R4a_krel_dynamic",
            "rnase_krel_score_R4b_krel_dynamic",
            "rnase_krel_score_R7_krel_dynamic",
            "rnase_score_R4a_dynamic",
            "rnase_score_R4b_dynamic",
            "rnase_score_R7_dynamic",
        ]

        def compute(missing):
            from tauso.populate.populate_rnase import populate_rnase_features

            return populate_rnase_features(self.data)

        self._step("RNase H", expected_features, compute)

    def calculate_off_target_single(self):
        """Off-target features for RNase H1 + cytoplasmic rRNA (18S/5.8S/28S/5S) and their sum.

        rRNA captures the RiboGreen/total-RNA assay confound the transcriptome off-target
        features miss (see rrna_targets.py). Needs the rRNA reference: run `tauso setup-rrna`.
        """
        from tauso.features.hybridization.off_target.rrna_targets import RRNA_ACCESSIONS

        cutoffs = list(RISEARCH_SCORE_CUTOFFS)
        rrna_species = list(RRNA_ACCESSIONS)
        targets = ["RNASEH1"] + rrna_species
        rrna_total_features = [f"off_target_single_rRNA_total_c{c}" for c in cutoffs]
        expected_features = [f"off_target_single_{g}_c{c}" for g in targets for c in cutoffs] + rrna_total_features

        def compute(missing):
            from tauso.features.hybridization.off_target.off_target_specific_gene import (
                off_target_single_gene_hybridization,
            )
            from tauso.features.hybridization.off_target.rrna_targets import get_rrna_loci

            # rRNA species are not in GRCh38; inject their RefSeq loci without clobbering real genes.
            gene_to_data_full = self.cache.get_full_gene_data()
            for name, locus in get_rrna_loci().items():
                gene_to_data_full.setdefault(name, locus)

            generated = []
            # All cutoffs for a gene are derived together per ASO batch.
            for target_gene in targets:
                self.data, feature_names = off_target_single_gene_hybridization(
                    self.data, target_gene, gene_to_data_full, cutoffs=cutoffs, n_jobs=self.cpus
                )
                generated.extend(feature_names)

            for cutoff in cutoffs:
                total_col = f"off_target_single_rRNA_total_c{cutoff}"
                species_cols = [f"off_target_single_{sp}_c{cutoff}" for sp in rrna_species]
                self.data[total_col] = self.data[species_cols].sum(axis=1)
                generated.append(total_col)
            return self.data, generated

        self._step("single-gene off-target", expected_features, compute)

    def calculate_on_target_site_features(self):
        """On-target site features against each ASO's canonical gene: total hybridization
        (Sum exp(-E/RT)) and log effective number of sites (target multiplicity), per score cutoff."""
        cutoffs = list(RISEARCH_SCORE_CUTOFFS)
        expected_features = [f"on_target_total_hybridization_{c}" for c in cutoffs] + [
            f"on_target_log_number_of_sites_{c}" for c in cutoffs
        ]

        def compute(missing):
            from tauso.features.hybridization.off_target.off_target_specific_gene import (
                add_on_target_site_features,
            )

            # On-target evaluates only against each row's canonical gene, so the lean dict suffices.
            gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())
            return add_on_target_site_features(self.data, gene_to_data, cutoffs=cutoffs, n_jobs=self.cpus)

        self._step("on-target site", expected_features, compute, save_only_missing=True)

    def calculate_mfe(self):
        """Calculates Minimum Free Energy (MFE) fold features."""

        # Each setting tuple (flank, window, step) under the name of the column it fills.
        settings_by_name = {mfe_feature_name(*setting): setting for setting in DEFAULT_SETTINGS}

        def compute(missing):
            self._check_dependencies([STRUCTURE_SENSE_START, STRUCTURE_SENSE_LENGTH])
            gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())
            # All settings are calculated in a single application.
            settings = [settings_by_name[name] for name in missing]
            return populate_mfe_features(self.data, gene_to_data, n_jobs=self.cpus, verbose=False, settings=settings)

        self._step("MFE", list(settings_by_name), compute)

    def calculate_access(self):
        """Calculates target-site accessibility fold features."""

        # Each setting tuple (flank, max_bp_span, open_len) under the name of the column it fills.
        settings_by_name = {access_feature_name(*setting): setting for setting in DEFAULT_ACCESS_SETTINGS}

        def compute(missing):
            self._check_dependencies([STRUCTURE_SENSE_START, STRUCTURE_SENSE_LENGTH])
            gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())
            # Shared flanks may be computed in a single pass.
            settings = [settings_by_name[name] for name in missing]
            return populate_access_features(self.data, gene_to_data, n_jobs=self.cpus, verbose=False, settings=settings)

        self._step("accessibility", list(settings_by_name), compute)

    def calculate_sequence_one_hot(self):
        """Calculates terminal one-hot encoded sequence features."""

        # The expected names come from the same helper the populate step uses, so the two
        # lists cannot drift apart.
        from tauso.populate.populate_sequence import one_hot_feature_names

        def compute(missing):
            from tauso.populate.populate_sequence import populate_sequence_one_hot_encoded

            return populate_sequence_one_hot_encoded(self.data, cpus=self.cpus)

        self._step("sequence one-hot", one_hot_feature_names(), compute, save_only_missing=True)

    def calculate_sequence_chemistry(self):
        """Calculates sequence chemistry features."""
        from tauso.populate.populate_sequence_chemistry import FEATURE_SPECS as CHEMISTRY_SPECS

        def compute(missing):
            from tauso.populate.populate_sequence_chemistry import populate_sequence_chemistry_features

            return populate_sequence_chemistry_features(self.data, features=missing, cpus=self.cpus)

        self._step("sequence chemistry", [name for name, _ in CHEMISTRY_SPECS], compute)

    def calculate_toxicity(self):
        """Calculates sequence-derived toxicity / liability features (tox_* family)."""
        from tauso.populate.populate_toxicity import FEATURE_SPECS as TOX_SPECS

        def compute(missing):
            from tauso.populate.populate_toxicity import populate_toxicity_features

            return populate_toxicity_features(self.data, features=missing, cpus=self.cpus)

        self._step("toxicity", [name for name, _ in TOX_SPECS], compute)

    def calculate_hybridization(self):
        """Calculates hybridization features."""
        from tauso.populate.populate_hybridization import HYBR_FEATURE_TO_CALCULATION

        def compute(missing):
            from tauso.populate.populate_hybridization import populate_hybridization

            return populate_hybridization(self.data, n_cores=self.cpus, features_to_run=missing)

        self._step("hybridization", list(HYBR_FEATURE_TO_CALCULATION), compute)

    def calculate_modification(self):
        """Calculates modification features."""
        from tauso.populate.populate_modification import MODIFICATION_FEATURE_TO_CALCULATION

        def compute(missing):
            from tauso.populate.populate_modification import populate_modifications

            return populate_modifications(self.data, n_cores=self.cpus, features_to_run=missing)

        self._step("modification", list(MODIFICATION_FEATURE_TO_CALCULATION), compute)

    def calculate_backbone_features(self):
        """Populates the PS-backbone features (mod_ps_* family) from PS_PATTERN."""
        from tauso.populate.populate_backbone import BACKBONE_FEATURES, populate_backbone_features

        self._step(
            "backbone",
            BACKBONE_FEATURES,
            lambda missing: populate_backbone_features(self.data, features=missing),
        )

    def calculate_interaction(self):
        """ASO self-fold (RNA parameters) gated to gymnotic uptake."""
        feature = "interaction_internal_fold_rna_gymnosis"

        def compute(missing):
            deps = ["seq_internal_fold_rna", "transfection_gymnosis"]
            self._load_features_into_data(deps)
            self._check_dependencies(deps)
            self.data[feature] = internal_fold_gymnosis(
                self.data["seq_internal_fold_rna"], self.data["transfection_gymnosis"]
            )
            return self.data, [feature]

        self._step("interaction", [feature], compute)

    def calculate_experimental_conditions(self):
        """Pass-through experimental-condition features: ASO dose and plating density.

        Both arrive on the loaded data (renamed from the raw OligoAI columns) and are stored as
        raw values -- XGBoost is monotone-invariant, so no log transform is applied.
        """

        def compute(missing):
            # Nothing to compute: the columns arrive on the loaded data, so the step only
            # checks they are there and hands them to the saver.
            self._check_dependencies(missing)
            return self.data, missing

        self._step("experimental-condition", [VOLUME_NM, DENSITY_CELLS_PER_WELL], compute)

    def calculate_ribo_seq(self):
        """Calculates ribosome profiling (Ribo-seq) features for both 40S and 80S subunits."""
        flanks = (0, 10, 20, 50, 100, 125, 150)
        how = "mean"
        tracks = ("40s", "80s")

        from ...features.context.ribo_seq import add_genomic_coordinates, feature_names, get_feature_prefix
        from ..populate_context import populate_ribo_seq

        expected_features = []
        for track in tracks:
            expected_features.extend(feature_names(flanks, how, prefix=get_feature_prefix(track)))

        def compute(missing):
            gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())
            self.data = add_genomic_coordinates(self.data, gene_to_data)
            generated = []
            for track in tracks:
                self.data, track_features = populate_ribo_seq(
                    "human", self.data, flanks=flanks, how=how, n_jobs=self.cpus, track=track
                )
                generated.extend(track_features)
            return self.data, generated

        self._step("Ribo-seq", expected_features, compute, save_only_missing=True)

    def calculate_cub(self):
        """Calculates all Codon Usage Bias (CUB) features: tAI, CAI, and ENC."""
        from tauso.populate.populate_codon_usage import populate_cai, populate_enc, populate_tai

        cds_windows = [20, 30, 40, 50, 60, 70]

        def expected(prefix):
            return [f"{prefix}_{flank}" for flank in cds_windows] + [f"{prefix}_global"]

        families = [("tAI", "tai_score"), ("CAI", "cai_score"), ("ENC", "enc_score")]
        if not any(self._get_missing_features(expected(prefix)) for _, prefix in families):
            logger.info("All Codon Usage Bias (CUB) features exist. Skipping.")
            return

        self._check_dependencies([STRUCTURE_SENSE_START])

        # The registry and the genomic context are shared by all three families, so they are
        # loaded once here rather than inside each step.
        registry = self.cache.get_gene_registry(self._get_unique_genes())
        self._ensure_genomic_context(cds_windows)

        # populate_tai takes no n_jobs; the other two do.
        self._step(
            "tAI",
            expected("tai_score"),
            lambda missing: populate_tai(self.data, cds_windows, registry),
            save_only_missing=True,
        )
        self._step(
            "CAI",
            expected("cai_score"),
            lambda missing: populate_cai(self.data, cds_windows, registry, n_jobs=self.cpus),
            save_only_missing=True,
        )
        self._step(
            "ENC",
            expected("enc_score"),
            lambda missing: populate_enc(self.data, cds_windows, registry, n_jobs=self.cpus),
            save_only_missing=True,
        )

    def calculate_off_target_general(self):
        """Calculates general off-target hybridization scores."""
        from tauso.features.hybridization.off_target.add_off_target_feat import AggregationMethod
        from tauso.populate.populate_off_target import serialize_feature_name

        # Define the parameter spaces
        methods = [AggregationMethod.BOLTZMANN_SUM]
        top_ns = list(OFF_TARGET_TOP_NS)
        cutoffs = list(RISEARCH_SCORE_CUTOFFS)

        # Generate all combinations dynamically
        configs = [(m, n, c) for m in methods for n in top_ns for c in cutoffs]

        expected_features = [serialize_feature_name(m, n, c, is_specific=False) for m, n, c in configs]

        def compute(missing):
            from tauso.populate.populate_off_target import populate_off_target_general

            gene_to_data = self.cache.get_full_gene_data()
            self._check_dependencies([CELL_LINE_DEPMAP])
            cell_lines_depmap = self.data[CELL_LINE_DEPMAP].dropna().unique().tolist()
            transcriptomes = self.cache.get_transcriptomes(cell_lines_depmap=cell_lines_depmap)

            # All (top_n, cutoff) features are derived together per chunk:
            # target FASTA built at max(needed_top_ns), smaller top_n derived by
            # gene-subset filter, cutoffs by score-filter on the streaming pyarrow
            # output. Only the top_n values whose features are still missing are
            # passed in, so we don't grow the target FASTA past what's needed.
            generated = []
            for method in methods:
                needed_top_ns = sorted(
                    {
                        n
                        for n in top_ns
                        for c in cutoffs
                        if serialize_feature_name(method, n, c, is_specific=False) in missing
                    }
                )
                if not needed_top_ns:
                    continue

                self.data, method_features = populate_off_target_general(
                    ASO_df=self.data,
                    gene_to_data=gene_to_data,
                    cell_line2data=transcriptomes,
                    top_n_list=needed_top_ns,
                    cutoff_list=cutoffs,
                    method=method,
                    n_jobs=self.cpus,
                )
                generated.extend(method_features)
            return self.data, generated

        self._step("general off-target", expected_features, compute, save_only_missing=True)

    def calculate_off_target_specific(self):
        """Calculates cell-line specific off-target hybridization scores."""
        from tauso.features.hybridization.off_target.add_off_target_feat import AggregationMethod
        from tauso.populate.populate_off_target import serialize_feature_name

        method = AggregationMethod.BOLTZMANN_SUM
        top_n_list = list(OFF_TARGET_TOP_NS)
        cutoff_list = list(RISEARCH_SCORE_CUTOFFS)

        # Generate expected combinations
        expected_features = [
            serialize_feature_name(method, n, c, is_specific=True) for n in top_n_list for c in cutoff_list
        ]

        def compute(missing):
            from tauso.populate.populate_off_target import populate_off_target_specific

            gene_to_data = self.cache.get_full_gene_data()
            self._check_dependencies([CELL_LINE_DEPMAP])
            cell_lines_depmap = self.data[CELL_LINE_DEPMAP].dropna().unique().tolist()
            transcriptomes = self.cache.get_transcriptomes(cell_lines_depmap=cell_lines_depmap)

            # All (top_n, cutoff) features are derived together per
            # (cell_line, chunk): target FASTA built at max(needed_top_ns), smaller
            # top_n derived by gene-subset filter, cutoffs by score-filter on the
            # streaming pyarrow output. Only the top_n values whose features are
            # still missing are passed in.
            needed_top_ns = sorted(
                {
                    n
                    for n in top_n_list
                    for c in cutoff_list
                    if serialize_feature_name(method, n, c, is_specific=True) in missing
                }
            )
            if not needed_top_ns:
                return self.data, []

            return populate_off_target_specific(
                ASO_df=self.data,
                gene_to_data=gene_to_data,
                cell_line2data=transcriptomes,
                top_n_list=needed_top_ns,
                cutoff_list=cutoff_list,
                method=method,
                n_jobs=self.cpus,
            )

        self._step("cell-line specific off-target", expected_features, compute, save_only_missing=True)

    def calculate_mrna_halflife(self):
        """Calculates mRNA stability and half-life features."""
        expected_features = ["halflife_value", "halflife_source", "halflife_cell_proxy"]

        def compute(missing):
            # Dependencies are checked before the provider is loaded, which is expensive.
            from tauso.data.consts import CANONICAL_GENE_NAME, CELL_LINE
            from tauso.populate.populate_mrna_halflife import populate_mrna_halflife_features

            self._check_dependencies([CANONICAL_GENE_NAME, CELL_LINE])
            provider = self.cache.get_halflife_provider()
            return populate_mrna_halflife_features(self.data, provider)

        self._step("mRNA half-life", expected_features, compute, save_only_missing=True)

    def calculate_rbp(self):
        """RBP motif-occupancy features."""
        flank_size = 5
        window_col = f"flank_sequence_{flank_size}"

        # Sentinel: if the summary feature exists on disk the block already completed.
        sentinel = [f"rbp_interaction_total_{flank_size}_generic"]

        def compute(missing):
            from tauso.populate.populate_rbp import (
                populate_complexity_features,
                populate_rbp_affinity_features,
            )

            # RBP uses only the +-5 pre-mRNA flank, which _ensure_genomic_context always builds;
            # it needs no CDS windows, so none are requested.
            self._ensure_genomic_context(cds_windows=[])
            rbp_map, pwm_db = self.cache.get_rbp_assets()
            # Empty gene->data map: every row falls back to the uniform [.25]*4 background, so no
            # per-transcript composition is used.
            uniform_background = {}
            data, ind_feats = populate_rbp_affinity_features(
                self.data, rbp_map, pwm_db, uniform_background, window_col, n_jobs=self.cpus
            )
            data, glob_feats = populate_complexity_features(data, ind_feats, suffix=str(flank_size), type="generic")
            return data, ind_feats + glob_feats

        self._step("RBP", sentinel, compute)

    def calculate_oligowalk(self):
        import shutil

        from tauso.features.competition.oligowalk_utils import populate_oligowalk

        if shutil.which("OligoWalk") is None:
            raise RuntimeError(
                "Missing dependency: 'rnastructure' is not installed or not in PATH. "
                "Please install it via Conda (e.g., conda install -c bioconda rnastructure)."
            )

        genes_u = self._get_unique_genes()
        gene_to_data_lean = self.cache.get_lean_gene(genes_u=genes_u)

        self.data, new_features = populate_oligowalk(self.data, gene_to_data_lean)

        if not self.data["error"].isnull().all():
            unique_errors = self.data["error"].dropna().unique()
            err_str = f"Found {len(unique_errors)} unique error types:"
            logger.warning(err_str)
            for err in unique_errors:
                line = f" - {err}"
                logger.warning(line)
                err_str += "\n" + line
            raise ValueError(err_str)

        for feature in new_features:
            if feature != "error":
                self._save_calculated_feature(feature_name=feature)

    def calculate_flank_features(self):
        """Base composition of the +/-20 nt pre-mRNA flanks around the ASO target site: ``flank_at_skew_20``
        = (A-T)/(A+T) and ``flank_gc_content_20`` = (G+C)/(A+C+G+T). NaN when the target is unlocated."""
        from tauso.features.flank_features import compute_flank_composition

        feats = ["flank_at_skew_20", "flank_gc_content_20"]

        def compute(missing):
            self._check_dependencies([STRUCTURE_SENSE_START, STRUCTURE_SENSE_LENGTH])
            gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())
            skew, gc = compute_flank_composition(
                self.data[STRUCTURE_SENSE_START].to_numpy(),
                self.data[STRUCTURE_SENSE_LENGTH].to_numpy(),
                self.data[CANONICAL_GENE_NAME].to_numpy(),
                gene_to_data,
            )
            self.data["flank_at_skew_20"] = skew
            self.data["flank_gc_content_20"] = gc
            return self.data, feats

        self._step("flank composition", feats, compute)

    def calculate_duplication(self):
        """Distinct near-full-length copies of the ASO target in its own pre-mRNA:
        ``on_target_duplication_exact`` (exact) and ``on_target_duplication_near1`` (<=1 mismatch).
        See :mod:`tauso.features.duplication_features`."""
        from tauso.data.consts import ASO_SEQUENCE
        from tauso.features.duplication_features import compute_duplications

        feats = ["on_target_duplication_exact", "on_target_duplication_near1"]

        def compute(missing):
            self._check_dependencies([STRUCTURE_SENSE_START])
            gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())
            gene_mrna = {
                g: dna_to_rna(gene_to_data[g].full_mrna)
                for g in gene_to_data
                if getattr(gene_to_data[g], "full_mrna", None)
            }
            exact, near1 = compute_duplications(
                self.data[STRUCTURE_SENSE_START].to_numpy(),
                self.data[CANONICAL_GENE_NAME].to_numpy(),
                self.data[ASO_SEQUENCE].astype(str).to_numpy(),
                gene_mrna,
            )
            self.data["on_target_duplication_exact"] = exact
            self.data["on_target_duplication_near1"] = near1
            return self.data, feats

        self._step("duplication", feats, compute)

    def calculate_all(self):
        """Executes the full calculation pipeline and times each step."""
        self._require_str_columns([PS_PATTERN, CHEMICAL_PATTERN])

        # 1. Define the pipeline as a list of functions (no parentheses!)
        pipeline_steps = [
            self.calculate_structure,
            self.calculate_cub,
            self.calculate_basic,
            self.calculate_experimental_conditions,
            self.calculate_sequence,
            self.calculate_expression,
            self.calculate_rnase,
            self.calculate_on_target_site_features,
            self.calculate_mfe,
            self.calculate_access,
            self.calculate_flank_features,
            self.calculate_duplication,
            self.calculate_sequence_one_hot,
            self.calculate_sequence_chemistry,
            self.calculate_toxicity,
            self.calculate_modification,
            self.calculate_hybridization,
            self.calculate_backbone_features,
            self.calculate_interaction,
            self.calculate_ribo_seq,
            self.calculate_off_target_general,
            self.calculate_off_target_single,
            self.calculate_off_target_specific,
            self.calculate_mrna_halflife,
            self.calculate_rbp,
        ]

        logger.info("=== Starting pipeline with %d steps ===", len(pipeline_steps))

        # Opt-in: a long cluster run that dies on memory can be re-run with this set to
        # see which step grew the frame and which columns carry it.
        profile_memory = os.environ.get("TAUSO_PROFILE_MEM") == "1"

        for step in pipeline_steps:
            # step.__name__ dynamically grabs the name of the function (e.g., 'calculate_cub')
            logger.info(f"[Calculator] Starting step: {step.__name__}")
            try:
                with Timer(name=step.__name__):
                    step()
                if profile_memory:
                    log_dataframe_memory(self.data, f"after {step.__name__}")
            except Exception:
                logger.error(f"\n[Calculator] Crashed during step: {step.__name__}\n")
                raise

        logger.info("[Calculator] Pipeline Complete")
