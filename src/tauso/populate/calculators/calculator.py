import logging
import os

import numpy as np
import pandas as pd
import psutil
import pyarrow.parquet as pq
from pympler import asizeof

from ...data.consts import (
    CANONICAL_GENE_NAME,
    CELL_LINE_DEPMAP,
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
    STRUCTURE_SENSE_DIST_TO_CLOSEST_START,
    STRUCTURE_SENSE_DIST_TO_CLOSEST_STOP,
    STRUCTURE_SENSE_DIST_TO_SPLICE_JUNCTION_EXONIC,
    STRUCTURE_SENSE_DIST_TO_SPLICE_JUNCTION_INTRONIC,
    STRUCTURE_SENSE_LENGTH,
    STRUCTURE_SENSE_MRNA_DIST_TO_CANONICAL_STOP,
    STRUCTURE_SENSE_MRNA_DIST_TO_CLOSEST_STOP,
    STRUCTURE_SENSE_START,
    STRUCTURE_SENSE_START_FROM_END,
    STRUCTURE_SENSE_START_FROM_END_NORM,
    STRUCTURE_SENSE_START_NORM,
    STRUCTURE_SENSE_TYPE,
    TRANSFECTION_RAW,
)
from ...features.interaction_features import protein_affinity_gymnosis
from ...timer import Timer
from ..feature_cache import cache_path_if_present, loose_shard_dir, save_feature_internal
from ..populate_context import (
    EXPRESSION_FEATURE_NAMES,
    populate_special_gene_expression,
    populate_target_expression,
    populate_transfection,
)
from ..populate_sequence import FEATURE_SPECS, populate_sequence_features
from ..populate_structure import get_populated_df_with_structure_features
from .cache import AssetCache

logger = logging.getLogger(__name__)


class Calculator:
    def __init__(
        self,
        data: pd.DataFrame,
        data_version: str = None,
        overwrite=False,
        cpus: int = None,
        cache: AssetCache = None,
        get_feature_dir=None,
    ):
        self.data = data
        self.cpus = cpus if cpus else 32
        self.data_version = data_version
        self.index = f"index_{data_version}" if data_version else None
        self.overwrite = overwrite
        self.get_feature_dir_func = get_feature_dir
        self.cache = cache if cache else AssetCache(genome="GRCh38")  # TODO: generalize for mice as well

        self._genes_u = None
        self._genes_u = self._get_unique_genes()
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
        feature_dir = self.get_feature_dir_func(self.data_version)
        cache_cols, cache = self._cache_columns()
        for feature in feature_names:
            if feature in self.data.columns:
                continue
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

        feature_dir = self.get_feature_dir_func(self.data_version)
        cache_cols, cache = self._cache_columns()
        missing = []

        for feature in expected_features:
            if self._resolve_feature_source(feature_dir, feature, cache_cols, cache) is not None:
                logger.debug("Skipping '%s': already present.", feature)
            else:
                missing.append(feature)

        return missing

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
        features = [name for name, _ in FEATURE_SPECS]
        missing = self._get_missing_features(features)

        if missing:
            logger.info("Computing %d sequence features...", len(missing))
            self.data, feature_names = populate_sequence_features(self.data, features=missing, cpus=self.cpus)

            for feature in feature_names:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All sequence features exist. Skipping.")

    def calculate_basic(self):
        """Calculates fast boolean/categorical and transfection features with dependency checking."""

        # ==========================================
        # 1. Basic Chemistry Features (chem_*_gen generation flags)
        # ==========================================
        # chem_1st_gen completes the {1st, 2nd, 3rd}-generation triple: 1st-gen is a pure
        # DNA/PS oligo (no sugar mods); 2nd-gen has 2'-MOE; 3rd-gen has LNA or cEt. They're
        # not mutually exclusive in general but are in this dataset's gapmer chemistries.
        expected_basic = ["chem_1st_gen", "chem_2nd_gen", "chem_3rd_gen"]
        missing_basic = self._get_missing_features(expected_basic)

        if missing_basic:
            logger.info("Computing %d basic features...", len(missing_basic))

            from tauso.data.consts import MODIFICATION_STRING

            self._check_dependencies([MODIFICATION_STRING])

            has_moe = self.data[MODIFICATION_STRING].str.contains("MOE", na=False)
            has_high_aff = self.data[MODIFICATION_STRING].str.contains("LNA|cEt", na=False)

            if "chem_1st_gen" in missing_basic:
                self.data["chem_1st_gen"] = (~(has_moe | has_high_aff)).astype(int)
            if "chem_2nd_gen" in missing_basic:
                self.data["chem_2nd_gen"] = has_moe.astype(int)
            if "chem_3rd_gen" in missing_basic:
                self.data["chem_3rd_gen"] = has_high_aff.astype(int)

            for feature in missing_basic:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All basic chemistry features exist. Skipping.")

        # ==========================================
        # 2. 5'-Terminal Nucleotide Features (term5p_*)
        # ==========================================
        # 5'-terminal base identity. RNase H1 has a documented sequence preference at the
        # cleavage site (Wu & Lima, JBC 2004; Lima et al., JBC 2004) and the 5' base of
        # the ASO drives that preference. Stored as DNA letters so seq[0] == 'U' is also
        # accepted as 'T' in term5p_is_t.
        expected_term5p = ["term5p_is_purine", "term5p_is_g", "term5p_is_t"]
        missing_term5p = self._get_missing_features(expected_term5p)

        if missing_term5p:
            logger.info("Computing %d 5'-terminal base features...", len(missing_term5p))

            from tauso.data.consts import ASO_SEQUENCE

            self._check_dependencies([ASO_SEQUENCE])

            seq_5p = self.data[ASO_SEQUENCE].str[0]
            if "term5p_is_purine" in missing_term5p:
                self.data["term5p_is_purine"] = seq_5p.isin(["A", "G"]).astype(int)
            if "term5p_is_g" in missing_term5p:
                self.data["term5p_is_g"] = (seq_5p == "G").astype(int)
            if "term5p_is_t" in missing_term5p:
                self.data["term5p_is_t"] = seq_5p.isin(["T", "U"]).astype(int)

            for feature in missing_term5p:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All 5'-terminal base features exist. Skipping.")

        # ==========================================
        # 3. Transfection Features
        # ==========================================
        expected_transfection = ["transfection_electroporation", "transfection_gymnosis", "transfection_lipofection"]
        missing_transfection = self._get_missing_features(expected_transfection)

        if missing_transfection:
            logger.info("Computing %d transfection features...", len(missing_transfection))

            # 2. Check Dependency
            self._check_dependencies([TRANSFECTION_RAW])

            self.data, _ = populate_transfection(self.data)

            # 4. Save
            for feature in missing_transfection:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All transfection features exist. Skipping.")

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
            STRUCTURE_SENSE_MRNA_DIST_TO_CANONICAL_STOP,
            STRUCTURE_SENSE_MRNA_DIST_TO_CLOSEST_STOP,
            STRUCTURE_SENSE_DIST_TO_SPLICE_JUNCTION_EXONIC,
            STRUCTURE_SENSE_DIST_TO_SPLICE_JUNCTION_INTRONIC,
        ]

        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing structure features...")

            genes_u = self._get_unique_genes()
            gene_to_data = self.cache.get_lean_gene(genes_u=genes_u)

            self.data = get_populated_df_with_structure_features(self.data, genes_u, gene_to_data)

            for feature in expected_features:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All structure features exist. Loading from disk...")
            self._load_features_into_data(expected_features)

    def calculate_expression(self):
        """Calculates mRNA expression features."""
        expected_features = EXPRESSION_FEATURE_NAMES
        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d expression features...", len(missing))

            self._check_dependencies([CELL_LINE_DEPMAP])
            cell_lines_depmap = self.data[CELL_LINE_DEPMAP].dropna().unique().tolist()

            transcriptomes = self.cache.get_transcriptomes(cell_lines_depmap=cell_lines_depmap)

            self.data, target_feats = populate_target_expression(self.data, transcriptomes)
            for feature in target_feats:
                self._save_calculated_feature(feature_name=feature)

            self.data, special_feats = populate_special_gene_expression(self.data, transcriptomes)
            for feature in special_feats:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All expression features exist. Skipping.")

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

        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d RNase H features...", len(missing))

            # Lazy import to keep initialization fast
            from tauso.populate.populate_rnase import populate_rnase_features

            # Compute the features
            self.data, generated_features = populate_rnase_features(self.data)

            # Save them securely
            for feature in generated_features:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All RNase H features exist. Skipping.")

    def calculate_off_target_single(self):
        """Off-target features for RNase H1 + cytoplasmic rRNA (18S/5.8S/28S/5S) and their sum.

        rRNA captures the RiboGreen/total-RNA assay confound the transcriptome off-target
        features miss (see rrna_targets.py). Needs the rRNA reference: run `tauso setup-rrna`.
        """
        from tauso.features.hybridization.off_target.rrna_targets import RRNA_ACCESSIONS

        cutoffs = [800, 1000, 1200]
        rrna_species = list(RRNA_ACCESSIONS)
        targets = ["RNASEH1"] + rrna_species
        rrna_total_features = [f"off_target_single_rRNA_total_c{c}" for c in cutoffs]
        expected_features = [f"off_target_single_{g}_c{c}" for g in targets for c in cutoffs] + rrna_total_features

        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d specific off-target features...", len(missing))

            from tauso.features.hybridization.off_target.off_target_specific_gene import (
                off_target_single_gene_hybridization,
            )
            from tauso.features.hybridization.off_target.rrna_targets import get_rrna_loci

            # rRNA species are not in GRCh38; inject their RefSeq loci without clobbering real genes.
            gene_to_data_full = self.cache.get_full_gene_data()
            for name, locus in get_rrna_loci().items():
                gene_to_data_full.setdefault(name, locus)

            # All cutoffs for a gene come from ONE RIsearch pass per ASO batch.
            for target_gene in targets:
                self.data, feature_names = off_target_single_gene_hybridization(
                    self.data, target_gene, gene_to_data_full, cutoffs=cutoffs, n_jobs=self.cpus
                )
                for feature_name in feature_names:
                    self._save_calculated_feature(feature_name=feature_name)

            for cutoff in cutoffs:
                total_col = f"off_target_single_rRNA_total_c{cutoff}"
                species_cols = [f"off_target_single_{sp}_c{cutoff}" for sp in rrna_species]
                self.data[total_col] = self.data[species_cols].sum(axis=1)
                self._save_calculated_feature(feature_name=total_col)
        else:
            logger.info("All specific off-target features exist. Skipping.")

    def calculate_on_target_hybridization(self):
        """Calculates on-target total hybridization features."""
        cutoffs = [800, 1000, 1200]
        expected_features = [f"on_target_total_hybridization_{c}" for c in cutoffs]

        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d on-target hybridization features...", len(missing))

            from tauso.features.hybridization.off_target.off_target_specific_gene import on_target_total_hybridization

            # Optimization: We can reuse the lean dictionary because on-target
            # only evaluates against the canonical gene of each row.
            gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())

            # All cutoffs come from ONE RIsearch pass per (gene, ASO batch).
            needed_cutoffs = [c for c in cutoffs if f"on_target_total_hybridization_{c}" in missing]
            if needed_cutoffs:
                self.data, generated_names = on_target_total_hybridization(
                    self.data, gene_to_data, cutoffs=needed_cutoffs, n_jobs=self.cpus
                )
                for feature_name in generated_names:
                    self._save_calculated_feature(feature_name=feature_name)
        else:
            logger.info("All on-target hybridization features exist. Skipping.")

    def calculate_mfe(self):
        """Calculates Minimum Free Energy (MFE) fold features."""
        # Lazy import to keep initialization fast
        from tauso.populate.populate_fold import DEFAULT_SETTINGS, populate_mfe_features

        # Dynamically generate the expected feature names from the settings tuples (flank, window, step)
        expected_features = [f"fold_mfe_win{w}_flank{f}_step{s}" for f, w, s in DEFAULT_SETTINGS]

        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d MFE features...", len(missing))

            self._check_dependencies([STRUCTURE_SENSE_START, STRUCTURE_SENSE_LENGTH])

            # Reuse the lean dictionary
            gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())

            # Compute every missing setting in one pass: populate_mfe_features groups
            # settings that share a (flank, window) so their folds are computed once,
            # and the parallel dispatch is paid a single time instead of per setting.
            missing_settings = [
                (f, w, s) for f, w, s in DEFAULT_SETTINGS if f"fold_mfe_win{w}_flank{f}_step{s}" in missing
            ]

            self.data, generated_features = populate_mfe_features(
                self.data, gene_to_data, n_jobs=self.cpus, verbose=False, settings=missing_settings
            )

            for feature in generated_features:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All MFE features exist. Skipping.")

    def calculate_sense_accessibility(self):
        """Calculates sense accessibility features."""
        from tauso.populate.populate_fold import DEFAULT_SENSE_CONFIGURATION, populate_sense_accessibility_batch

        # Dynamically generate expected feature names based on the configs
        expected_features = []
        for config in DEFAULT_SENSE_CONFIGURATION:
            seeds_str = "-".join(map(str, config["seeds"]))
            expected_features.append(
                f"fold_access_{config['flank']}flank_{config['access']}access_{seeds_str}seed_sizes"
            )

        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing sense accessibility features...")

            # Reuse the lean dictionary securely
            gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())

            for config, expected_name in zip(DEFAULT_SENSE_CONFIGURATION, expected_features):
                if expected_name in missing:
                    c_flank = config["flank"]
                    c_access = config["access"]
                    c_seeds = config["seeds"]

                    logger.info("Running: Flank=%d, Access=%d, Seeds=%s...", c_flank, c_access, c_seeds)

                    self.data, generated_name = populate_sense_accessibility_batch(
                        self.data,
                        gene_to_data,
                        batch_size=1000,
                        flank_size=c_flank,
                        access_size=c_access,
                        seed_sizes=c_seeds,
                        n_jobs=self.cpus,
                    )

                    self._save_calculated_feature(feature_name=generated_name)
                    logger.info("Saved: %s", generated_name)
        else:
            logger.info("All sense accessibility features exist. Skipping.")

    def calculate_sequence_one_hot(self):
        """Calculates one-hot encoded sequence features."""

        # 1. Determine the expected names BEFORE running the heavy function
        # We need the max sequence length from the data to predict the names
        from tauso.data.consts import ASO_SEQUENCE

        max_len = int(self.data[ASO_SEQUENCE].str.len().max())

        expected_features = []
        for pos in range(max_len):
            for nuc in ["A", "C", "G", "T"]:
                expected_features.append(f"ohe_pos{pos}_{nuc}")

        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing sequence one-hot features...")
            from tauso.populate.populate_sequence import populate_sequence_one_hot_encoded

            # Pass max_len so it strictly aligns with our expected features
            self.data, generated_features = populate_sequence_one_hot_encoded(
                self.data, max_len=max_len, cpus=self.cpus
            )
            for feature in generated_features:
                if feature in missing:
                    self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All sequence one-hot features exist. Skipping.")

    def calculate_sequence_chemistry(self):
        """Calculates sequence chemistry features."""
        from tauso.populate.populate_sequence_chemistry import FEATURE_SPECS

        expected_features = [name for name, _ in FEATURE_SPECS]
        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d sequence chemistry features...", len(missing))
            from tauso.populate.populate_sequence_chemistry import populate_sequence_chemistry_features

            # Pass only the missing features to avoid redundant calculations
            self.data, generated_features = populate_sequence_chemistry_features(
                self.data, features=missing, cpus=self.cpus
            )

            for feature in generated_features:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All sequence chemistry features exist. Skipping.")

    def calculate_toxicity(self):
        """Calculates sequence-derived toxicity / liability features (tox_* family)."""
        from tauso.populate.populate_toxicity import FEATURE_SPECS as TOX_SPECS

        expected_features = [name for name, _ in TOX_SPECS]
        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d toxicity features...", len(missing))
            from tauso.populate.populate_toxicity import populate_toxicity_features

            self.data, generated_features = populate_toxicity_features(self.data, features=missing, cpus=self.cpus)

            for feature in generated_features:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All toxicity features exist. Skipping.")

    def calculate_hybridization(self):
        """Calculates hybridization features."""
        from tauso.populate.populate_hybridization import HYBR_FEATURE_TO_CALCULATION

        # Dynamically grab expected features (row-wise + derived)
        expected_features = list(HYBR_FEATURE_TO_CALCULATION.keys())

        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d hybridization features...", len(missing))
            from tauso.populate.populate_hybridization import populate_hybridization

            # Pass `missing` so it ONLY computes the deltas
            self.data, generated_features = populate_hybridization(
                self.data, n_cores=self.cpus, features_to_run=missing
            )

            for feature in generated_features:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All hybridization features exist. Skipping.")

    def calculate_modification(self):
        """Calculates modification features."""
        from tauso.populate.populate_modification import MODIFICATION_FEATURE_TO_CALCULATION

        # Dynamically grab expected features
        expected_features = list(MODIFICATION_FEATURE_TO_CALCULATION.keys())
        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d modification features...", len(missing))
            from tauso.populate.populate_modification import populate_modifications

            # Pass `missing` to `features_to_run` so it ONLY computes the deltas
            self.data, generated_features = populate_modifications(
                self.data, n_cores=self.cpus, features_to_run=missing
            )

            for feature in generated_features:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All modification features exist. Skipping.")

    def calculate_backbone_features(self):
        """Calculates features derived from the PS_PATTERN backbone (mod_ps_* family).

        PS_PATTERN encodes one inter-nucleotide BOND per character ('*' = PS,
        'd' = PO), so it has length len(CHEMICAL_PATTERN) - 1 for a regular
        oligo. Each bond is attributed to the region of its 5' nucleotide,
        matching the convention used in get_dna_rna_dg_region.
        """
        expected_features = [
            "mod_ps_po_percentage",
            "mod_ps_end_score",
            "mod_ps_max_consecutive_po",
            "mod_ps_wing5_count",
            "mod_ps_gap_count",
            "mod_ps_wing3_count",
            "mod_ps_frac_mod",
            "mod_ps_frac_dna",
        ]
        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d backbone features...", len(missing))
            import re

            from tauso.common.modifications import get_longest_dna_gap
            from tauso.data.consts import CHEMICAL_PATTERN, PS_PATTERN

            self._check_dependencies([PS_PATTERN])

            if "mod_ps_max_consecutive_po" in missing:
                self.data["mod_ps_max_consecutive_po"] = self.data[PS_PATTERN].apply(
                    lambda x: max((len(chunk) for chunk in x.split("*")), default=0) if isinstance(x, str) else 0
                )

            if "mod_ps_end_score" in missing:
                self.data["mod_ps_end_score"] = self.data[PS_PATTERN].apply(
                    lambda x: (
                        len(re.match(r"^\**", x).group()) + len(re.search(r"\**$", x).group())
                        if isinstance(x, str) and x
                        else 0
                    )
                )

            if "mod_ps_po_percentage" in missing:
                total_po = self.data[PS_PATTERN].apply(lambda x: x.count("d") if isinstance(x, str) else 0)
                total_ps = self.data[PS_PATTERN].apply(lambda x: x.count("*") if isinstance(x, str) else 0)
                backbone_length = total_po + total_ps
                self.data["mod_ps_po_percentage"] = (total_po / backbone_length.replace(0, pd.NA)).fillna(0)

            need_placement = any(c in missing for c in ("mod_ps_wing5_count", "mod_ps_gap_count", "mod_ps_wing3_count"))
            need_interaction = any(c in missing for c in ("mod_ps_frac_mod", "mod_ps_frac_dna"))

            if need_placement or need_interaction:
                self._check_dependencies([CHEMICAL_PATTERN])

            # Routing flag: a row is a "real" gapmer iff it has both flanks AND a deoxy gap.
            # mod_ps_wing5/3_count and mod_ps_frac_mod are NaN on non-gapmer rows because their
            # value is 0 by construction there, which would leak zero-variance noise into
            # split selection.
            def _is_gapmer_pat(pat):
                if not isinstance(pat, str) or not pat:
                    return False
                start, end, gap_len = get_longest_dna_gap(pat)
                return gap_len > 0 and start > 0 and end < len(pat)

            is_gapmer_series = (
                self.data[CHEMICAL_PATTERN].map(_is_gapmer_pat) if need_placement or need_interaction else None
            )

            if need_placement:

                def _ps_placement(chem_pat, ps_pat):
                    if not isinstance(chem_pat, str) or not isinstance(ps_pat, str):
                        return 0, 0, 0
                    gap_start, gap_end, gap_len = get_longest_dna_gap(chem_pat)
                    if gap_len == 0:
                        return 0, 0, 0
                    return (
                        ps_pat[:gap_start].count("*"),
                        ps_pat[gap_start:gap_end].count("*"),
                        ps_pat[gap_end:].count("*"),
                    )

                triples = self.data.apply(lambda r: _ps_placement(r[CHEMICAL_PATTERN], r[PS_PATTERN]), axis=1)
                if "mod_ps_gap_count" in missing:
                    # Meaningful for all chemistries (= total PS for all-DNA, 0 for all-modified)
                    self.data["mod_ps_gap_count"] = triples.map(lambda t: t[1]).astype(int)
                if "mod_ps_wing5_count" in missing:
                    wing5 = triples.map(lambda t: float(t[0]))
                    self.data["mod_ps_wing5_count"] = wing5.where(is_gapmer_series, np.nan)
                if "mod_ps_wing3_count" in missing:
                    wing3 = triples.map(lambda t: float(t[2]))
                    self.data["mod_ps_wing3_count"] = wing3.where(is_gapmer_series, np.nan)

            if need_interaction:

                def _frac_with_ps(chem_pat, ps_pat, marker_chars):
                    if not isinstance(chem_pat, str) or not isinstance(ps_pat, str):
                        return 0.0
                    # Each PS_PATTERN bond is attributed to the region of its 5' nucleotide
                    n = min(len(chem_pat), len(ps_pat))
                    qualifying, ps = 0, 0
                    for i in range(n):
                        if chem_pat[i] in marker_chars:
                            qualifying += 1
                            if ps_pat[i] == "*":
                                ps += 1
                    return ps / qualifying if qualifying else 0.0

                if "mod_ps_frac_mod" in missing:
                    # Fraction of high-affinity-sugar positions (M/C/L) that carry a PS bond
                    # on their 3' side. NaN on non-gapmer rows (no mods or no wings).
                    raw = self.data.apply(
                        lambda r: _frac_with_ps(r[CHEMICAL_PATTERN], r[PS_PATTERN], {"M", "C", "L"}),
                        axis=1,
                    ).astype(float)
                    self.data["mod_ps_frac_mod"] = raw.where(is_gapmer_series, np.nan)
                if "mod_ps_frac_dna" in missing:
                    # Fraction of deoxy ('d') positions that carry a PS bond on their 3' side.
                    # Meaningful for all chemistries (pure-DNA = total PS / total bonds).
                    self.data["mod_ps_frac_dna"] = self.data.apply(
                        lambda r: _frac_with_ps(r[CHEMICAL_PATTERN], r[PS_PATTERN], {"d"}),
                        axis=1,
                    ).astype(float)

            for feature in missing:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All backbone features exist. Skipping.")

    def calculate_interaction(self):
        """Cross-feature interaction features."""
        feature = "interaction_protein_affinity_gymnosis"
        if not self._get_missing_features([feature]):
            logger.info("%s exists. Skipping.", feature)
            return

        from tauso.data.consts import CHEMICAL_PATTERN, PS_PATTERN

        self._load_features_into_data(["seq_internal_fold"])
        self._check_dependencies([PS_PATTERN, CHEMICAL_PATTERN, "seq_internal_fold"])

        if "transfection_gymnosis" in self.data.columns:
            gymnosis = self.data["transfection_gymnosis"]
        else:
            logger.warning("transfection_gymnosis not in data; %s will be 0 for all rows.", feature)
            gymnosis = None
        self.data[feature] = protein_affinity_gymnosis(
            self.data[PS_PATTERN], self.data[CHEMICAL_PATTERN], self.data["seq_internal_fold"], gymnosis
        )
        self._save_calculated_feature(feature_name=feature)

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
        missing = self._get_missing_features(expected_features)

        if not missing:
            logger.info("All Ribo-seq features exist. Skipping.")
            return

        logger.info("Computing %d Ribo-seq features across %d tracks...", len(missing), len(tracks))
        gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())
        self.data = add_genomic_coordinates(self.data, gene_to_data)

        for track in tracks:
            self.data, generated_features = populate_ribo_seq(
                "human",
                self.data,
                flanks=flanks,
                how=how,
                n_jobs=self.cpus,
                track=track,
            )
            for feature in generated_features:
                if feature in missing:
                    self._save_calculated_feature(feature_name=feature)

    def calculate_cub(self):
        """Calculates all Codon Usage Bias (CUB) features: tAI, CAI, and ENC."""
        cds_windows = [20, 30, 40, 50, 60, 70]

        # 1. Define expected features for all three
        expected_tai = [f"tai_score_{flank}" for flank in cds_windows] + ["tai_score_global"]
        expected_cai = [f"cai_score_{flank}" for flank in cds_windows] + ["cai_score_global"]
        expected_enc = [f"enc_score_{flank}" for flank in cds_windows] + ["enc_score_global"]

        # 2. Check disk for missing files
        missing_tai = self._get_missing_features(expected_tai)
        missing_cai = self._get_missing_features(expected_cai)
        missing_enc = self._get_missing_features(expected_enc)

        # 3. If nothing is missing, skip everything
        if not any([missing_tai, missing_cai, missing_enc]):
            logger.info("All Codon Usage Bias (CUB) features exist. Skipping.")
            return

        logger.info("Computing Codon Usage Bias (CUB) features...")

        self._check_dependencies([STRUCTURE_SENSE_START])

        # 4. Load the shared heavy dependencies EXACTLY ONCE
        registry = self.cache.get_gene_registry(self._get_unique_genes())
        self._ensure_genomic_context(cds_windows)

        # 5. Execute conditionally based on what is missing
        if missing_tai:
            logger.info("Computing %d tAI features...", len(missing_tai))
            from tauso.populate.populate_codon_usage import populate_tai

            self.data, generated_features = populate_tai(self.data, cds_windows, registry)
            for feature in generated_features:
                if feature in missing_tai:
                    self._save_calculated_feature(feature_name=feature)

        if missing_cai:
            logger.info("Computing %d CAI features...", len(missing_cai))
            from tauso.populate.populate_codon_usage import populate_cai

            self.data, generated_features = populate_cai(self.data, cds_windows, registry, n_jobs=self.cpus)
            for feature in generated_features:
                if feature in missing_cai:
                    self._save_calculated_feature(feature_name=feature)

        if missing_enc:
            logger.info("Computing %d ENC features...", len(missing_enc))
            from tauso.populate.populate_codon_usage import populate_enc

            self.data, generated_features = populate_enc(self.data, cds_windows, registry, n_jobs=self.cpus)
            for feature in generated_features:
                if feature in missing_enc:
                    self._save_calculated_feature(feature_name=feature)

    def calculate_off_target_general(self):
        """Calculates general off-target hybridization scores."""
        from tauso.features.hybridization.off_target.add_off_target_feat import AggregationMethod
        from tauso.populate.populate_off_target import serialize_feature_name

        # Define the parameter spaces
        methods = [AggregationMethod.BOLTZMANN_SUM]
        top_ns = [50, 100, 200]
        cutoffs = [800, 1000, 1200]

        # Generate all combinations dynamically
        configs = [(m, n, c) for m in methods for n in top_ns for c in cutoffs]

        expected_features = [serialize_feature_name(m, n, c, is_specific=False) for m, n, c in configs]
        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d general off-target features...", len(missing))
            from tauso.populate.populate_off_target import populate_off_target_general

            # Load the heavy dictionaries (happens instantly if already in memory)
            gene_to_data = self.cache.get_full_gene_data()

            self._check_dependencies([CELL_LINE_DEPMAP])
            cell_lines_depmap = self.data[CELL_LINE_DEPMAP].dropna().unique().tolist()

            transcriptomes = self.cache.get_transcriptomes(cell_lines_depmap=cell_lines_depmap)

            # All (top_n, cutoff) features come from ONE RIsearch pass per chunk:
            # target FASTA built at max(needed_top_ns), smaller top_n derived by
            # gene-subset filter, cutoffs by score-filter on the streaming pyarrow
            # output. Only the top_n values whose features are still missing are
            # passed in, so we don't grow the target FASTA past what's needed.
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

                self.data, generated_features = populate_off_target_general(
                    ASO_df=self.data,
                    gene_to_data=gene_to_data,
                    cell_line2data=transcriptomes,
                    top_n_list=needed_top_ns,
                    cutoff_list=cutoffs,
                    method=method,
                    n_jobs=self.cpus,
                )

                for feature in generated_features:
                    if feature in missing:
                        self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All general off-target features exist. Skipping.")

    def calculate_off_target_specific(self):
        """Calculates cell-line specific off-target hybridization scores."""
        from tauso.features.hybridization.off_target.add_off_target_feat import AggregationMethod
        from tauso.populate.populate_off_target import serialize_feature_name

        method = AggregationMethod.BOLTZMANN_SUM
        top_n_list = [50, 100, 200]
        cutoff_list = [800, 1000, 1200]

        # Generate expected combinations
        expected_features = [
            serialize_feature_name(method, n, c, is_specific=True) for n in top_n_list for c in cutoff_list
        ]

        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d specific off-target features...", len(missing))
            from tauso.populate.populate_off_target import populate_off_target_specific

            gene_to_data = self.cache.get_full_gene_data()
            self._check_dependencies([CELL_LINE_DEPMAP])
            cell_lines_depmap = self.data[CELL_LINE_DEPMAP].dropna().unique().tolist()

            transcriptomes = self.cache.get_transcriptomes(cell_lines_depmap=cell_lines_depmap)

            # All (top_n, cutoff) features come from ONE RIsearch pass per
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
            if needed_top_ns:
                self.data, generated_features = populate_off_target_specific(
                    ASO_df=self.data,
                    gene_to_data=gene_to_data,
                    cell_line2data=transcriptomes,
                    top_n_list=needed_top_ns,
                    cutoff_list=cutoff_list,
                    method=method,
                    n_jobs=self.cpus,
                )

                for feature in generated_features:
                    if feature in missing:
                        self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All specific off-target features exist. Skipping.")

    def calculate_mrna_halflife(self):
        """Calculates mRNA stability and half-life features."""
        expected_features = ["halflife_value", "halflife_source", "halflife_cell_proxy"]
        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d mRNA half-life features...", len(missing))

            # Check dependencies BEFORE loading heavy objects
            from tauso.data.consts import CANONICAL_GENE_NAME, CELL_LINE

            self._check_dependencies([CANONICAL_GENE_NAME, CELL_LINE])

            from tauso.populate.populate_mrna_halflife import populate_mrna_halflife_features

            # 1. Lazy load the provider
            provider = self.cache.get_halflife_provider()

            # 2. Compute
            self.data, generated_features = populate_mrna_halflife_features(self.data, provider)

            # 3. Save only what was missing
            for feature in generated_features:
                if feature in missing:
                    self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All mRNA half-life features exist. Skipping.")

    def calculate_rbp(self):
        """RBP motif-occupancy features."""
        flank_size = 5
        window_col = f"flank_sequence_{flank_size}"

        # Sentinel: if the summary feature exists on disk the block already completed.
        expected = [f"rbp_interaction_total_{flank_size}_generic"]
        if not self._get_missing_features(expected):
            logger.info("All RBP features exist. Skipping.")
            return

        logger.info("Computing RBP affinity pipeline...")
        # RBP uses only the +-5 pre-mRNA flank, which _ensure_genomic_context always builds; it
        # needs no CDS windows, so none are requested.
        self._ensure_genomic_context(cds_windows=[])
        rbp_map, pwm_db = self.cache.get_rbp_assets()

        from tauso.populate.populate_rbp import (
            populate_complexity_features,
            populate_rbp_affinity_features,
        )

        # Empty gene->data map: every row falls back to the uniform [.25]*4 background, so no
        # per-transcript composition is used.
        uniform_background = {}
        self.data, ind_feats = populate_rbp_affinity_features(
            self.data, rbp_map, pwm_db, uniform_background, window_col, n_jobs=self.cpus
        )
        self.data, glob_feats = populate_complexity_features(
            self.data, ind_feats, suffix=str(flank_size), type="generic"
        )
        for feature in ind_feats + glob_feats:
            self._save_calculated_feature(feature_name=feature)

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

    def calculate_all(self):
        """Executes the full calculation pipeline and times each step."""

        # 1. Define the pipeline as a list of functions (no parentheses!)
        pipeline_steps = [
            self.calculate_structure,
            self.calculate_cub,
            self.calculate_basic,
            self.calculate_sequence,
            self.calculate_expression,
            self.calculate_rnase,
            self.calculate_on_target_hybridization,
            self.calculate_mfe,
            self.calculate_sense_accessibility,
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

        for step in pipeline_steps:
            # step.__name__ dynamically grabs the name of the function (e.g., 'calculate_cub')
            logger.info(f"[Calculator] Starting step: {step.__name__}")
            try:
                with Timer(name=step.__name__):
                    step()
                    self.monitor_internal_objects(f"AFTER {step.__name__}")
            except Exception as _:
                logger.error(f"\n[Calculator] Crashed during step: {step.__name__}\n")
                raise

        logger.info("[Calculator] Pipeline Complete")

    def monitor_internal_objects(self, label=""):
        logger.debug(f"\n>>> INTERNAL MEMORY REPORT: {label} <<<")

        # 1. Inspect self.data (The DataFrame)
        if hasattr(self, "data") and hasattr(self.data, "memory_usage"):
            usage = self.data.memory_usage(deep=True).sum() / (1024**2)
            logger.debug(f"  [DF] self.data: {usage:.2f} MB")
            # Print top 3 heaviest columns
            col_usage = self.data.memory_usage(deep=True) / (1024**2)
            logger.debug(f"       Heaviest Cols: {col_usage.sort_values(ascending=False).head(3).to_dict()}")

        # 2. Inspect self.cache (The AssetCache)
        if hasattr(self, "cache"):
            # asizeof is recursive; it's better for complex objects like caches
            cache_size = asizeof.asizeof(self.cache) / (1024**2)
            logger.debug(f"  [CACHE] self.cache: {cache_size:.2f} MB")

        # 3. Check for large local variables in the current scope
        logger.debug(f"  [PROC] Total Process RSS: {psutil.Process(os.getpid()).memory_info().rss / (1024**2):.1f} MB")
