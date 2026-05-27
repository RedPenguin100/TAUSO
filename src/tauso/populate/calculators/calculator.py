import logging
import os

import pandas as pd
import psutil
from pympler import asizeof

from ...data.consts import (
    CANONICAL_GENE,
    CELL_LINE_DEPMAP,
    SENSE_3UTR,
    SENSE_5UTR,
    SENSE_EXON,
    SENSE_INTRON,
    SENSE_LENGTH,
    SENSE_START,
    SENSE_START_FROM_END,
    SENSE_TYPE,
    SENSE_UTR,
)
from ...timer import Timer
from ..feature_cache import save_feature_internal
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

    def _load_features_into_data(self, feature_names: list):
        """Reads saved feature CSVs back into self.data for columns needed as in-memory dependencies."""
        feature_dir = self.get_feature_dir_func(self.data_version)
        for feature in feature_names:
            if feature in self.data.columns:
                continue
            file_path = os.path.join(feature_dir, f"{feature}.csv")
            feat_df = pd.read_csv(file_path)
            self.data = self.data.merge(feat_df[[self.index, feature]], on=self.index, how="left")
            logger.debug("Loaded '%s' from disk into DataFrame.", feature)

    def _get_missing_features(self, expected_features: list) -> list:
        if self.overwrite:
            return expected_features

        feature_dir = self.get_feature_dir_func(self.data_version)
        missing = []

        for feature in expected_features:
            file_path = os.path.join(feature_dir, f"{feature}.csv")
            if os.path.exists(file_path):
                logger.debug("Skipping '%s': CSV already exists.", feature)
            else:
                missing.append(feature)

        return missing

    def _get_unique_genes(self):
        """Lazy getter for the unique genes list."""
        if self._genes_u is None:
            self._genes_u = list(set(self.data[CANONICAL_GENE]))

            if "RNASEH1" not in self._genes_u:
                self._genes_u.append("RNASEH1")
        return self._genes_u

    def _ensure_genomic_context(self, cds_windows: list):
        """Ensures context columns exist in self.data before running codon algorithms."""
        if not self._context_added:
            logger.info("Adding external mRNA and genomic context columns...")
            from tauso.algorithms.genomic_context_windows import add_external_mrna_and_context_columns

            mapper = self.cache.get_gene_mapper()
            registry = self.cache.get_gene_registry(self._get_unique_genes())

            flank_sizes_premrna = [20, 30, 40, 50, 60, 70]  # From your snippet

            self.data = add_external_mrna_and_context_columns(
                df=self.data,
                mapper=mapper,
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
        # 1. Basic Chemistry & Cell Line Features
        # ==========================================
        expected_basic = ["is_hepa", "chem_2nd_gen", "chem_3rd_gen"]
        missing_basic = self._get_missing_features(expected_basic)

        if missing_basic:
            logger.info("Computing %d basic features...", len(missing_basic))

            # Check dependencies BEFORE attempting to compute
            from tauso.data.consts import CELL_LINE_DEPMAP_PROXY, HEPA_PROXIES, MODIFICATION

            self._check_dependencies([CELL_LINE_DEPMAP_PROXY, MODIFICATION])

            if "is_hepa" in missing_basic:
                self.data["is_hepa"] = self.data[CELL_LINE_DEPMAP_PROXY].isin(HEPA_PROXIES).astype(int)

            if "chem_2nd_gen" in missing_basic:
                self.data["chem_2nd_gen"] = self.data[MODIFICATION].str.contains("MOE", na=False).astype(int)

            if "chem_3rd_gen" in missing_basic:
                self.data["chem_3rd_gen"] = self.data[MODIFICATION].str.contains("LNA|cEt", na=False).astype(int)

            for feature in missing_basic:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All basic chemistry/hepa features exist. Skipping.")

        # ==========================================
        # 2. Transfection Features
        # ==========================================
        expected_transfection = ["Electroporation", "Gymnosis", "Lipofection", "Other"]
        missing_transfection = self._get_missing_features(expected_transfection)

        if missing_transfection:
            logger.info("Computing %d transfection features...", len(missing_transfection))

            # 2. Check Dependency
            self._check_dependencies(["transfection_method"])

            self.data, _ = populate_transfection(self.data)

            # 4. Save
            for feature in missing_transfection:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All transfection features exist. Skipping.")

    def calculate_structure(self):
        expected_features = [
            SENSE_START,
            SENSE_START_FROM_END,
            SENSE_LENGTH,
            SENSE_EXON,
            SENSE_INTRON,
            SENSE_UTR,
            SENSE_3UTR,
            SENSE_5UTR,
            SENSE_TYPE,
            f"n_{SENSE_EXON}",
            f"n_{SENSE_INTRON}",
            f"n_{SENSE_3UTR}",
            f"n_{SENSE_5UTR}",
            f"n_{SENSE_UTR}",
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
            "RNaseH1_Krel_dinucleotide_score_R4a_krel_dinuc_dynamic",
            "RNaseH1_Krel_dinucleotide_score_R4b_krel_dinuc_dynamic",
            "RNaseH1_Krel_dinucleotide_score_R7_krel_dinuc_dynamic",
            "RNaseH1_score_dinucleotide_R4a_dinuc_dynamic",
            "RNaseH1_score_dinucleotide_R4b_dinuc_dynamic",
            "RNaseH1_score_dinucleotide_R7_dinuc_dynamic",
            "RNaseH1_Krel_score_R4a_krel_dynamic",
            "RNaseH1_Krel_score_R4b_krel_dynamic",
            "RNaseH1_Krel_score_R7_krel_dynamic",
            "RNaseH1_score_R4a_dynamic",
            "RNaseH1_score_R4b_dynamic",
            "RNaseH1_score_R7_dynamic",
            "RNaseH1_Potency_TCCC",
            "RNaseH1_Potency_TTCC",
            "RNaseH1_Potency_TCTC",
            "RNaseH1_Inefficacy_GGGG",
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

            self._check_dependencies([SENSE_START, SENSE_LENGTH])

            # Reuse the lean dictionary
            gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())

            # Compute every missing setting in one pass: populate_mfe_features groups
            # settings that share a (flank, window) so their folds are computed once,
            # and the parallel dispatch is paid a single time instead of per setting.
            missing_settings = [(f, w, s) for f, w, s in DEFAULT_SETTINGS if f"fold_mfe_win{w}_flank{f}_step{s}" in missing]

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
            expected_features.append(f"fold_access_{config['flank']}flank_{config['access']}access_{seeds_str}seed_sizes")

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
        from tauso.data.consts import SEQUENCE

        max_len = int(self.data[SEQUENCE].str.len().max())

        expected_features = []
        for pos in range(max_len):
            for nuc in ["A", "C", "G", "T"]:
                expected_features.append(f"OHE_pos{pos}_{nuc}")

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
        """Calculates features derived from the PS_PATTERN backbone."""
        expected_features = ["is_all_ps", "po_percentage", "ps_end_score", "max_consecutive_PO"]
        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d backbone features...", len(missing))
            import re

            from tauso.data.consts import PS_PATTERN

            # Ensure assign_chemistry/assign_backbone has populated PS_PATTERN
            self._check_dependencies([PS_PATTERN])

            if "is_all_ps" in missing:
                self.data["is_all_ps"] = self.data[PS_PATTERN].apply(
                    lambda x: 1 if isinstance(x, str) and x and x.replace("*", "") == "" else 0
                )

            if "max_consecutive_PO" in missing:
                self.data["max_consecutive_PO"] = self.data[PS_PATTERN].apply(
                    lambda x: max((len(chunk) for chunk in x.split("*")), default=0) if isinstance(x, str) else 0
                )

            if "ps_end_score" in missing:
                self.data["ps_end_score"] = self.data[PS_PATTERN].apply(
                    lambda x: (
                        len(re.match(r"^\**", x).group()) + len(re.search(r"\**$", x).group())
                        if isinstance(x, str) and x
                        else 0
                    )
                )

            if "po_percentage" in missing:
                # Compute intermediates securely without cluttering self.data
                total_po = self.data[PS_PATTERN].apply(lambda x: x.count("d") if isinstance(x, str) else 0)
                total_ps = self.data[PS_PATTERN].apply(lambda x: x.count("*") if isinstance(x, str) else 0)
                backbone_length = total_po + total_ps

                # Safe division, avoids ZeroDivisionError
                import pandas as pd

                self.data["po_percentage"] = (total_po / backbone_length.replace(0, pd.NA)).fillna(0)

            for feature in missing:
                self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All backbone features exist. Skipping.")

    def calculate_ribo_seq(self):
        """Calculates ribosome profiling (Ribo-seq) features."""
        flanks = (0, 10, 20, 50, 100, 125, 150)
        how = "mean"

        from ...features.context.ribo_seq import add_genomic_coordinates, feature_names
        from ..populate_context import populate_ribo_seq

        expected_features = feature_names(flanks, how)
        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d Ribo-seq features...", len(missing))
            mapper = self.cache.get_gene_mapper()
            self.data = add_genomic_coordinates(self.data, mapper)
            self.data, generated_features = populate_ribo_seq(
                "human", self.data, flanks=flanks, how=how, n_jobs=self.cpus
            )
            for feature in generated_features:
                if feature in missing:
                    self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All Ribo-seq features exist. Skipping.")

    def calculate_cub(self):
        """Calculates all Codon Usage Bias (CUB) features: tAI, CAI, and ENC."""
        cds_windows = [20, 30, 40, 50, 60, 70]

        # 1. Define expected features for all three
        expected_tai = [f"tAI_score_{flank}_CDS" for flank in cds_windows] + ["tAI_score_global_CDS"]
        expected_cai = [f"CAI_score_{flank}_CDS" for flank in cds_windows] + ["CAI_score_global_CDS"]
        expected_enc = [f"ENC_score_{flank}_CDS" for flank in cds_windows] + ["ENC_score_global_CDS"]

        # 2. Check disk for missing files
        missing_tai = self._get_missing_features(expected_tai)
        missing_cai = self._get_missing_features(expected_cai)
        missing_enc = self._get_missing_features(expected_enc)

        # 3. If nothing is missing, skip everything
        if not any([missing_tai, missing_cai, missing_enc]):
            logger.info("All Codon Usage Bias (CUB) features exist. Skipping.")
            return

        logger.info("Computing Codon Usage Bias (CUB) features...")

        self._check_dependencies([SENSE_START])

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
        methods = [AggregationMethod.ARTM]
        top_ns = [25, 50, 100, 200]
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

            # All cutoffs for a given (method, top_n) come from ONE RIsearch pass per
            # chunk, so group the work by (method, top_n) and derive every missing
            # cutoff together instead of re-running RIsearch per cutoff.
            for method in methods:
                for top_n in top_ns:
                    needed_cutoffs = [
                        c for c in cutoffs if serialize_feature_name(method, top_n, c, is_specific=False) in missing
                    ]
                    if not needed_cutoffs:
                        continue

                    self.data, generated_features = populate_off_target_general(
                        ASO_df=self.data,
                        gene_to_data=gene_to_data,
                        cell_line2data=transcriptomes,
                        top_n_list=[top_n],
                        cutoff_list=needed_cutoffs,
                        method=method,
                        n_jobs=self.cpus,
                    )

                    for feature in generated_features:
                        self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All general off-target features exist. Skipping.")

    def calculate_off_target_specific(self):
        """Calculates cell-line specific off-target hybridization scores."""
        from tauso.features.hybridization.off_target.add_off_target_feat import AggregationMethod
        from tauso.populate.populate_off_target import serialize_feature_name

        method = AggregationMethod.ARTM
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

            # All cutoffs for a given top_n come from ONE RIsearch pass per (cell_line,
            # chunk, shard), so derive every missing cutoff together per top_n.
            for top_n in top_n_list:
                needed_cutoffs = [
                    c for c in cutoff_list if serialize_feature_name(method, top_n, c, is_specific=True) in missing
                ]
                if not needed_cutoffs:
                    continue

                self.data, generated_features = populate_off_target_specific(
                    ASO_df=self.data,
                    gene_to_data=gene_to_data,
                    cell_line2data=transcriptomes,
                    top_n_list=[top_n],
                    cutoff_list=needed_cutoffs,
                    method=method,
                    n_jobs=self.cpus,
                )

                for feature in generated_features:
                    self._save_calculated_feature(feature_name=feature)
        else:
            logger.info("All specific off-target features exist. Skipping.")

    def calculate_mrna_halflife(self):
        """Calculates mRNA stability and half-life features."""
        expected_features = ["mRNA_HalfLife", "HalfLife_Source", "Mapped_Cell_Proxy"]
        missing = self._get_missing_features(expected_features)

        if missing:
            logger.info("Computing %d mRNA half-life features...", len(missing))

            # Check dependencies BEFORE loading heavy objects
            from tauso.data.consts import CANONICAL_GENE, CELL_LINE

            self._check_dependencies([CANONICAL_GENE, CELL_LINE])

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
        """Calculates all RBP Affinity, Interaction, Functional, and Regional features."""
        flank_size = 50
        window_col = f"flank_sequence_{flank_size}"

        # 1. Define "Sentinel" features.
        # If these summary features exist, we know the block successfully completed.
        expected_interaction = [f"RBP_interaction_total_{flank_size}_expression"]
        expected_affinity = [f"RBP_interaction_total_{flank_size}_generic"]
        expected_functional = [f"RBP_interaction_stabilizer_{flank_size}"]
        expected_regions = [
            f"RBP_interaction_total_left",
            f"RBP_interaction_total_core",
            f"RBP_interaction_total_right",
        ]

        # Check disk
        missing_int = self._get_missing_features(expected_interaction)
        missing_aff = self._get_missing_features(expected_affinity)
        missing_func = self._get_missing_features(expected_functional)
        missing_reg = self._get_missing_features(expected_regions)

        if not any([missing_int, missing_aff, missing_func, missing_reg]):
            logger.info("All RBP features exist. Skipping.")
            return

        logger.info("Computing RBP Pipeline...")

        # --- THE FIX: Ensure context exists instead of just checking for it ---
        # We pass the standard cds_windows so it matches what CUB needs later/earlier
        self._ensure_genomic_context(cds_windows=[20, 30, 40, 50, 60, 70])

        # 2. Lazy load the shared heavy assets
        rbp_map, pwm_db = self.cache.get_rbp_assets()
        unique_cell_lines = self.data[CELL_LINE_DEPMAP].dropna().unique().tolist()
        expression_matrix = self.cache.get_rbp_expression_matrix(unique_cell_lines)
        gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())

        from tauso.features.rbp.rbp_annotations import rbp_role_map_strict
        from tauso.populate.populate_rbp import (
            populate_complexity_features,
            populate_functional_features,
            populate_rbp_affinity_features,
            populate_rbp_interaction_features,
            populate_rbp_region,
        )

        # ==========================================
        # BLOCK A: Interaction Features
        # ==========================================
        if missing_int:
            logger.info("Processing Interaction Features...")
            self.data, ind_feats = populate_rbp_interaction_features(
                self.data, rbp_map, pwm_db, expression_matrix, gene_to_data, window_col, n_jobs=self.cpus
            )
            self.data, glob_feats = populate_complexity_features(
                self.data, ind_feats, suffix=str(flank_size), type="expression"
            )
            for feature in ind_feats + glob_feats:
                self._save_calculated_feature(feature_name=feature)

        # ==========================================
        # BLOCK B: Affinity Features
        # ==========================================
        if missing_aff:
            logger.info("Processing Affinity Features...")
            self.data, ind_feats = populate_rbp_affinity_features(
                self.data, rbp_map, pwm_db, gene_to_data, window_col, n_jobs=self.cpus
            )
            self.data, glob_feats = populate_complexity_features(
                self.data, ind_feats, suffix=str(flank_size), type="generic"
            )
            for feature in ind_feats + glob_feats:
                self._save_calculated_feature(feature_name=feature)

        # ==========================================
        # BLOCK C: Functional Features
        # ==========================================
        if missing_func:
            logger.info("Processing Functional Features...")
            self.data, functional_features = populate_functional_features(
                self.data, rbp_role_map_strict, flank_size=flank_size
            )
            for feature in functional_features:
                self._save_calculated_feature(feature_name=feature)

        # ==========================================
        # BLOCK D: Regional Features
        # ==========================================
        if missing_reg:
            logger.info("Processing Regional Features...")
            import warnings

            import pandas as pd
            from pandas.errors import PerformanceWarning

            from tauso.features.rbp.rbp_features import create_positional_sequence_columns

            warnings.simplefilter(action="ignore", category=PerformanceWarning)
            warnings.simplefilter(action="ignore", category=pd.errors.SettingWithCopyWarning)

            self.data = create_positional_sequence_columns(self.data, window_col, flank_size=flank_size)

            # Determine which specific regions are actually missing based on sentinels
            regions_to_run = [feat.split("_")[-1] for feat in missing_reg]

            for region in regions_to_run:
                # Fix Memory Fragmentation exactly as you specified
                self.data = self.data.copy()

                self.data, new_features = populate_rbp_region(
                    df=self.data,
                    region_name=region,
                    rbp_map=rbp_map,
                    pwm_db=pwm_db,
                    expr_matrix=expression_matrix,
                    role_map=rbp_role_map_strict,
                    gene_to_data=gene_to_data,
                    n_jobs=self.cpus,
                )

                for feature in new_features:
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
            self.calculate_modification,
            self.calculate_hybridization,
            self.calculate_backbone_features,
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
