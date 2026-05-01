import os

import pandas as pd

from tauso.data.consts import CANONICAL_GENE, CELL_LINE_DEPMAP
from tauso.features.feature_extraction import save_feature_internal
from tauso.populate.calculators.cache import AssetCache
from tauso.populate.populate_context import populate_transfection
from tauso.populate.populate_sequence import FEATURE_SPECS, populate_sequence_features
from tauso.populate.populate_structure import get_populated_df_with_structure_features
from tauso.timer import Timer


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

        if self.get_feature_dir_func is not None:
            print(f"WARNING get_feature_dir_func is None. To save features, please pass a function to the calculator.")

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

    def _get_missing_features(self, expected_features: list) -> list:
        if self.overwrite:
            return expected_features

        feature_dir = self.get_feature_dir_func(self.data_version)
        missing = []

        for feature in expected_features:
            file_path = os.path.join(feature_dir, f"{feature}.csv")
            if os.path.exists(file_path):
                print(f"Skipping '{feature}': CSV already exists.")
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
            print("Adding external mRNA and genomic context columns...")
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
            print(f"Computing {len(missing)} sequence features...")
            self.data, feature_names = populate_sequence_features(self.data, features=missing, cpus=self.cpus)

            for feature in feature_names:
                self._save_calculated_feature(feature_name=feature)
        else:
            print("All sequence features exist. Skipping.")

    def calculate_basic(self):
        """Calculates fast boolean/categorical and transfection features with dependency checking."""

        # ==========================================
        # 1. Basic Chemistry & Cell Line Features
        # ==========================================
        expected_basic = ["is_hepa", "chem_2nd_gen", "chem_3rd_gen"]
        missing_basic = self._get_missing_features(expected_basic)

        if missing_basic:
            print(f"Computing {len(missing_basic)} basic features...")

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
            print("All basic chemistry/hepa features exist. Skipping.")

        # ==========================================
        # 2. Transfection Features
        # ==========================================
        expected_transfection = ["Electroporation", "Gymnosis", "Lipofection", "Other"]
        missing_transfection = self._get_missing_features(expected_transfection)

        if missing_transfection:
            print(f"Computing {len(missing_transfection)} transfection features...")

            # 2. Check Dependency
            self._check_dependencies(["transfection_method"])

            self.data, _ = populate_transfection(self.data)

            # 4. Save
            for feature in missing_transfection:
                self._save_calculated_feature(feature_name=feature)
        else:
            print("All transfection features exist. Skipping.")

    def calculate_structure(self):
        expected_features = [
            "sense_start",
            "sense_start_from_end",
            "sense_length",
            "sense_exon",
            "sense_intron",
            "sense_utr",
            "sense_type",
        ]

        missing = self._get_missing_features(expected_features)

        if missing:
            print("Computing structure features...")

            genes_u = self._get_unique_genes()
            gene_to_data = self.cache.get_lean_gene(genes_u=genes_u)

            self.data = get_populated_df_with_structure_features(self.data, genes_u, gene_to_data)

            for feature in expected_features:
                self._save_calculated_feature(feature_name=feature)
        else:
            print("All structure features exist. Skipping.")

    def calculate_expression(self):
        """Calculates mRNA expression features."""
        # TODO: Replace with the actual list of feature names generated by populate_mrna_expression
        expected_features = ["target_expression", "rnase_expression"]
        missing = self._get_missing_features(expected_features)

        if missing:
            print(f"Computing {len(missing)} expression features...")

            self._check_dependencies([CELL_LINE_DEPMAP])
            cell_lines_depmap = self.data[CELL_LINE_DEPMAP].dropna().unique().tolist()

            transcriptomes = self.cache.get_transcriptomes(cell_lines_depmap=cell_lines_depmap)

            from tauso.populate.populate_context import populate_mrna_expression

            self.data, generated_features = populate_mrna_expression(self.data, transcriptomes)

            # Use the dynamically generated features to save,
            # ensuring we catch exactly what the function produced
            for feature in generated_features:
                self._save_calculated_feature(feature_name=feature)
        else:
            print("All expression features exist. Skipping.")

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
            print(f"Computing {len(missing)} RNase H features...")

            # Lazy import to keep initialization fast
            from tauso.populate.populate_rnase import populate_rnase_features

            # Compute the features
            self.data, generated_features = populate_rnase_features(self.data)

            # Save them securely
            for feature in generated_features:
                self._save_calculated_feature(feature_name=feature)
        else:
            print("All RNase H features exist. Skipping.")

    def calculate_off_target_single(self):
        """Calculates specific gene off-target features using the full gene dictionary."""
        # TODO: Replace with the exact string names generated by off_target_specific_seq_pandarallel
        expected_features = [
            f"off_target_single_RNASEH1_c0",
            f"off_target_single_RNASEH1_c1200",
            f"off_target_single_ACTB_c1200",
            f"off_target_single_ACTB_c1200",
        ]

        missing = self._get_missing_features(expected_features)

        if missing:
            print("Computing specific off-target features...")

            from tauso.features.hybridization_off_target.off_target_specific_gene import (
                off_target_specific_seq_pandarallel,
            )

            # Lazy load the BIG dictionary
            gene_to_data_full = self.cache.get_full_gene_data()
            cutoffs = [1200, 0]
            targets = ["RNASEH1", "ACTB"]

            for target_gene in targets:
                for cutoff in cutoffs:
                    self.data, feature_name = off_target_specific_seq_pandarallel(
                        self.data, target_gene, gene_to_data_full, cutoff=cutoff, n_jobs=self.cpus, verbose=True
                    )
                    self._save_calculated_feature(feature_name=feature_name)
        else:
            print("All specific off-target features exist. Skipping.")

    def calculate_on_target_hybridization(self):
        """Calculates on-target total hybridization features."""
        cutoffs = [1200, 0]
        expected_features = [f"on_target_total_hybridization_{c}" for c in cutoffs]

        missing = self._get_missing_features(expected_features)

        if missing:
            print(f"Computing {len(missing)} on-target hybridization features...")

            from tauso.features.hybridization_off_target.off_target_specific_gene import on_target_total_hybridization

            # Optimization: We can reuse the lean dictionary because on-target
            # only evaluates against the canonical gene of each row.
            gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())

            for cutoff in cutoffs:
                feature_name = f"on_target_total_hybridization_{cutoff}"
                if feature_name in missing:
                    self.data, generated_name = on_target_total_hybridization(
                        self.data, gene_to_data, cutoff=cutoff, n_jobs=self.cpus, verbose=True
                    )
                    print(f"Generated name: {generated_name}")
                    print(f"Feature name: {feature_name}")
                    self._save_calculated_feature(feature_name=feature_name)
        else:
            print("All on-target hybridization features exist. Skipping.")

    def calculate_mfe(self):
        """Calculates Minimum Free Energy (MFE) fold features."""
        # Lazy import to keep initialization fast
        from tauso.populate.populate_fold import DEFAULT_SETTINGS, populate_mfe_features

        # Dynamically generate the expected feature names from the settings tuples (flank, window, step)
        expected_features = [f"mfe_win{w}_flank{f}_step{s}" for f, w, s in DEFAULT_SETTINGS]

        missing = self._get_missing_features(expected_features)

        from tauso.data.consts import SENSE_LENGTH, SENSE_START

        if missing:
            print(f"Computing {len(missing)} MFE features...")

            self._check_dependencies([SENSE_START, SENSE_LENGTH])

            # Reuse the lean dictionary
            gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())

            for setting in DEFAULT_SETTINGS:
                f, w, s = setting
                feature_name = f"mfe_win{w}_flank{f}_step{s}"

                # Only run the heavy computation if this specific setting is missing
                if feature_name in missing:
                    self.data, generated_features = populate_mfe_features(
                        self.data, gene_to_data, n_jobs=self.cpus, verbose=False, settings=[setting]
                    )

                    for feature in generated_features:
                        self._save_calculated_feature(feature_name=feature)
        else:
            print("All MFE features exist. Skipping.")

    def calculate_sense_accessibility(self):
        """Calculates sense accessibility features."""
        from tauso.populate.populate_fold import DEFAULT_SENSE_CONFIGURATION, populate_sense_accessibility_batch

        # Dynamically generate expected feature names based on the configs
        expected_features = []
        for config in DEFAULT_SENSE_CONFIGURATION:
            seeds_str = "_".join(map(str, config["seeds"]))
            expected_features.append(f"access_{config['flank']}flank_{config['access']}access_{seeds_str}seed_sizes")

        missing = self._get_missing_features(expected_features)

        if missing:
            print("Computing sense accessibility features...")

            # Reuse the lean dictionary securely
            gene_to_data = self.cache.get_lean_gene(self._get_unique_genes())

            for config, expected_name in zip(DEFAULT_SENSE_CONFIGURATION, expected_features):
                if expected_name in missing:
                    c_flank = config["flank"]
                    c_access = config["access"]
                    c_seeds = config["seeds"]

                    print(f"Running: Flank={c_flank}, Access={c_access}, Seeds={c_seeds}...")

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
                    print(f"Saved: {generated_name}")
        else:
            print("All sense accessibility features exist. Skipping.")

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
            print("Computing sequence one-hot features...")
            from tauso.populate.populate_sequence import populate_sequence_one_hot_encoded

            # Pass max_len so it strictly aligns with our expected features
            self.data, generated_features = populate_sequence_one_hot_encoded(
                self.data, max_len=max_len, cpus=self.cpus
            )
            print("After one hot populate ")
            for feature in generated_features:
                if feature in missing:
                    self._save_calculated_feature(feature_name=feature)
        else:
            print("All sequence one-hot features exist. Skipping.")

    def calculate_sequence_chemistry(self):
        """Calculates sequence chemistry features."""
        from tauso.populate.populate_sequence_chemistry import FEATURE_SPECS

        expected_features = [name for name, _ in FEATURE_SPECS]
        missing = self._get_missing_features(expected_features)

        if missing:
            print(f"Computing {len(missing)} sequence chemistry features...")
            from tauso.populate.populate_sequence_chemistry import populate_sequence_chemistry_features

            # Pass only the missing features to avoid redundant calculations
            self.data, generated_features = populate_sequence_chemistry_features(
                self.data, features=missing, cpus=self.cpus
            )

            for feature in generated_features:
                self._save_calculated_feature(feature_name=feature)
        else:
            print("All sequence chemistry features exist. Skipping.")

    def calculate_hybridization(self):
        """Calculates hybridization features."""
        from tauso.populate.populate_hybridization import HYBR_FEATURE_TO_CALCULATION

        # Dynamically grab expected features
        expected_features = list(HYBR_FEATURE_TO_CALCULATION.keys())
        if "DNA_HYBR_DIFF" not in expected_features:
            expected_features.append("DNA_HYBR_DIFF")

        missing = self._get_missing_features(expected_features)

        if missing:
            print(f"Computing {len(missing)} hybridization features...")
            from tauso.populate.populate_hybridization import populate_hybridization

            # Pass `missing` so it ONLY computes the deltas
            self.data, generated_features = populate_hybridization(
                self.data, n_cores=self.cpus, features_to_run=missing
            )

            for feature in generated_features:
                self._save_calculated_feature(feature_name=feature)
        else:
            print("All hybridization features exist. Skipping.")

    def calculate_modification(self):
        """Calculates modification features."""
        from tauso.populate.populate_modification import MODIFICATION_FEATURE_TO_CALCULATION

        # Dynamically grab expected features
        expected_features = list(MODIFICATION_FEATURE_TO_CALCULATION.keys())
        missing = self._get_missing_features(expected_features)

        if missing:
            print(f"Computing {len(missing)} modification features...")
            from tauso.populate.populate_modification import populate_modifications

            # Pass `missing` to `features_to_run` so it ONLY computes the deltas
            self.data, generated_features = populate_modifications(
                self.data, n_cores=self.cpus, features_to_run=missing
            )

            for feature in generated_features:
                self._save_calculated_feature(feature_name=feature)
        else:
            print("All modification features exist. Skipping.")

    def calculate_backbone_features(self):
        """Calculates features derived from the PS_PATTERN backbone."""
        expected_features = ["is_all_ps", "po_percentage", "ps_end_score", "max_consecutive_PO"]
        missing = self._get_missing_features(expected_features)

        if missing:
            print(f"Computing {len(missing)} backbone features...")
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
            print("All backbone features exist. Skipping.")

    def calculate_ribo_seq(self):
        """Calculates ribosome profiling (Ribo-seq) features."""
        # Dynamically predict expected feature names
        flanks = (0, 10, 20, 50, 100, 125, 150)
        how = "mean"

        expected_features = [f"ribo_gene_{how}"]
        for f in flanks:
            expected_features.append(f"ribo_f{f}_{how}")
            if f > 0:
                expected_features.append(f"ribo_upstream_{f}_{how}")
                expected_features.append(f"ribo_downstream_{f}_{how}")

        missing = self._get_missing_features(expected_features)

        if missing:
            print(f"Computing {len(missing)} Ribo-seq features...")
            from tauso.features.context.ribo_seq import add_genomic_coordinates
            from tauso.populate.populate_context import populate_ribo_seq

            mapper = self.cache.get_gene_mapper()
            self.data = add_genomic_coordinates(self.data, mapper)

            # Pass the arguments to match our expected names
            self.data, generated_features = populate_ribo_seq("human", self.data, flanks=flanks, how=how)

            for feature in generated_features:
                if feature in missing:
                    self._save_calculated_feature(feature_name=feature)
        else:
            print("All Ribo-seq features exist. Skipping.")

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
            print("All Codon Usage Bias (CUB) features exist. Skipping.")
            return

        print("Computing Codon Usage Bias (CUB) features...")

        from tauso.data.consts import SENSE_START

        self._check_dependencies([SENSE_START])

        # 4. Load the shared heavy dependencies EXACTLY ONCE
        registry = self.cache.get_gene_registry(self._get_unique_genes())
        self._ensure_genomic_context(cds_windows)

        # 5. Execute conditionally based on what is missing
        if missing_tai:
            print(f"  -> Computing {len(missing_tai)} tAI features...")
            from tauso.populate.populate_codon_usage import populate_tai

            self.data, generated_features = populate_tai(self.data, cds_windows, registry)
            for feature in generated_features:
                if feature in missing_tai:
                    self._save_calculated_feature(feature_name=feature)

        if missing_cai:
            print(f"  -> Computing {len(missing_cai)} CAI features...")
            from tauso.populate.populate_codon_usage import populate_cai

            self.data, generated_features = populate_cai(self.data, cds_windows, registry, n_jobs=self.cpus)
            for feature in generated_features:
                if feature in missing_cai:
                    self._save_calculated_feature(feature_name=feature)

        if missing_enc:
            print(f"  -> Computing {len(missing_enc)} ENC features...")
            from tauso.populate.populate_codon_usage import populate_enc

            self.data, generated_features = populate_enc(self.data, cds_windows, registry, n_jobs=self.cpus)
            for feature in generated_features:
                if feature in missing_enc:
                    self._save_calculated_feature(feature_name=feature)

    def calculate_off_target_general(self):
        """Calculates general off-target hybridization scores."""
        from tauso.features.hybridization_off_target.add_off_target_feat import AggregationMethod
        from tauso.features.hybridization_off_target.off_target_feature import serialize_feature_name

        # Define the parameter spaces
        methods = [AggregationMethod.ARTM, AggregationMethod.MECH]
        top_ns = [25, 50, 100, 200]
        cutoffs = [1200, 1000, 800]

        # Generate all combinations dynamically
        configs = [(m, n, c) for m in methods for n in top_ns for c in cutoffs]

        expected_features = [serialize_feature_name(m, n, c, is_specific=False) for m, n, c in configs]
        missing = self._get_missing_features(expected_features)

        if missing:
            print(f"Computing {len(missing)} general off-target features...")
            from tauso.features.hybridization_off_target.off_target_feature import populate_off_target_general

            # Load the heavy dictionaries (happens instantly if already in memory)
            gene_to_data = self.cache.get_full_gene_data()

            self._check_dependencies([CELL_LINE_DEPMAP])
            cell_lines_depmap = self.data[CELL_LINE_DEPMAP].dropna().unique().tolist()

            transcriptomes = self.cache.get_transcriptomes(cell_lines_depmap=cell_lines_depmap)

            for method, top_n, cutoff in configs:
                feat_name = serialize_feature_name(method, top_n, cutoff, is_specific=False)

                # Only execute if this specific configuration is missing
                if feat_name in missing:
                    self.data, generated_features = populate_off_target_general(
                        ASO_df=self.data,
                        gene_to_data=gene_to_data,
                        cell_line2data=transcriptomes,
                        top_n_list=[top_n],
                        cutoff_list=[cutoff],
                        method=method,
                        n_cores=self.cpus,
                    )

                    for feature in generated_features:
                        self._save_calculated_feature(feature_name=feature)
        else:
            print("All general off-target features exist. Skipping.")

    def calculate_off_target_specific(self):
        """Calculates cell-line specific off-target hybridization scores."""
        from tauso.features.hybridization_off_target.add_off_target_feat import AggregationMethod
        from tauso.features.hybridization_off_target.off_target_feature import serialize_feature_name

        method = AggregationMethod.ARTM
        top_n_list = [50, 100, 200]
        cutoff_list = [800, 1000, 1200]

        # Generate expected combinations
        expected_features = [
            serialize_feature_name(method, n, c, is_specific=True) for n in top_n_list for c in cutoff_list
        ]

        missing = self._get_missing_features(expected_features)

        if missing:
            print(f"Computing {len(missing)} specific off-target features...")
            from tauso.features.hybridization_off_target.off_target_feature import populate_off_target_specific

            gene_to_data = self.cache.get_full_gene_data()
            self._check_dependencies([CELL_LINE_DEPMAP])
            cell_lines_depmap = self.data[CELL_LINE_DEPMAP].dropna().unique().tolist()

            transcriptomes = self.cache.get_transcriptomes(cell_lines_depmap=cell_lines_depmap)

            for top_n in top_n_list:
                for cutoff in cutoff_list:
                    feat_name = serialize_feature_name(method, top_n, cutoff, is_specific=True)

                    if feat_name in missing:
                        self.data, generated_features = populate_off_target_specific(
                            ASO_df=self.data,
                            gene_to_data=gene_to_data,
                            cell_line2data=transcriptomes,
                            top_n_list=[top_n],
                            cutoff_list=[cutoff],
                            method=method,
                            n_cores=self.cpus,
                        )

                        for feature in generated_features:
                            self._save_calculated_feature(feature_name=feature)
        else:
            print("All specific off-target features exist. Skipping.")

    def calculate_mrna_halflife(self):
        """Calculates mRNA stability and half-life features."""
        expected_features = ["mRNA_HalfLife", "HalfLife_Source", "Mapped_Cell_Proxy"]
        missing = self._get_missing_features(expected_features)

        if missing:
            print(f"Computing {len(missing)} mRNA half-life features...")

            # Check dependencies BEFORE loading heavy objects
            from tauso.data.consts import CANONICAL_GENE, CELL_LINE

            self._check_dependencies([CANONICAL_GENE, CELL_LINE])

            from tauso.features.context.mrna_halflife import populate_mrna_halflife_features

            # 1. Lazy load the provider
            provider = self.cache.get_halflife_provider()

            # 2. Compute
            self.data, generated_features = populate_mrna_halflife_features(self.data, provider)

            # 3. Save only what was missing
            for feature in generated_features:
                if feature in missing:
                    self._save_calculated_feature(feature_name=feature)
        else:
            print("All mRNA half-life features exist. Skipping.")

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
            print("All RBP features exist. Skipping.")
            return

        print("Computing RBP Pipeline...")

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
            print("\n--- Processing Interaction Features ---")
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
            print("\n--- Processing Affinity Features ---")
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
            print("\n--- Processing Functional Features ---")
            self.data, functional_features = populate_functional_features(
                self.data, rbp_role_map_strict, flank_size=flank_size
            )
            for feature in functional_features:
                self._save_calculated_feature(feature_name=feature)

        # ==========================================
        # BLOCK D: Regional Features
        # ==========================================
        if missing_reg:
            print("\n--- Processing Regional Features ---")
            import warnings

            import pandas as pd
            from pandas.errors import PerformanceWarning

            from tauso.features.rbp.RBP_features import create_positional_sequence_columns

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
            msg = f"Found {len(unique_errors)} unique error types:"
            print(msg)
            err_str = msg
            for err in unique_errors:
                msg = f" - {err}"
                print(msg)
                err_str += "\n" + msg
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

        print(f"=== Starting pipeline with {len(pipeline_steps)} steps ===")

        for step in pipeline_steps:
            # step.__name__ dynamically grabs the name of the function (e.g., 'calculate_cub')
            print("Starting step: ", step.__name__)
            try:
                with Timer(name=step.__name__):
                    step()
            except Exception as _:
                print(f"\nCalculator crashed during step: {step.__name__}\n")
                raise

        print("=== Pipeline Complete ===")
