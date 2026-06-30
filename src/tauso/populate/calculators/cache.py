import logging

from ...genome.LocusInfo import LocusInfo
from ...genome.read_human_genome import get_locus_to_data_dict

logger = logging.getLogger(__name__)


class AssetCache:
    def __init__(self, genome="GRCh38"):
        self.genome = genome

        self._gene_to_data_full = None
        self.custom_gene = None  #  not None if user decides to pass a FASTA file

        # RBP
        self._rbp_map = None
        self._pwm_db = None

        self._halflife_provider = None
        self._transcriptomes = None

        self._gene_to_data_lean = None
        self._genes_u = None

        self._gene_registry = None
        self._context_added = False

    def preload(self):
        self.get_full_gene_data()
        self.get_halflife_provider()

        self.get_rbp_assets()

        # Omitting gene_to_data_lean, get_transcriptomes, get_gene_registry
        # because they require more context

    def get_genome(self):
        return self._genome

    def get_full_gene_data(self):
        """Lazy getter for the full genome data dictionary (no subset)."""
        if self._gene_to_data_full is None:
            logger.info("Loading FULL genome data dictionary into memory (happens once)...")

            self._gene_to_data_full = get_locus_to_data_dict(include_introns=True, genome=self.genome)
            if self.custom_gene is not None:
                name, locus_info = self.custom_gene
                self._gene_to_data_full[name] = locus_info
        return self._gene_to_data_full

    def get_rbp_assets(self):
        """Lazy loader for ATTRACT RBP map and PWM database."""
        if self._rbp_map is None or self._pwm_db is None:
            logger.info("Loading ATTRACT RBP data into memory (happens once)...")
            from tauso.features.rbp.load_rbp import load_attract_data

            self._rbp_map, self._pwm_db = load_attract_data()
        return self._rbp_map, self._pwm_db

    def get_halflife_provider(self):
        """Lazy loader for the mRNA Half-Life Provider."""
        if self._halflife_provider is None:
            logger.info("Loading mRNA Half-Life Oracle into memory (happens once)...")
            from tauso.features.context.mrna_halflife import HalfLifeProvider, load_halflife_mapping

            mapping = load_halflife_mapping()
            self._halflife_provider = HalfLifeProvider(mapping)

        return self._halflife_provider

    def get_transcriptomes(self, cell_lines_depmap):
        """Lazy loader for the FULL transcriptomes (required for off-target scanning)."""
        if self._transcriptomes is None:
            logger.info("Loading FULL transcriptomes into memory (happens once)...")
            import os
            from pathlib import Path

            from tauso.common.gtf import filter_gtf_genes
            from tauso.data.data import get_data_dir, load_gtf_db
            from tauso.features.codon_usage.find_cai_reference import load_cell_line_gene_expression
            from tauso.features.expression.general_expression import get_general_expression_of_genes

            db = load_gtf_db()
            valid_genes = filter_gtf_genes(db, filter_mode="non_mt")

            data_dir = get_data_dir()
            expression_dir = os.path.join(data_dir, "processed_expression")

            self._transcriptomes = load_cell_line_gene_expression(
                cell_lines_depmap,
                valid_genes,
                expression_dir=expression_dir,
            )

            mean_exp_data = get_general_expression_of_genes(
                Path(data_dir) / "OmicsExpressionTPMLogp1HumanAllGenesStranded.parquet", valid_genes
            )
            self._transcriptomes["general"] = mean_exp_data

        return self._transcriptomes

    def get_lean_gene(self, genes_u):
        if self._gene_to_data_lean is not None and genes_u == self._genes_u:
            # if self._genes_u is None, then the data object is None as well and we skip the if statement
            return self._gene_to_data_lean

        if genes_u != self._genes_u:
            logger.warning("Genes %s differ from cached genes %s, regenerating cache.", genes_u, self._genes_u)

        logger.info("Loading genome data dictionary into memory (happens once)...")
        self._gene_to_data_lean = get_locus_to_data_dict(
            include_introns=True, gene_subset=genes_u, genome=self.genome, canonical_only=True
        )
        if self.custom_gene is not None:
            name, locus_info = self.custom_gene
            self._gene_to_data_lean[name] = locus_info
        self._genes_u = genes_u

        return self._gene_to_data_lean

    def get_gene_registry(self, genes_u):
        """Lazy loader for the built gene sequence registry."""
        if self._gene_registry is not None and genes_u == self._genes_u:
            return self._gene_registry

        if genes_u != self._genes_u:
            logger.warning("Genes %s differ from cached genes %s, regenerating cache.", genes_u, self._genes_u)

        logger.info("Building gene sequence registry (happens once)...")
        from tauso.genome.TranscriptMapper import build_gene_sequence_registry

        gene_to_data = self.get_lean_gene(genes_u)

        self._gene_registry = build_gene_sequence_registry(genes=genes_u, gene_to_data=gene_to_data)
        return self._gene_registry

    def set_custom_gene(self, name, sequence):
        """Register a user-provided gene sequence, replacing any previous custom gene and
        invalidating the cached gene data so the new sequence is used on the next access
        (otherwise a reused cache would silently keep scoring the first gene)."""
        self.custom_gene = (name, LocusInfo(seq=sequence))
        self._gene_to_data_full = None
        self._gene_registry = None
