class LocusInfo:
    def __init__(self, seq=None):
        # Memory Optimization: We no longer store lists of strings for exons/introns.
        # These are computed on the fly via @property decorators.
        self.exon_indices = []
        self.intron_indices = []
        self.stop_codons = []
        self.five_prime_utr = ""
        self.three_prime_utr = ""
        self.full_mrna = None
        self.gene_start = None
        self.gene_end = None
        self.strand = None
        self.gene_type = None
        self.utr_indices = []

        if seq is not None:
            self.exon_indices = [(0, len(seq))]
            self.gene_start = 0
            self.gene_end = len(seq)
            self.strand = "+"
            self.gene_type = "unknown"
            self.full_mrna = seq

    def _get_sequence_slice(self, start, end):
        """Helper to extract a sequence from the full_mrna based on genomic coordinates."""
        if not self.full_mrna or self.gene_start is None or self.gene_end is None:
            return ""

        if self.strand == "+":
            # Relative to the 5' end of the + strand (gene_start)
            rel_start = start - self.gene_start
            rel_end = end - self.gene_start
        else:
            # Relative to the 5' end of the - strand (gene_end)
            # Because full_mrna is reverse complemented, we map from the right side.
            rel_start = self.gene_end - end
            rel_end = self.gene_end - start

        return self.full_mrna[rel_start:rel_end]

    @property
    def exons(self):
        return [self._get_sequence_slice(start, end) for start, end in self.exon_indices]

    @property
    def introns(self):
        return [self._get_sequence_slice(start, end) for start, end in self.intron_indices]

    @property
    def exon_concat(self):
        return "".join(self.exons)

    def __repr__(self):
        print("LocusInfo:")
        # We use a custom print to avoid aggressively computing properties just for a representation
        for field, value in self.__dict__.items():
            print(f"  {field}: {value}")
        return ""
