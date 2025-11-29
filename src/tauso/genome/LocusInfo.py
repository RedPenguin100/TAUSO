class LocusInfo:
    def __init__(self):
        self.exons = []
        self.introns = []
        self.exon_indices = []
        self.intron_indices = []
        self.stop_codons = []
        self.five_prime_utr = ""
        self.three_prime_utr = ""
        self.exon_concat = None
        self.full_mrna = None
        self.cds_start = None
        self.cds_end = None
        self.strand = None
        self.gene_type = None
        self.utr_indices = []
