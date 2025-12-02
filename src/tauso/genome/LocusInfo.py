class LocusInfo:
    def __init__(self, seq=None):
        self.exons = []
        self.introns = []
        self.exon_indices = []
        self.intron_indices = []
        self.stop_codons = []
        self.five_prime_utr = ""
        self.three_prime_utr = ""
        self.exon_concat = None
        self.full_mrna = None
        self.gene_start = None
        self.gene_end = None
        self.strand = None
        self.gene_type = None
        self.utr_indices = []

        if seq is not None:
            self.exons = [seq]
            self.exon_indices = [(0, len(seq) - 1)]
            self.gene_start = 0
            self.gene_end = len(seq) - 1
            self.strand = "+"
            self.gene_type = "unknown"
            self.exon_concat = seq
            self.full_mrna = seq

    def __repr__(self):
        print("LocusInfo:")
        for field, value in self.__dict__.items():
            print(f"  {field}: {value}")