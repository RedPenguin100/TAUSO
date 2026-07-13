"""Curated corrections to the patent-provided target gene, applied in 2_assign_canonical_gene.

Each ASO's gene comes from its patent label (target_gene); these entries fix known label errors
found by curate_gene_labels.py. Genes mapped to NA are unrecoverable and become NaN (dropped in
3_process_data). Genes mapped to themselves keep the patent's original over a conflicting alignment.
"""
import pandas as pd

MANUAL_CANONICAL_MAPPING = {
    # Original is wrong (compilation error)
    "SFRS10": "TRA2B",      # outdated alias
    "SENP2": "SUMO2",       # patent targeted SUMO2, confused with the enzyme SENP2 that cuts it
    "SRC": "CSK",           # patent targeted C-SRC kinase (CSK); compiler stopped at SRC
    "NRF1": "NKRF",         # patent targeted NF-kappa-B-repressing factor (NKRF), not NRF1
    "BRCA2": "N4BP2L2",     # BRCA2-region transcription unit CG005, later renamed N4BP2L2
    "MAPK6": "MAPK12",      # patent targeted ERK6 (MAPK12); confused with ERK3
    "NFKBIL1": "TONSL",     # IKB-R, old name for NFKBIL2, officially TONSL
    "PTPRJ": "PTPRF",       # patent target LAR (alias of PTPRF), mixed up with PTPRJ
    "THRAP3": "ZNHIT3",     # patent targeted TRIP3 (alias of ZNHIT3)
    # Original is wrong (patent / mislabel)
    "UBE3A": "SNHG14",      # patent discusses UBE3A regulation but targets SNHG14
    "RASD1": "RASD2",       # patent annotated incorrectly
    # Target named only in target_mrna (target_gene left blank); recovered via the fallback in
    # 2_assign_canonical_gene. Each verified by aligning the ASOs to the named gene body.
    "PKK": "KLKB1",                                 # plasma prekallikrein = KLKB1 (444/467 align)
    "GCCR": "NR3C1",                                # glucocorticoid receptor = NR3C1 (148/148)
    "estrogen-responsive finger protein": "TRIM25",  # EFP = TRIM25 (73/78)
    "KOX1": "ZNF10",                                # KOX1 = ZNF10 (72/78)
    "PPAR binding protein": "MED1",                 # PBP/TRAP220 = MED1 (72/77)
    "fetoprotein transcription factor": "NR5A2",    # FTF/LRH-1 = NR5A2 (66/75)
    "ABC transporter MHC 1": "TAP1",               # TAP1 (73/78)
    "endothelial differentiation gene 2": "LPAR1",  # EDG2 = LPAR1 (67/72)
    "jumonji": "JARID2",                            # jumonji = JARID2 (33/37)
    "E2-EPF": "UBE2S",                             # E2-EPF ubiquitin carrier = UBE2S (32/37)
    "FBP-interacting repressor": "PUF60",           # FIR = PUF60 (32/37)
    "jerky-like 1": "JRKL",                        # JRKL (37/37)
    "HIP-1 protein interactor": "IFT57",            # HIPPI = IFT57 (75/78)
    "zinedin": "STRN4",                            # zinedin = STRN4 (68/78)
    "cytokine-inducible kinase": "PLK3",            # FNK/PRK = PLK3 (61/78; 0/78 to PIM1)
    # The most popular alignment disagrees with the original; we keep the original
    "ANGPT2": "ANGPT2",
    "ANGPTL3": "ANGPTL3",
    "BAG5": "BAG5",
    "FNTB": "FNTB",
    "IGF2": "IGF2",
    "IL18": "IL18",
    "MMP1": "MMP1",
    "MRPS16": "MRPS16",
    "PLN": "PLN",
    "S1PR2": "S1PR2",
    "STAT5B": "STAT5B",
    # Patent errors, unrecoverable -> dropped
    "RHO": pd.NA,           # patent targets a mutation, not the reference sequence
    "NAIP": pd.NA,          # old, mis-annotated genome
}
