import re

SEQUENCE = "Sequence"
INHIBITION = "Inhibition(%)"
CANONICAL_GENE = "Canonical Gene Name"
CELL_LINE_ORGANISM = "Cell line organism"
VOLUME = "ASO_volume(nM)"
CHEMICAL_PATTERN = "Chemical_Pattern"
TREATMENT_PERIOD = "Treatment_Period(hours)"
CELL_LINE = "Cell_line"
TRANSFECTION = "Transfection"
DENSITY = "Density(cells/well)"
DENSITY_UPDATED = "Density(cells_per_well)"  # Avoiding /
PRIMER_PROBE = "Primer_probe_set"
MODIFICATION = "Modification"
SENSE_START = "sense_start"
SENSE_LENGTH = "sense_length"
SENSE_TYPE = "sense_type"
SENSE_SEQUENCE = "sense_sequence"
PRE_MRNA_SEQUENCE = "pre_mrna_sequence"
SENSE_AVG_ACCESSIBILITY = "sense_avg_accessibility"
CELL_LINE_DEPMAP_PROXY = "Cell_Line_Depmap_Proxy"
CELL_LINE_DEPMAP = "Cell_Line_Depmap"

# OligoAI syntax
SUGAR_MODS = 'sugar_mods'
BACKBONE_MODS = 'backbone_mods'

ISIS = "ISIS"
LINKAGE = "Linkage"
LINKAGE_LOCATION = "Linkage_Location"
TARGET_GENE = "Target_gene"
SMILES = "Smiles"

PS_PATTERN = "ps_pattern"


def standardize_cell_line_name(name: str) -> str:
    """
    Standardizes cell line names to alphanumeric uppercase.
    Example: 'A-431' -> 'A431'
    """
    if not name or not isinstance(name, str):
        return str(name)
    if name == "U-251" or name == "U251":
        return "U251MG"
    return re.sub(r"[^a-zA-Z0-9]", "", name).upper()


HEPG2 = "HepG2"  # Note: DepMap usually lists this as HepG2, but keeping your key 'HepG'
SNU_449 = "SNU-449"
HELA = "HeLa"
A431 = "A431"
SK_MEL_28 = "SK-MEL-28"
SH_SY5Y = "SH-SY5Y"
U251 = "U251"
U_251_MG = "U-251 MG"
H929 = "H929"
KMS11 = "KMS11"
NCI_H460 = "NCI-H460"
SK_N_AS = "SK-N-AS"
SK_N_SH = "SK-N-SH"
KARPAS_229 = "KARPAS-229"

# 1. Map Input Names -> Canonical "Reasonable" Proxy
CELL_LINE_TO_DEPMAP_PROXY_DICT = {
    # --- Direct Matches & Normalization ---
    "A-431": "A-431",
    "A431": "A-431",
    "A-549": "A549",
    "A549": "A549",
    "A459": "A549",  # Typo Fixed
    "A172": "A-172",
    "G-361": "G-361",
    "H929": "NCI-H929",
    "HeLa": "HeLa",
    "Hela": "HeLa",
    "Hep3B": "Hep 3B2.1-7",
    "HepB3": "Hep 3B2.1-7",  # Typo Fixed
    "HepG2": "Hep G2",
    'HepG2/Hep3B': 'Hep G2',                # Split: Mapping to dominant line
    "Huh7": "HuH-7",
    "HK-2": "HK-2",
    "Jurkat": "JURKAT",
    "K-562": "K-562",
    "KARPAS-229": "Karpas-299",  # Typo Fixed
    "KMS11": "KMS-11",
    "LNCaP": "LNCaP clone FGC",
    "MCF7": "MCF7",
    "MM.1R": "MM.1S",  # Proxy: Parental line
    "MM.1S": "MM.1S",
    "NCI-H460": "NCI-H460",
    "PC3": "PC-3",
    "SH-SY5Y": "SH-SY5Y",
    "SH-SY-5Y": "SH-SY5Y",  # Typo Fixed
    "SK-MEL-28": "SK-MEL-28",
    "SK-N-AS": "SK-N-AS",
    "SK-N-SH": "SK-N-SH",
    "SKOV3": "SK-OV-3",
    "SNU-449": "SNU-449",
    "SW872": "SW872",
    "T-24": "T24",
    "T24": "T24",
    "THP-1": "THP-1",
    "U251": "U-251 MG",
    "U-251 MG": "U-251 MG",
    "VCaP": "VCaP",
    # --- Experimental Contexts (User Defined) ---
    "Angptl2/Actin": "SK-N-AS",
    "SK cells asyn": "SK-N-AS",
    "PAC neurons asyn": "SH-SY5Y",
    "Human Neuronal Cell": "SH-SY5Y",
    "iCell GABANeurons": "SH-SY5Y",           # Neuronal proxy
    "ReproNeuro": "SH-SY5Y",                 # Neuronal proxy
    "ReproNeuro Neurons": "SH-SY5Y",
    "ReproNeuro neurons (ReproCELL)": "SH-SY5Y",
    "ReproNeuro Neurons (ReproCELL)" : "SH-SY5Y",
    "54-2": None,                            # Keep as None unless you confirm lineage

    # --- No Valid Proxy (Primary Cells / Distinct Lineages) ---
    "HepaRG": None,
    "HuVEC": None,
    "HUVEC": None,
    "hSKM": None,
    "hSKMc": None,
    "SCA2-04": None,
    "differentiated human adipocytes": None,
    "iCell cardiomyocytes2": None,
    "iCell cardiomyocytes (R1017)": None,
}

HEPA_PROXIES = {"Hep 3B2.1-7", "Hep G2", "HuH-7", "SNU-449", "HepaRG"}


# 2. Map Canonical Proxy -> DepMap ID
CELL_LINE_TO_DEPMAP = {
    "A-172": "ACH-000558",
    "A-431": "ACH-001328",
    "A549": "ACH-000681",
    "G-361": "ACH-000572",
    "HK-2": "ACH-001087",
    "HeLa": "ACH-001086",
    "Hep 3B2.1-7": "ACH-000625",
    "Hep G2": "ACH-000739",
    "HuH-7": "ACH-000480",
    "JURKAT": "ACH-000995",
    "K-562": "ACH-000551",
    "Karpas-299": "ACH-000053",
    "KMS-11": "ACH-000714",
    "LNCaP clone FGC": "ACH-000977",
    "MCF7": "ACH-000019",
    "MM.1S": "ACH-000763",
    "NCI-H460": "ACH-000463",
    "NCI-H929": "ACH-000050",
    "PC-3": "ACH-000090",
    "SH-SY5Y": "ACH-001188",
    "SK-MEL-28": "ACH-000615",
    "SK-N-AS": "ACH-000260",
    "SK-N-SH": "ACH-000149",
    "SK-OV-3": "ACH-000811",
    "SNU-449": "ACH-000420",
    "SW872": "ACH-002310",
    "T24": "ACH-000018",
    "THP-1": "ACH-000146",
    "U-251 MG": "ACH-000232",
    "VCaP": "ACH-000115",
}

CELL_LINE_MAPPING = {
    # Direct Matches & Synonyms
    "HepG2": HEPG2,
    "HepaRG": HEPG2,  # Liver -> Liver proxy
    "SNU-449": SNU_449,
    "Hela": HELA,
    "A431": A431,
    "A-431": A431,
    "SK-MEL-28": SK_MEL_28,
    "SH-SY5Y": SH_SY5Y,
    "U251": U251,
    "U-251 MG": U251,  # Map both to the canonical U251
    "H929": H929,
    "KMS11": KMS11,
    "NCI-H460": NCI_H460,
    "SK-N-AS": SK_N_AS,
    "SK-N-SH": SK_N_SH,
    "KARPAS-229": "KARPAS299",  # NOT A TYPO, a proxy!
    # Neuronal / Experimental Contexts
    "Angptl2/Actin": SK_N_AS,
    "SK cells asyn": SK_N_AS,  # "SK cells" usually refers to SK-N-AS in this context
    "PAC neurons asyn": SH_SY5Y,  # SH-SY5Y is the standard neuronal proxy here
    "Human Neuronal Cell": SH_SY5Y,  # Generic neuronal -> SH-SY5Y proxy
    # Blood / Immune Proxies
    "MM.1R": H929,  # Myeloma -> H929 (Myeloma) proxy
}

HUMAN = "human"
MONKEY = "monkey"
MOUSE = "mouse"
UNKNOWN = "unknown"


cell_line_to_organism = {
    # --- Human ---
    "A431": HUMAN,
    "A-431": HUMAN,
    "Angptl2/Actin": HUMAN,
    "CC-2580": HUMAN,
    "H929": HUMAN,
    "Hela": HUMAN,
    "HepG2": HUMAN,
    "HepaRG": HUMAN,
    "Human IPS": HUMAN,
    "Human Neuronal Cell": HUMAN,
    "KARPAS-229": HUMAN,
    "KMS11": HUMAN,
    "MM.1R": HUMAN,
    "NCI-H460": HUMAN,
    "PAC neurons asyn": HUMAN,
    "SH-SY5Y": HUMAN,
    "SK cells asyn": HUMAN,
    "SK-MEL-28": HUMAN,
    "SNU-449": HUMAN,
    "U251": HUMAN,
    # --- Mouse ---
    "Mouse Primary Neuronal Cell": MOUSE,
    "Oscillation": MOUSE,
    "Tubulin": MOUSE,
    # --- Monkey ---
    "LLC-MK2 monkey": MONKEY,  # Green Monkey
}
