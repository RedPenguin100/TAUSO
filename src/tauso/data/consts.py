import re

SEQUENCE = "Sequence"
INHIBITION = "Inhibition(%)"
CANONICAL_GENE = "Canonical Gene Name"
CELL_LINE_ORGANISM = "Cell line organism"
VOLUME = "ASO_volume(nM)"
CHEMICAL_PATTERN = "Chemical_Pattern"  # LNA / MOE / cEt
TREATMENT_PERIOD = "Treatment_Period(hours)"
CELL_LINE = "Cell_line"
TRANSFECTION = "Transfection"
DENSITY = "Density(cells/well)"
DENSITY_UPDATED = "Density(cells_per_well)"  # Avoiding /
PRIMER_PROBE = "Primer_probe_set"
MODIFICATION = "Modification"  # MMMdddMMM or CCCdddCCC etc
SENSE_START = "sense_start"
SENSE_START_FROM_END = "sense_start_from_end"
SENSE_LENGTH = "sense_length"
SENSE_TYPE = "sense_type"
SENSE_SEQUENCE = "sense_sequence"
SENSE_EXON = "sense_exon"
SENSE_INTRON = "sense_intron"
SENSE_UTR = "sense_utr"
SENSE_3UTR = "sense_3utr"
SENSE_5UTR = "sense_5utr"
SENSE_CDS = "sense_cds"
SENSE_START_NORM = "sense_start_norm"
SENSE_START_FROM_END_NORM = "sense_start_from_end_norm"
SENSE_DIST_TO_CANONICAL_STOP = "sense_dist_to_canonical_stop"
SENSE_DIST_TO_CLOSEST_STOP = "sense_dist_to_closest_stop"
SENSE_DIST_TO_CANONICAL_START = "sense_dist_to_canonical_start"
SENSE_DIST_TO_CLOSEST_START = "sense_dist_to_closest_start"
SENSE_MRNA_DIST_TO_CANONICAL_STOP = "sense_mrna_dist_to_canonical_stop"
SENSE_MRNA_DIST_TO_CLOSEST_STOP = "sense_mrna_dist_to_closest_stop"
PRE_MRNA_SEQUENCE = "pre_mrna_sequence"
SENSE_AVG_ACCESSIBILITY = "sense_avg_accessibility"
CELL_LINE_DEPMAP_PROXY = "Cell_Line_Depmap_Proxy"
CELL_LINE_DEPMAP = "Cell_Line_Depmap"

# OligoAI syntax
SUGAR_MODS = "sugar_mods"
BACKBONE_MODS = "backbone_mods"

ISIS = "ISIS"
LINKAGE = "Linkage"
LINKAGE_LOCATION = "Linkage_Location"
TARGET_GENE = "Target_gene"
SMILES = "Smiles"

PS_PATTERN = "ps_pattern"


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

# Map dataset cell-line spellings -> canonical DepMap proxy. Lookup is
# case- and punctuation-insensitive (see _norm_cell_line_key), so a single
# entry covers "HeLa" / "Hela" / "HELA" / " hela ", "A-431" / "A431", etc.
# A value of None means "no honest proxy" (primary cells, iPSC-derived models,
# near-isogenic mismatches) and resolves to a NaN DepMap ID downstream.
CELL_LINE_TO_DEPMAP_PROXY_DICT = {
    # --- Cancer lines with valid DepMap proxies ---
    "A-431": "A-431",
    "A-549": "A549",
    "A459": "A549",  # Typo
    "A172": "A-172",
    "G-361": "G-361",
    "H929": "NCI-H929",
    "HeLa": "HeLa",
    "Hep3B": "Hep 3B2.1-7",
    "HepB3": "Hep 3B2.1-7",  # Typo
    "HepG2": "Hep G2",
    "HepG2/Hep3B": "Hep G2",  # Combined experiment; map to dominant line
    "Huh7": "HuH-7",
    "HK-2": "HK-2",
    "Jurkat": "JURKAT",
    "K-562": "K-562",
    "KMS-11": "KMS-11",
    "LNCaP": "LNCaP clone FGC",
    "MCF7": "MCF7",
    "MM.1R": "MM.1S",  # Proxy: parental dex-sensitive line
    "MM.1S": "MM.1S",
    "NCI-H460": "NCI-H460",
    "PC3": "PC-3",
    "SH-SY5Y": "SH-SY5Y",
    "SK-MEL-28": "SK-MEL-28",
    "SK-N-AS": "SK-N-AS",
    "SK-N-SH": "SK-N-SH",
    "SKOV3": "SK-OV-3",
    "SNU-449": "SNU-449",
    "SW872": "SW872",
    "T24": "T24",
    "THP-1": "THP-1",
    "U251": "U-251 MG",
    "VCaP": "VCaP",
    # --- Experimental labels mapped by domain knowledge ---
    "Angptl2/Actin": "SK-N-AS",
    "SK cells asyn": "SK-N-AS",
    # --- Explicitly blocked: no honest cancer-line proxy ---
    # iPSC-derived and primary neuronal models, primary cells, lineages too far
    # from any DepMap line. iCell GABANeurons and ReproNeuro all target UBE3A
    # via gymnosis; assigning SH-SY5Y would import a wrong transcriptome.
    "PAC neurons asyn": None,
    "Human Neuronal Cell": None,
    "iCell GABANeurons": None,
    # ReproNeuro appears under 4 spellings for the same 78 ASOs; the "(ReproCELL)"
    # suffix variants don't collapse with the base name under punctuation-strip,
    # so each spelling needs its own blocked-list entry to surface as "Generic".
    "ReproNeuro": None,
    "ReproNeuro Neurons": None,
    "ReproNeuro Neurons (ReproCELL)": None,
    "ReproNeuro neurons (ReproCELL)": None,
    "54-2": None,  # Lineage unconfirmed
    "HepaRG": None,
    "HUVEC": None,
    "hSKM": None,
    "hSKMc": None,
    "SCA2-04": None,
    "differentiated human adipocytes": None,
    "iCell cardiomyocytes2": None,
    "iCell cardiomyocytes (R1017)": None,
    "HEK293T": None,  # Tempting to map to HEK293 but they diverge enough that the proxy lies.
    "KARPAS-229": None,  # B-cell line, not in DepMap. Karpas-299 is a different lineage (T-cell ALCL).
}


def _norm_cell_line_key(raw):
    """Lower-case + strip every non-alphanumeric so HeLa/Hela/HELA/'A-431'/A431 collide."""
    if not isinstance(raw, str):
        return None
    return re.sub(r"[^a-z0-9]", "", raw.lower())


_PROXY_LOOKUP_NORM = {_norm_cell_line_key(k): v for k, v in CELL_LINE_TO_DEPMAP_PROXY_DICT.items()}


def resolve_depmap_proxy(raw):
    """Return the canonical DepMap proxy spelling for a dataset cell-line name.

    Lookup is case- and punctuation-insensitive. Returns:
      - the canonical proxy string when one is registered,
      - None when the entry is explicitly blocked (primary cell / iPSC / etc.)
        OR the name is unknown -- treated the same so unknown spellings
        propagate as a NaN DepMap ID rather than silently using the raw input.
    """
    key = _norm_cell_line_key(raw)
    if key is None or key not in _PROXY_LOOKUP_NORM:
        return None
    return _PROXY_LOOKUP_NORM[key]


def resolve_depmap_id(raw):
    """Resolve a raw cell-line name straight to its DepMap ACH-... ID, or None."""
    proxy = resolve_depmap_proxy(raw)
    if proxy is None:
        return None
    return CELL_LINE_TO_DEPMAP.get(proxy)


def standardize_cell_line_name(raw: str) -> str:
    """Return an UPPERCASE-stripped tag used by the codon-usage scorer cache.

    Distinct from resolve_depmap_proxy: that one returns the canonical proxy
    spelling for DepMap-merge, this one returns a cache-key tag for scorer
    lookups. Both go through the same case/punctuation-tolerant proxy table.
    Explicitly-blocked entries (HepaRG, HEK293T, ...) and missing input collapse
    to "Generic" so the scorer falls back to the cross-cell-line consensus.
    """
    if not raw or not isinstance(raw, str) or raw.strip().lower() == "none":
        return "Generic"

    key = _norm_cell_line_key(raw)
    if key in _PROXY_LOOKUP_NORM:
        proxy = _PROXY_LOOKUP_NORM[key]
        if proxy is None:
            return "Generic"  # Explicitly blocked proxy (e.g., HepaRG, HEK293T)
        return re.sub(r"[^a-zA-Z0-9]", "", proxy).upper()

    # Unknown name: fall back to a stripped tag built from the raw input.
    return re.sub(r"[^a-zA-Z0-9]", "", raw).upper()


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
