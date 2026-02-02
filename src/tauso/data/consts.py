import re



SEQUENCE = 'Sequence'
INHIBITION = 'Inhibition(%)'
CANONICAL_GENE = 'Canonical Gene Name'
CELL_LINE_ORGANISM = 'Cell line organism'
VOLUME = 'ASO_volume(nM)'
CHEMICAL_PATTERN = 'Chemical_Pattern'
TREATMENT_PERIOD = 'Treatment_Period(hours)'
CELL_LINE = 'Cell_line'
TRANSFECTION = 'Transfection'
DENSITY = 'Density(cells/well)'
DENSITY_UPDATED = 'Density(cells_per_well)' # Avoiding /
PRIMER_PROBE = 'Primer_probe_set'
MODIFICATION = 'Modification'
SENSE_START = 'sense_start'
SENSE_LENGTH = 'sense_length'
SENSE_TYPE = 'sense_type'
SENSE_SEQUENCE = 'sense_sequence'
PRE_MRNA_SEQUENCE = 'pre_mrna_sequence'
SENSE_AVG_ACCESSIBILITY = 'sense_avg_accessibility'

CELL_LINE_DEPMAP = 'Cell_Line_Depmap'

ISIS = 'ISIS'
LINKAGE = 'Linkage'
LINKAGE_LOCATION = 'Linkage_Location'
TARGET_GENE = 'Target_gene'
SMILES = 'Smiles'


def standardize_cell_line_name(name: str) -> str:
    """
    Standardizes cell line names to alphanumeric uppercase.
    Example: 'A-431' -> 'A431'
    """
    if not name or not isinstance(name, str):
        return str(name)
    if name == "U-251" or name == "U251":
        return "U251MG"
    return re.sub(r'[^a-zA-Z0-9]', '', name).upper()

HEPG2        = 'HepG2'      # Note: DepMap usually lists this as HepG2, but keeping your key 'HepG'
SNU_449     = 'SNU-449'
HELA        = 'HeLa'
A431        = 'A431'
SK_MEL_28   = 'SK-MEL-28'
SH_SY5Y     = 'SH-SY5Y'
U251        = 'U251'
U_251_MG    = 'U-251 MG'
H929        = 'H929'
KMS11       = 'KMS11'
NCI_H460    = 'NCI-H460'
SK_N_AS     = 'SK-N-AS'
SK_N_SH     = 'SK-N-SH'

CELL_LINE_TO_DEPMAP = {
    'HepG2': 'ACH-000739',
    'SNU-449': 'ACH-000420',
    'HeLa': 'ACH-001086',
    'A431': 'ACH-001328',
    'A-431': 'ACH-001328',
    'SK-MEL-28': 'ACH-000615',
    'SH-SY5Y': 'ACH-001188',
    'U251': 'ACH-000232',
    'U-251 MG': 'ACH-000232',
    'H929': 'ACH-000050',
    'KMS11': 'ACH-000714',
    'NCI-H460': 'ACH-000463',
    'SK-N-AS': 'ACH-000260',
    'SK-N-SH': 'ACH-000149',
}

CELL_LINE_MAPPING = {
    # Direct Matches & Synonyms
    'HepG2':     HEPG2,
    'HepaRG':    HEPG2,       # Liver -> Liver proxy
    'SNU-449':   SNU_449,
    'Hela':      HELA,
    'A431':      A431,
    'A-431':      A431,
    'SK-MEL-28': SK_MEL_28,
    'SH-SY5Y':   SH_SY5Y,
    'U251':      U251,
    'U-251 MG':  U251,       # Map both to the canonical U251
    'H929':      H929,
    'KMS11':     KMS11,
    'NCI-H460':  NCI_H460,
    'SK-N-AS':   SK_N_AS,
    'SK-N-SH':   SK_N_SH,

    # Neuronal / Experimental Contexts
    'Angptl2/Actin':        SK_N_AS,
    'SK cells asyn':        SK_N_AS,  # "SK cells" usually refers to SK-N-AS in this context
    'PAC neurons asyn':     SH_SY5Y,  # SH-SY5Y is the standard neuronal proxy here
    'Human Neuronal Cell':  SH_SY5Y,  # Generic neuronal -> SH-SY5Y proxy

    # Blood / Immune Proxies
    'MM.1R':                H929,     # Myeloma -> H929 (Myeloma) proxy
}

HUMAN = 'human'
MONKEY = 'monkey'
MOUSE = 'mouse'
UNKNOWN = 'unknown'


cell_line_to_organism = {
    # --- Human ---
    'A431': HUMAN,
    'A-431' : HUMAN,
    'Angptl2/Actin': HUMAN,
    'CC-2580': HUMAN,
    'H929': HUMAN,
    'Hela': HUMAN,
    'HepG2': HUMAN,
    'HepaRG': HUMAN,
    'Human IPS': HUMAN,
    'Human Neuronal Cell': HUMAN,
    'KARPAS-229': HUMAN,
    'KMS11': HUMAN,
    'MM.1R': HUMAN,
    'NCI-H460': HUMAN,
    'PAC neurons asyn': HUMAN,
    'SH-SY5Y': HUMAN,
    'SK cells asyn': HUMAN,
    'SK-MEL-28': HUMAN,
    'SNU-449': HUMAN,
    'U251': HUMAN,

    # --- Mouse ---
    'Mouse Primary Neuronal Cell': MOUSE,
    'Oscillation': MOUSE,
    'Tubulin': MOUSE,

    # --- Monkey ---
    'LLC-MK2 monkey': MONKEY, # Green Monkey
}