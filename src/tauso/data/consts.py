import re

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
    'HepG': 'ACH-000739',
    'SNU-449': 'ACH-000420',
    'HeLa': 'ACH-001086',
    'A431': 'ACH-001328',
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