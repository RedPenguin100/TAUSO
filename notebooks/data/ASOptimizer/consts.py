"""ASOptimizer-specific column-name constants.

Parked here while the ASOptimizer pipeline lives in
``notebooks/data/ASOptimizer/`` as a decoupled module. These columns
originate from ASOptimizer source files and don't fit TAUSO's canonical
naming convention. When the ASOptimizer integration is reworked, these
constants either get folded into ``src/tauso/data/consts.py`` (canonical
TAUSO names) or stay here as ASOptimizer-side input names that get
renamed by an ASOptimizer rename map analogous to
``notebooks.preprocessing.process_oligo_data_rename``.
"""

SMILES = "Smiles"
LINKAGE_LOCATION = "Linkage_Location"
