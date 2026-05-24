from enum import Enum


class Interaction(Enum):
    RNA_RNA = "RNA_RNA"
    RNA_DNA_NO_WOBBLE = "RNA_DNA_NO_WOBBLE"
    DNA_DNA = "DNA_DNA"
    MODIFIED = "MODIFIED"
