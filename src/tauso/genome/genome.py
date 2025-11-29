import os
import gffutils
from pyfaidx import Fasta
from platformdirs import user_data_dir
from ..consts import APP_NAME


def get_paths():
    """Returns a dict of paths. Validates existence."""
    d = user_data_dir(APP_NAME)
    paths = {
        "fasta": os.path.join(d, "GRCh38.fa"),
        "gtf": os.path.join(d, "GRCh38.gtf"),
        "db": os.path.join(d, "GRCh38.db")
    }

    # Simple check
    if not os.path.exists(paths['db']):
        raise FileNotFoundError(
            "Tauso data not found. Please run: 'tauso setup-genome'"
        )
    return paths


def load_genome():
    """Returns a pyfaidx.Fasta object (Acts just like genomepy.Genome)."""
    paths = get_paths()
    # Fasta() is efficient and lazy-loading
    return Fasta(paths['fasta'])


def load_db():
    """Returns the gffutils database."""
    paths = get_paths()
    return gffutils.FeatureDB(paths['db'])