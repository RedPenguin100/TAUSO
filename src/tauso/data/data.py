import os

import gffutils
from pyfaidx import Fasta
from platformdirs import user_data_dir

APP_NAME = "tauso"

def get_data_dir():
    # 1. Check for Docker/Render Environment Variable override
    env_path = os.environ.get('TAUSO_DATA_DIR')

    if env_path:
        d = env_path
    else:
        # 2. Fallback to standard User Data Dir (e.g., ~/.local/share/tauso)
        d = user_data_dir(APP_NAME)

    # Ensure the directory exists (works for both paths)
    os.makedirs(d, exist_ok=True)
    return d

def get_paths(genome="GRCh38"):
    """
    Returns paths for a specific genome assembly.
    Default: GRCh38
    """
    d = get_data_dir()
    return {
        "dir": d,
        "fasta": os.path.join(d, f"{genome}.fa"),
        "gtf": os.path.join(d, f"{genome}.gtf"),
        "db": os.path.join(d, f"{genome}.db")
    }

def load_db(genome="GRCh38"):
    """Returns the gffutils database for the specified genome."""
    paths = get_paths(genome)
    if not os.path.exists(paths['db']):
        raise FileNotFoundError(f"Database for {genome} not found. Run 'tauso setup-genome --genome {genome}'")
    return gffutils.FeatureDB(paths['db'])

def load_genome(genome="GRCh38"):
    """Returns the pyfaidx Fasta object."""
    paths = get_paths(genome)
    if not os.path.exists(paths['fasta']):
        raise FileNotFoundError(f"FASTA for {genome} not found. Run 'tauso setup-genome --genome {genome}'")
    return Fasta(paths['fasta'])