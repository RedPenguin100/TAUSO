"""Build the OligoAI training data: download -> split -> canonical gene -> index -> process.

Orchestrates the numbered steps in this folder. Outputs go to raw_data/ (gitignored): the feature
pipeline reads the indexed table, the model reads the averaged table. Step 3 needs the genome
(`tauso setup-genome`); pass --skip-process to stop after indexing.

    python notebooks/data/OligoAI/setup_data.py [--cpus N] [--skip-process]
"""
import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = [
    "1_download.py",
    "1_5_assign_split.py",
    "2_assign_canonical_gene.py",
    "2_5_index.py",
    "3_process_data.py",
]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cpus", type=int, default=1, help="CPUs for the structure step (step 3).")
    p.add_argument("--skip-process", action="store_true", help="Stop after indexing (skip step 3).")
    args = p.parse_args()

    steps = STEPS[:-1] if args.skip_process else STEPS
    for name in steps:
        cmd = [sys.executable, str(HERE / name)]
        if name == "3_process_data.py":
            cmd += ["--cpus", str(args.cpus)]
        print(f"--- {name} ---", flush=True)
        subprocess.run(cmd, check=True)
    print("Data ready under raw_data/")


if __name__ == "__main__":
    main()
