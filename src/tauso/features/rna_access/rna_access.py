import os
import shlex
import subprocess
import tempfile
import uuid

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..._raccess.core import find_raccess


def _parse_position_line(line):
    """'pos<TAB>seed,energy;seed,energy;...' -> (pos, {seed_size: energy})."""
    pos, energies = line.split("\t")
    pairs = (entry.split(",") for entry in energies.split(";") if entry)
    return int(pos), {int(seed): float(energy) for seed, energy in pairs}


def _parse_record(block):
    """One '>'-delimited raccess record -> (record_id, per-position energy DataFrame)."""
    header, *lines = block.splitlines()
    positions, energies = zip(*(_parse_position_line(line) for line in lines if line))
    return header.rstrip(), pd.DataFrame(list(energies), index=list(positions))


def parse_raccess_output(text):
    """Parse raccess stdout into {record_id: DataFrame(index=pos, columns=seed sizes)}.

    Missing (pos, seed) pairs become NaN, courtesy of the DataFrame constructor.
    """
    return dict(_parse_record(block) for block in text.split(">") if block)


class RNAAccess:
    """Thin wrapper around the `run_raccess` executable.

    `calculate` returns one DataFrame per input sequence: columns are seed
    (segment) sizes, the index is the 0-based position, and the value at row
    `pos` is the opening energy of the segment [pos, pos + seed_size - 1].
    """

    USED_RT = 0.61633008  # [kcal/mol]

    def __init__(self, segment_sizes=None, max_span=None, exe_path=None):
        self.segment_sizes = segment_sizes
        self.max_span = max_span
        self.uuid_str = None
        self.exe_path = exe_path if exe_path else find_raccess()

    def set_uuid_for_web(self, uuid_str):
        self.uuid_str = uuid_str

    @staticmethod
    def to_seq_rec(seq_id):
        rna_id, seq = seq_id
        return SeqRecord(Seq(seq), id=rna_id, name=rna_id, description="")

    def calculate(self, seq_id_list):
        """Run raccess on (rna_id, seq) pairs; return {rna_id: per-position energy DataFrame}.

        Temp files live in /dev/shm (RAM) when available and are always cleaned up.
        """
        run_id = self.uuid_str or uuid.uuid4().hex
        ram_dir = "/dev/shm" if os.path.isdir("/dev/shm") else tempfile.gettempdir()
        seq_path = os.path.join(ram_dir, f"{run_id}_trig_seq.fa")
        out_path = os.path.join(ram_dir, f"{run_id}_trig_raccess.txt")

        try:
            SeqIO.write(map(self.to_seq_rec, seq_id_list), seq_path, "fasta")

            seeds = ",".join(map(str, self.segment_sizes))
            cmd = (
                f"{self.exe_path} -outfile={out_path} -seqfile={seq_path} "
                f"-access_len={seeds} -max_span={self.max_span}"
            )
            subprocess.run(
                shlex.split(cmd),
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )

            with open(out_path, "r") as f:
                return parse_raccess_output(f.read())
        finally:
            for path in (seq_path, out_path):
                if os.path.exists(path):
                    os.remove(path)


if __name__ == "__main__":
    YEAST_M_CHERRY = "ATGTCTAAGGGGGAAGAAGACAATATGGCGATTATTAAAGAGTTTATGAGATTTAAAGTACATATGGAAGGAAGTGTTAATGGTCACGAGTTTGAGATCGAAGGTGAAGGTGAAGGTCGTCCATATGAGGGTACGCAAACAGCAAAACTAAAGGTGACTAAAGGGGGACCATTACCTTTCGCTTGGGATATACTGTCACCACAATTCATGTACGGATCGAAAGCTTACGTAAAGCACCCGGCCGACATTCCTGATTATTTAAAGTTGTCTTTCCCTGAAGGGTTCAAATGGGAAAGAGTTATGAATTTTGAGGATGGAGGTGTTGTGACGGTAACTCAAGATTCATCTTTGCAAGATGGCGAATTCATTTATAAAGTTAAATTGAGAGGAACTAACTTTCCAAGCGATGGTCCAGTCATGCAAAAAAAGACCATGGGCTGGGAAGCTAGCTCAGAACGGATGTACCCGGAAGACGGCGCATTAAAGGGAGAGATCAAGCAGCGACTTAAGTTAAAAGATGGCGGGCATTATGATGCAGAAGTAAAGACAACCTACAAAGCCAAAAAACCCGTGCAGCTGCCTGGTGCGTATAATGTTAACATAAAACTAGACATTACATCCCACAACGAAGACTACACTATAGTCGAACAATACGAAAGGGCAGAAGGTAGACATTCGACAGGTGGTATGGATGAGTTGTATAAATAA".replace(
        "T", "U"
    )

    ra = RNAAccess([6], 120)
    res = ra.calculate([("trigger", YEAST_M_CHERRY)])
    print(res)
