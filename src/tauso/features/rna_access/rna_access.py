import os
import shlex
import subprocess
import uuid
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..._raccess.core import find_raccess
from .file_util import FileUtil


def parse_line(line_str, d_empty):
    pos, multi_data = line_str.split('\t')
    pos = int(pos)
    multi_data_list = list(filter(len, multi_data.split(';')))
    d = dict(s.split(',') for s in multi_data_list)

    d = {int(k): float(v) for k, v in d.items()}

    # removed not needed DataFrame already added Nan where no value
    # d = {k: d.get(k, float('nan')) for k in d_empty}

    return pos, d


def parse_single(data, segment_sizes):
    d_empty = {k: '' for k in segment_sizes}

    lines = data.splitlines()

    id_str = lines[0].rstrip()

    ind_rec_list = list(map(lambda line: parse_line(line, d_empty), filter(len, lines[1:])))

    # df = pd.DataFrame(info_list, index=, columns=segment_sizes)
    indexes = list(zip(*ind_rec_list))[0]
    records = list(zip(*ind_rec_list))[1]
    df = pd.DataFrame(records, index=indexes)

    return id_str, df


def parse(data, segment_sizes):
    seq_res_list = list(filter(len, data.split('>')))
    res = map(lambda seq_red: parse_single(seq_red, segment_sizes), seq_res_list)
    return dict(res)


# TODO review usage in parallel since it is writing temporal file needs uuid prefix for file and maybe delete it after
#  usage
class RNAAccess(object):
    # column in dataframe are segment sizes
    # pos 0-based coordinate
    # segment [pos, pos + segment_size - 1]
    USED_RT = 0.61633008  # [kcal/mol]

    def __init__(self, segment_sizes=None, max_span=None, exe_path=None):
        self.segment_sizes = segment_sizes
        self.max_span = max_span

        self.uuid_str = None
        if exe_path is None:
            self.exe_path = find_raccess()
        else:
            self.exe_path = exe_path

        out_dir = FileUtil.get_output_dir()
        p = Path(out_dir).expanduser().resolve()
        p.mkdir(parents=True, exist_ok=True)

    def set_uuid_for_web(self, uuid_str):
        self.uuid_str = uuid_str

    @staticmethod
    def to_seq_rec(seq_id):
        seq_rec = SeqRecord(
            Seq(seq_id[1]),
            id=seq_id[0],
            name=seq_id[0],
            description="",
        )
        return seq_rec

    def calculate(self, seq_id_list):
        if not self.uuid_str:
            uuid_str = str(uuid.uuid4().hex)
        else:
            uuid_str = self.uuid_str

        seq_file_name = 'trig_seq.fa'
        seq_file_name = uuid_str + '_' + seq_file_name
        seq_path = FileUtil.get_output_path(seq_file_name)

        out_file_name = 'trig_raccess.txt'
        out_file_name = uuid_str + '_' + out_file_name
        out_path = FileUtil.get_output_path(out_file_name)

        try:
            seq_rec_list = map(self.to_seq_rec, seq_id_list)
            SeqIO.write(seq_rec_list, seq_path, "fasta")

            segment_sizes_str = ','.join(map(str, self.segment_sizes))
            cmd = (f"{self.exe_path} -outfile={out_path} -seqfile={seq_path} "
                   f"-access_len={segment_sizes_str} -max_span={self.max_span}")

            command = shlex.split(cmd)
            p = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

            with open(out_path, 'r') as f:
                data = f.read()
            res = parse(data, self.segment_sizes)
            return res
        finally:
            if os.path.exists(seq_path):
                os.remove(seq_path)
            if os.path.exists(out_path):
                os.remove(out_path)


if __name__ == "__main__":
    # FileUtil.set_root_dir()

    YEAST_M_CHERRY = 'ATGTCTAAGGGGGAAGAAGACAATATGGCGATTATTAAAGAGTTTATGAGATTTAAAGTACATATGGAAGGAAGTGTTAATGGTCACGAGTTTGAGATCGAAGGTGAAGGTGAAGGTCGTCCATATGAGGGTACGCAAACAGCAAAACTAAAGGTGACTAAAGGGGGACCATTACCTTTCGCTTGGGATATACTGTCACCACAATTCATGTACGGATCGAAAGCTTACGTAAAGCACCCGGCCGACATTCCTGATTATTTAAAGTTGTCTTTCCCTGAAGGGTTCAAATGGGAAAGAGTTATGAATTTTGAGGATGGAGGTGTTGTGACGGTAACTCAAGATTCATCTTTGCAAGATGGCGAATTCATTTATAAAGTTAAATTGAGAGGAACTAACTTTCCAAGCGATGGTCCAGTCATGCAAAAAAAGACCATGGGCTGGGAAGCTAGCTCAGAACGGATGTACCCGGAAGACGGCGCATTAAAGGGAGAGATCAAGCAGCGACTTAAGTTAAAAGATGGCGGGCATTATGATGCAGAAGTAAAGACAACCTACAAAGCCAAAAAACCCGTGCAGCTGCCTGGTGCGTATAATGTTAACATAAAACTAGACATTACATCCCACAACGAAGACTACACTATAGTCGAACAATACGAAAGGGCAGAAGGTAGACATTCGACAGGTGGTATGGATGAGTTGTATAAATAA'.replace(
        'T', 'U')
    g_seq = YEAST_M_CHERRY

    g_ra = RNAAccess([6], 120)
    # g_seq_id_list = [('trigger', g_seq[:60]), ('trigger2', g_seq[:60])]
    g_seq_id_list = [('trigger', g_seq)]
    g_res = g_ra.calculate(g_seq_id_list)
    print(g_res)
