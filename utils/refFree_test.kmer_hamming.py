#!/usr/bin/env python3


import sys
import gzip
import pandas as pd
import numpy as np

from itertools import product


class kmer_encoder:
    def __init__(self):
        """"""
        self.read_seq_mat = []
        self.read_seq_len = 0

        self.kmer_win = np.array([], dtype=np.float32)

    def file_buffer(self, read_file):
        """"""
        with gzip.open(read_file, 'rt') as FIN:
            for _ in FIN:
                line = FIN.readline().strip()

                self.read_seq_len = len(line)
                line = line[5:self.read_seq_len - 5]
                self.read_seq_len = len(line)

                self.read_seq_mat.append(line)

                FIN.readline()
                FIN.readline()

    def kmer_hamming(self):
        """"""
        pass


if __name__ == '__main__':
    pass
