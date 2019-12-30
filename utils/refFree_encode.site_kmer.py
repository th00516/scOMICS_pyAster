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

    def encode_site_kmer(self, K):
        """"""
        total_kmer_freq = []

        for i in range(0, self.read_seq_len - (K - 1)):
            kmer_freqs = {''.join(_): 0 for _ in product(('A', 'T', 'C', 'G'), repeat=K)}

            kmer_counts = dict(pd.value_counts([_[i:i + K] for _ in self.read_seq_mat], sort=False))

            total_kmer_count = sum(kmer_counts.values())

            for kmer_seq in kmer_freqs.keys():
                if kmer_seq in kmer_counts:
                    kmer_freqs[kmer_seq] = kmer_counts[kmer_seq] / total_kmer_count

            kmer = [kmer_freqs[_] for _ in sorted(kmer_freqs.keys())]

            total_kmer_freq.append(kmer)

        self.kmer_win = np.asarray(total_kmer_freq)

    def encode_kmer(self, K):
        """"""
        total_kmer_freq = []

        kmer_freqs = {''.join(_): 0 for _ in product(('A', 'T', 'C', 'G'), repeat=K)}

        for i in range(0, self.read_seq_len - (K - 1)):
            kmer_counts = dict(pd.value_counts([_[i:i + K] for _ in self.read_seq_mat], sort=False))

            total_kmer_count = sum(kmer_counts.values())

            for kmer_seq in kmer_freqs.keys():
                if kmer_seq in kmer_counts:
                    kmer_freqs[kmer_seq] = kmer_counts[kmer_seq] / total_kmer_count

            kmer = [kmer_freqs[_] for _ in sorted(kmer_freqs.keys())]

            total_kmer_freq.append(kmer)

        self.kmer_win = np.asarray(total_kmer_freq)


if __name__ == '__main__':
    encoder = kmer_encoder()
    encoder.file_buffer(sys.argv[1])
    encoder.encode_site_kmer(6)

    np.save(sys.argv[1] + '.npy', np.around(encoder.kmer_win / np.max(encoder.kmer_win) * 255))
