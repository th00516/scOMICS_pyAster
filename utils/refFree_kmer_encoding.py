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
        if K % 2 == 0:
            print('K should be an odd', file=sys.stderr)
            exit(1)

        total_kmer_freq = []

        for i in range(0, self.read_seq_len - (K - 1)):
            kmer_freqs = {''.join(_): 0 for _ in product(('A', 'T', 'C', 'G'), repeat=K)}

            kmer_counts = dict(pd.value_counts([_[i:i + K] for _ in self.read_seq_mat], sort=False))

            for kmer_seq in kmer_counts.keys():
                if 'N' not in kmer_seq:
                    kmer_freqs[kmer_seq] = kmer_counts[kmer_seq]

            kmer = [kmer_freqs[_] for _ in sorted(kmer_freqs.keys())]

            total_kmer_freq.append(kmer)

        self.kmer_win = np.asarray(total_kmer_freq)

    def encode_kmer(self, K):
        """"""
        if K % 2 == 0:
            print('K should be an odd', file=sys.stderr)
            exit(1)

        total_kmer_freq = []

        kmer_freqs = {''.join(_): 0 for _ in product(('A', 'T', 'C', 'G'), repeat=K)}

        for i in range(0, self.read_seq_len - (K - 1)):
            kmer_counts = dict(pd.value_counts([_[i:i + K] for _ in self.read_seq_mat], sort=False))

            for kmer_seq in kmer_counts.keys():
                if 'N' not in kmer_seq:
                    kmer_freqs[kmer_seq] += kmer_counts[kmer_seq]

        kmer = [kmer_freqs[_] for _ in sorted(kmer_freqs.keys())]

        total_kmer_freq.append(kmer)

        self.kmer_win = np.asarray(total_kmer_freq)

    def kmer_filtering(self, f_min, f_max):
        """"""
        self.kmer_win[self.kmer_win < f_min] = 0
        self.kmer_win[self.kmer_win > f_max] = 0


if __name__ == '__main__':
    encoder = kmer_encoder()
    encoder.file_buffer(sys.argv[2])

    if sys.argv[1] == '--site_kmer':
        encoder.encode_site_kmer(7)

    elif sys.argv[1] == '--kmer':
        encoder.encode_kmer(9)

    else:
        print('Invalid option', file=sys.stderr)

    encoder.kmer_filtering(20, 1000)

    np.save(sys.argv[2] + '.npy', encoder.kmer_win)
