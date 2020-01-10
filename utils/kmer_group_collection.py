#!/usr/bin/env python3


import sys
import os.path
import gzip
import numpy as np


class kmer:
    def __init__(self):
        """"""
        self.seq_mat = []
        self.seq_len = 0
        self.kmer_set = {}
        self.kmer_group_set = set()

    def pack_reads(self, read_file):
        """"""
        with gzip.open(read_file, 'rt') as FIN:
            N = 1
            for _ in FIN:
                if N == 100000:
                    break

                line = FIN.readline().strip()

                self.seq_len = len(line)
                line = line[5:self.seq_len - 5]
                self.seq_len = len(line)

                self.seq_mat.append(line)

                FIN.readline()
                FIN.readline()

                N += 1

    def kmer_group_stat(self, K):
        """"""
        for seq in self.seq_mat:
            for n in range(0, self.seq_len - K + 1):
                if 'N' in seq[n:n + K]:
                    continue

                if seq[n:n + K] not in self.kmer_set:
                    self.kmer_set.update({seq[n:n + K]: 1})

                else:
                    self.kmer_set[seq[n:n + K]] += 1

        ##############################################################
        # keep k-mer that locate in 90%~98% of its ordered frequency #
        ##############################################################
        F_min, F_max = np.percentile([_ for _ in self.kmer_set.values()], [60, 95])
        ##############################################################

        for k_seq in self.kmer_set.keys()[::-1]:
            if F_min <= self.kmer_set[k_seq] <= F_max:
                if k_seq not in self.kmer_group_set:
                    self.kmer_group_set.add(k_seq)

                else:



if __name__ == '__main__':
    k = kmer()
    k.pack_reads(sys.argv[1])
    k.kmer_group_stat(int(sys.argv[-1]))
