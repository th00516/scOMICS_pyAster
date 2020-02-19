#!/usr/bin/env python3


import sys
import os.path
import gzip
import pandas as pd
import numpy as np
import numba as nb


class kmer_cluster:
    def __init__(self):
        """"""
        self.seq_len = 0

        self.kmer_set = set()
        self.kmer_group = []
        self.kmer_group_mat = pd.DataFrame()

    def gen_highF_kmer(self, read_files, K):
        """"""
        kmer_box = {}

        for f in read_files:
            with gzip.open(f, 'rt') as FIN:
                N = 1
                for _ in FIN:
                    if N == 1000000:
                        break

                    seq = FIN.readline().strip()

                    self.seq_len = len(seq)
                    seq = seq[5:self.seq_len - 5]
                    self.seq_len = len(seq)

                    for n in range(0, self.seq_len - K + 1):
                        if 'N' in seq[n:n + K]:
                            continue

                        if seq[n:n + K] not in kmer_box:
                            kmer_box.update({seq[n:n + K]: 1})

                        else:
                            kmer_box[seq[n:n + K]] += 1

                    FIN.readline()
                    FIN.readline()

                    N += 1

        F_min, F_max = np.percentile([_ for _ in kmer_box.values()], [90, 95])

        self.kmer_set.update([_ for _ in kmer_box.keys() if F_min <= kmer_box[_] <= F_max])

    @nb.njit()
    def gen_group_mat(self):
        """"""
        kmer_group_mat = []

        for R in self.kmer_set:
            for C in self.kmer_set:
                if R[1:] == C[:-1]:
                    kmer_group_mat.append((R, C, 1))

                else:
                    kmer_group_mat.append((R, C, 0))

        self.kmer_group_mat = pd.DataFrame(np.asarray(kmer_group_mat))
        self.kmer_group_mat = self.kmer_group_mat.pivot(index=0, columns=1, values=2)


if __name__ == '__main__':
    kmer = kmer_cluster()
    kmer.gen_highF_kmer(sys.argv[1:], 7)
    kmer.gen_group_mat()

    kmer.kmer_group_mat.to_csv('debug.csv')
