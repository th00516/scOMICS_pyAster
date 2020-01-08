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
        self.kmer_motif_type_set = {}
        self.kmer_motif_set = {}

        self.k_type_maker = {'A': 1, 'C': 2, 'G': 4, 'T': 8}

    def reshape_reads(self, read_file):
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

    def kmer_motif_stat(self, K):
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
        # keep k-mer that locate in 90%~95% of its ordered frequency #
        ##############################################################
        F_min, F_max = np.percentile([_ for _ in self.kmer_set.values()], [90, 95])
        ##############################################################

        for k_seq in self.kmer_set.keys():
            if F_min <= self.kmer_set[k_seq] <= F_max:
                k_seq_arr = [_ for _ in k_seq]

                k_seqL = ''.join(k_seq_arr[:int(int(sys.argv[-1]) / 2)])
                k_seqR = ''.join(k_seq_arr[int(int(sys.argv[-1]) / 2) + 1:])

                k_seq_motif = ''.join((k_seqL, 'N', k_seqR))

                k_mak = self.k_type_maker[k_seq_arr[int(K / 2)]]

                if k_seq_motif not in self.kmer_motif_type_set:
                    self.kmer_motif_type_set.update({k_seq_motif: k_mak})
                    self.kmer_motif_set.update({k_seq_motif: self.kmer_set[k_seq]})

                else:
                    self.kmer_motif_type_set[k_seq_motif] += k_mak
                    self.kmer_motif_set[k_seq_motif] += self.kmer_set[k_seq]


if __name__ == '__main__':
    k = kmer()
    k.reshape_reads(sys.argv[1])
    k.kmer_motif_stat(int(sys.argv[-1]))

    T_OU = open(sys.argv[1] + '.snp_idx.type_stat', 'wt')
    F_OU = open(sys.argv[1] + '.snp_idx.freq_stat', 'wt')

    for i in sorted(k.kmer_motif_type_set.keys()):
        x = k.kmer_motif_type_set[i]

        ###############################
        # count number of 1 in binary #
        ###############################
        count = 0
        while x:
            x = x & (x - 1)
            count += 1
        ###############################

        ####################################################
        # only homozygous & dualistic SNP should be stored #
        ####################################################
        if count <= 2:
            i = [_ for _ in i]

            iL = ''.join(i[:int(int(sys.argv[-1]) / 2)])
            iR = ''.join(i[int(int(sys.argv[-1]) / 2) + 1:])

            score_t = k.kmer_motif_type_set[''.join((iL, 'N', iR))]
            score_f = k.kmer_motif_set[''.join((iL, 'N', iR))]

            print(''.join((iL, 'N', iR)) + ',' + os.path.basename(sys.argv[1]) + ',' + str(score_t), file=T_OU)
            print(''.join((iL, 'N', iR)) + ',' + os.path.basename(sys.argv[1]) + ',' + str(score_f), file=F_OU)
        ####################################################

    T_OU.close()
    F_OU.close()
