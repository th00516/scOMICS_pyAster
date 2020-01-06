#!/usr/bin/env python3


import sys
import os.path
import gzip


class kmer:
    def __init__(self):
        """"""
        self.seq_mat = []
        self.seq_len = 0
        self.kmer_set = set()
        self.kmer_motif_set = {}

        self.k_typ = {'A': 1, 'C': 2, 'G': 4, 'T': 8}

    def reshape_reads(self, read_file):
        """"""
        with gzip.open(read_file, 'rt') as FIN:
            for _ in FIN:
                line = FIN.readline().strip()

                self.seq_len = len(line)
                line = line[5:self.seq_len - 5]
                self.seq_len = len(line)

                self.seq_mat.append(line)

                FIN.readline()
                FIN.readline()

    def kmer_motif_stat(self, K):
        """"""
        ###################################
        # only one copy would be assessed #
        ###################################
        self.kmer_set.update([seq[n:n + K] for n in range(0, self.seq_len - K + 1)
                              for seq in self.seq_mat if 'N' not in seq[n:n + K]])
        ###################################

        for k_seq in self.kmer_set:
            k_seq = [_ for _ in k_seq]
            k_sco = self.k_typ[k_seq[int(K / 2)]]

            k_seq[int(K / 2)] = 'N'
            k_seq = ''.join(k_seq)

            if k_seq not in self.kmer_motif_set:
                self.kmer_motif_set.update({k_seq: k_sco})

            else:
                self.kmer_motif_set[k_seq] += k_sco


if __name__ == '__main__':
    k = kmer()
    k.reshape_reads(sys.argv[1])
    k.kmer_motif_stat(int(sys.argv[-1]))

    with gzip.open(sys.argv[1] + '.snp_idx.gz', 'wt') as OU:
        for i in sorted(k.kmer_motif_set.keys()):
            x = k.kmer_motif_set[i]

            ###############################
            # count number of 1 in binary #
            ###############################
            count = 0
            while x:
                x = x & (x - 1)
                count += 1
            ###############################

            #######################################
            # only dualistic SNP should be stored #
            #######################################
            if count <= 2:
                print(i + '\t' + os.path.basename(sys.argv[1]) + '\t' + str(k.kmer_motif_set[i]), file=OU)
            #######################################
