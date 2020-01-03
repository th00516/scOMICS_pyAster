#!/usr/bin/env python3


import sys
import gzip
import multiprocessing


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
        self.kmer_set.update([seq[n:n + K] for n in range(0, self.seq_len - K + 1)
                              for seq in self.seq_mat if 'N' not in seq[n:n + K]])

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
    def kmer_handle(read_file, K):
        k = kmer()
        k.reshape_reads(read_file)
        k.kmer_motif_stat(int(K))

        return k.kmer_motif_set

    pool = multiprocessing.Pool(processes=2)
    results = []

    for i in range(0, 2):
        results.append(pool.apply_async(kmer_handle, (sys.argv[1], sys.argv[-1])))

    pool.close()
    pool.join()

    k1 = results[0].get()
    k2 = results[1].get()

    for i in sorted(set(list(k1.keys()) + list(k2.keys()))):
        k1[i] = k1[i] if i in k1 else 0
        k2[i] = k2[i] if i in k2 else 0

        print(i, k1[i], k2[i])
