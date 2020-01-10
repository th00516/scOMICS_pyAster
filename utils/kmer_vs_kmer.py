#!/usr/bin/env python3


# import sys
import numpy as np


class kmer_vs_kmer:
    def __init__(self):
        self.kmer_lst1 = np.array([])
        self.kmer_lst2 = np.array([])
        self.kmer = np.array([])

    def test(self, arr1, arr2):
        self.kmer_lst1 = np.array(arr1).reshape([-1, 1])
        self.kmer_lst2 = np.array(arr2)

        self.kmer = self.kmer_lst1 * self.kmer_lst2

        print(self.kmer)


if __name__ == '__main__':
    kmer_group = kmer_vs_kmer()
    kmer_group.test([1, 2, 3, 1, 5, 6], [2, 2, 6, 3, 4, 6])
