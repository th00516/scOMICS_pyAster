#!/usr/bin/env python3


import sys
import scipy.spatial as ss
import numpy as np

from itertools import product


class exp_matrix:
    def __init__(self, vec, K):
        """"""
        self.vec = vec
        self.K = K
        self.euc_mat = np.array([])

    def binary_expression(self):
        """"""
        self.vec[self.vec > 0] = 1

    def generate_euclidean_matrix_model(self):
        """"""
        kmer = [''.join(_) for _ in product(('A', 'T', 'C', 'G'), repeat=4)]
        self.euc_mat = [_ for _ in product(kmer, repeat=2)]
