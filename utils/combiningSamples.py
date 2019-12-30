#!/usr/bin/env python3


import sys
import pandas as pd
import numpy as np


def hpf(kmer_win, PoH_pass):
    """"""
    row, col = kmer_win.shape

    masker_window = np.ones((row, col), dtype=np.bool_)
    masker_window[int(row / 2) - int(row * (1 - PoH_pass)):int(row / 2) + int(row * (1 - PoH_pass)),
                  int(col / 2) - int(col * (1 - PoH_pass)):int(col / 2) + int(col * (1 - PoH_pass))] = 0

    kmer_win = abs(np.fft.ifft2(np.fft.ifftshift(np.fft.fftshift(np.fft.fft2(kmer_win)) * masker_window)))

    return kmer_win


if __name__ == '__main__':
    img_arr = []

    for mat in sys.argv[1:]:
        mat = np.load(mat)

        mat = hpf(mat, 0.6)

        mat = mat.reshape((1, mat.size))

        img_arr.extend(mat)

    D = pd.DataFrame(np.asarray(img_arr).transpose())
    D.to_csv('all_samples.csv', encoding='utf-8')
