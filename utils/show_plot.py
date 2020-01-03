#!/usr/bin/env python3


import sys
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages


if __name__ == '__main__':
    mat = np.load(sys.argv[1]).reshape((1, -1))
    mat = mat / np.max(mat) * 255

    pdf = PdfPages('eva/' + sys.argv[1] + '.pdf')

    P = plt.figure()

    F = P.add_subplot(111)
    F.scatter(range(0, np.nonzero(mat[0])[0].size), mat[0][np.nonzero(mat[0])[0]])

    plt.xlabel('kmer type')
    plt.ylabel('kmer freq')

    P.show()

    pdf.savefig()
    pdf.close()
