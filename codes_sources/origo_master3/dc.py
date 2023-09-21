#!/usr/bin/env python
"""Implementation of Causal Inference on Discrete Data via Estimating
Distance Correlations.
Link: http://www.mitpressjournals.org/doi/pdf/10.1162/NECO_a_00820
"""
from collections import defaultdict
from math import sqrt
import numpy as np
from scipy.spatial.distance import pdist, squareform

if __name__ == "__main__":
    from utils import normalise, pad_sequence

else:
    from .utils import normalise, pad_sequence



def dc(X, Y):
    prob_X, marg_X, prob_Y, marg_Y = distributions(X, Y)
    dXtoY = dCor(prob_X, marg_Y)
    dYtoX = dCor(prob_Y, marg_X)
    return (dXtoY, dYtoX)


def dcov(X, Y):
    n = X.shape[0]
    XY = np.multiply(X, Y)
    cov = sqrt(XY.sum()) / n
    return cov


def dvar(X):
    return sqrt(np.sum(X ** 2 / X.shape[0] ** 2))


def cent_dist(X):
    M = squareform(pdist(X))    # distance matrix
    rmean = M.mean(axis=1)
    cmean = M.mean(axis=0)
    gmean = rmean.mean()
    Rmean = np.tile(rmean, (M.shape[0], 1)).transpose()
    Cmean = np.tile(cmean, (M.shape[1], 1))
    Gmean = np.tile(gmean, M.shape)
    CM = M - Rmean - Cmean + Gmean
    return CM


def dCor(X, Y):
    A = cent_dist(X)
    B = cent_dist(Y)

    dcov_AB = dcov(A, B)
    dvar_A = dvar(A)
    dvar_B = dvar(B)

    dcor = 0.0
    if dvar_A > 0.0 and dvar_B > 0.0:
        dcor = dcov_AB / sqrt(dvar_A * dvar_B)

    return dcor


def distributions(X, Y):
    N = len(X)
    unq_X = set(map(tuple, X))
    unq_Y = set(map(tuple, Y))
    idx = list(range(N))
    idx_X = dict(list(zip(unq_X, idx)))
    idx_Y = dict(list(zip(unq_Y, idx)))

    freq_XY = np.zeros((len(unq_X), len(unq_Y)))
    for i in range(N):
        ix = idx_X[tuple(X[i])]
        iy = idx_Y[tuple(Y[i])]
        freq_XY[ix, iy] += 1

    freq_X = np.sum(freq_XY, axis=1)[np.newaxis]
    freq_Y = np.sum(freq_XY, axis=0)[np.newaxis]
    prob_X = (freq_X / np.sum(freq_X)).transpose()
    prob_Y = (freq_Y / np.sum(freq_Y)).transpose()

    freqs_X = np.tile(freq_X.transpose(), (1, len(unq_Y)))
    freqs_Y = np.tile(freq_Y, (len(unq_X), 1))
    marg_X = (freq_XY / freqs_X).transpose()
    marg_Y = freq_XY / freqs_Y
    return prob_X, marg_X, prob_Y, marg_Y


if __name__ == "__main__":
    X = [[2, 3], [2, 3], [2, 4], [2], [2], [3], [3], [3, 4], [2, 3]]
    Y = [[1], [1], [1], [1], [0], [0], [0], [0], [1]]

    print(dc(X, Y))
