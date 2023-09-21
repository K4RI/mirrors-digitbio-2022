#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Here we test causal direction of multivariate tübingen cause-effect pairs.
"""
import itertools
import os
import numpy as np

from univariate_discretizer import UnivariateIPDiscretizer
from multivariate_discretizer import *
from score import score_ergo_origo_dc, score_ergo_origo_dc_old
from store import Store
from utils import concatenate, normalise


def pair0071():
    tubingen_pairs_dir = os.path.join(os.path.dirname(
        __file__), "data", "tubingen-pairs")
    pair_path = os.path.join(tubingen_pairs_dir, "pair0071.txt")
    data = np.loadtxt(pair_path)
    dim_X, dim_Y = (6, 2)
    X = [data[:, i].tolist() for i in range(dim_X)]
    Y = [data[:, dim_X + i].tolist() for i in range(dim_Y)]

    discretizer = UnivariateIPDiscretizer(
        X[0], list(range(len(X[0]))), num_bins_X=2, num_bins_Y=2)
    ax, xd, ay, yd = discretizer.discretize()

    Xd = [xd]
    aX = [ax]

    next_alphabet = max(ax) + 1
    for data in X[1:]:
        data = [next_alphabet if val == 1 else '' for val in data]
        Xd.append(data)
        aX.append([next_alphabet])
        next_alphabet += 1

    Yd = []
    aY = []
    for data in Y:
        data = [next_alphabet if val == 1 else '' for val in data]
        Yd.append(data)
        aY.append([next_alphabet])
        next_alphabet += 1

    aX = list(itertools.chain(*aX))
    aY = list(itertools.chain(*aY))

    rows = concatenate(Xd, Yd)
    store = Store(None, rows, aX + aY)

    from dc import dc
    print("distance corr ::", dc(store.project(aX), store.project(aY)))

    print("pair     = 71")
    print("dim_X    = %d" % dim_X)
    print("dim_Y    = %d" % dim_Y)
    ergo, origo, dc = score_ergo_origo_dc(aX, aY, store)
    print("truth    = X → Y")

    if ergo.X_to_Y < ergo.Y_to_X:
        print("ergo     = X → Y")
    elif ergo.X_to_Y == ergo.Y_to_X:
        print("ergo     = X ~ Y")
    else:
        print("ergo     = Y → X")

    if origo.X_to_Y < origo.Y_to_X:
        print("origo    = X → Y")
    elif origo.X_to_Y == origo.Y_to_X:
        print("origo    = X ~ Y")
    else:
        print("origo    = Y → X")


pair0071()

truths = ["Y", "Y", "X", "Y", "X", "X", "X"]
dims = [(4, 4), (1, 3), (3, 2), (16, 16), (3, 7),
        (4, 5), (6, 1)]  # chemnitz, stock
pair_nums = [52, 53, 54, 55, 0, 100, 101]  # chemnitz, stock, car
tubingen_pairs_dir = os.path.join(
    os.path.dirname(__file__), "data", "tubingen-pairs")

for i, pair_num in enumerate(pair_nums):
    pair_path = os.path.join(tubingen_pairs_dir, "pair00%d.txt" % pair_num)
    data = np.loadtxt(pair_path)

    dim_X, dim_Y = dims[i]
    # do sth else with pair0071
    print()
    print("pair     = %d" % pair_num)
    print("data_dim =", data.shape)
    print("dim_X    = %d" % dim_X)
    print("dim_Y    = %d" % dim_Y)

    X = [data[:, j].tolist() for j in range(dim_X)]
    Y = [data[:, dim_X + j].tolist() for j in range(dim_Y)]

    discretizer = IPMvDiscretizer(X, Y)
    aX, Xd, aY, Yd = discretizer.discretize()

    aX = list(itertools.chain(*aX))
    aY = list(itertools.chain(*aY))

    rows = concatenate(Xd, Yd)
    store = Store(None, rows, aX + aY)
    ergo, origo, dc = score_ergo_origo_dc(aX, aY, store)

    truth = truths[i]
    if truth == "X":
        print("truth    = X → Y")
    else:
        print("truth    = Y → X")

    if ergo.X_to_Y < ergo.Y_to_X:
        print("ergo     = X → Y")
    elif ergo.X_to_Y == ergo.Y_to_X:
        print("ergo     = X ~ Y")
    else:
        print("ergo     = Y → X")

    if origo.X_to_Y < origo.Y_to_X:
        print("origo    = X → Y")
    elif origo.X_to_Y == origo.Y_to_X:
        print("origo    = X ~ Y")
    else:
        print("origo    = Y → X")

    if dc.X_to_Y < dc.Y_to_X:
        print("dc       = X → Y")
    elif dc.X_to_Y == dc.Y_to_X:
        print("dc       = X ~ Y")
    else:
        print("dc       = Y → X")
