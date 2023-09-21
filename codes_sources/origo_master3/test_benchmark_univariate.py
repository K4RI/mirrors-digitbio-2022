#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Here we test causal direction of univariate tübingen cause-effect pairs.
"""

import glob
import os
import re

import numpy as np
from scipy.stats import binom
import random

from univariate_discretizer import UnivariateIPDiscretizer, UnivariateNMLDiscretizer
from score import score_ergo_origo_dc
from store import Store
from utils import normalise, progress, reverse_argsort


XY_TRUTH_PATTERN = re.compile("x\s*-+\s*-*\s*>\s*y", re.IGNORECASE)
YX_TRUTH_PATTERN = re.compile("y\s*-+\s*-*\s*>\s*x", re.IGNORECASE)


def get_ground_truths_of_tubingen_pairs():
    tubingen_pairs_dir = os.path.join(os.path.dirname(
        __file__), "data", "tubingen-pairs")
    tubingen_pairs = glob.glob(
        tubingen_pairs_dir + os.sep + "pair[0-9][0-9][0-9][0-9]_des.txt")

    truths = []
    for pair in tubingen_pairs:
        with open(pair, "r") as fp:
            desc = fp.read()
            if XY_TRUTH_PATTERN.findall(desc):
                truths.append(("X", "Y"))
            else:
                truths.append(("Y", "X"))
                # assert len(YX_TRUTH_PATTERN.findall(desc)) > 0
    return truths


def get_all_tubingen_pairs():
    tubingen_pairs_dir = os.path.join(os.path.dirname(
        __file__), "data", "tubingen-pairs")
    tubingen_pairs = glob.glob(
        tubingen_pairs_dir + os.sep + "pair[0-9][0-9][0-9][0-9].txt")

    return tubingen_pairs


def load_tubingen_pair(pair_path):
    data = np.loadtxt(pair_path)
    X = data[:, 0]
    Y = data[:, 1]
    return X, Y


def load_tubingen_pairs():
    tubingen_pairs = get_all_tubingen_pairs()
    for pair in tubingen_pairs:
        yield load_tubingen_pair(pair)


def test_tubingen_pairs():
    epsilon = 0.0
    level = 0.05
    truths = get_ground_truths_of_tubingen_pairs()
    multivariate_pairs = [52, 53, 54, 55, 71]
    num_pairs = len(truths) - len(multivariate_pairs)

    num_correct = 0
    num_wrong = 0
    nsample = 0
    num_indecisive = 0
    res_ergo, res_origo, res_dc = [], [], []
    diffs_ergo, diffs_origo, diffs_dc = [], [], []

    progress(0, 95)
    for i, data in enumerate(load_tubingen_pairs()):
        if i + 1 in multivariate_pairs:
            continue

        X, Y = data
        # if i+1 in [65, 66, 67]:
        #     X, Y = preprocess_dc_style(X, 100), preprocess_dc_style(Y, 100)
        # else:
        #     X, Y = preprocess_dc_style(X), preprocess_dc_style(Y)
        discretizer = UnivariateIPDiscretizer(X, Y)
        aX, Xd, aY, Yd = discretizer.discretize()
        rows = np.column_stack((Xd, Yd)).tolist()
        store = Store(None, rows, aX + aY)
        ergo, origo, dc = score_ergo_origo_dc(aX, aY, store)

        nsample += 1
        diffs_ergo.append(abs(ergo.X_to_Y - ergo.Y_to_X))
        diffs_origo.append(abs(origo.X_to_Y - origo.Y_to_X))
        diffs_dc.append(abs(dc.X_to_Y - dc.Y_to_X))

        if ergo.X_to_Y < ergo.Y_to_X:
            cause_ergo = "X"
        elif ergo.X_to_Y > ergo.Y_to_X:
            cause_ergo = "Y"
        else:
            cause_ergo = ""

        if origo.X_to_Y < origo.Y_to_X:
            cause_origo = "X"
        elif origo.X_to_Y > origo.Y_to_X:
            cause_origo = "Y"
        else:
            cause_origo = ""

        if dc.X_to_Y < dc.Y_to_X:
            cause_dc = "X"
        elif dc.X_to_Y > dc.Y_to_X:
            cause_dc = "Y"
        else:
            cause_dc = ""

        true_cause = truths[i][0]
        if cause_ergo == "":
            res_ergo.append(random.choice([True, False]))
        elif cause_ergo == true_cause:
            res_ergo.append(True)
        else:
            res_ergo.append(False)

        if cause_origo == "":
            res_origo.append(random.choice([True, False]))
        elif cause_origo == true_cause:
            res_origo.append(True)
        else:
            res_origo.append(False)

        if cause_dc == "":
            res_dc.append(random.choice([True, False]))
        elif cause_dc == true_cause:
            res_dc.append(True)
        else:
            res_dc.append(False)

        progress(nsample, 95)

    print()
    indices_ergo = reverse_argsort(diffs_ergo)
    indices_origo = reverse_argsort(diffs_origo)
    indices_dc = reverse_argsort(diffs_dc)

    diffs_ergo = [diffs_ergo[i] for i in indices_ergo]
    diffs_origo = [diffs_origo[i] for i in indices_origo]
    diffs_dc = [diffs_dc[i] for i in indices_dc]

    res_ergo = [res_ergo[i] for i in indices_ergo]
    res_origo = [res_origo[i] for i in indices_origo]
    res_dc = [res_dc[i] for i in indices_dc]

    dec_rate = np.arange(0.01, 1.01, 0.01)
    accs_ergo, accs_origo, accs_dc = [], [], []
    fp = open("results/dec_rate_benchmark.dat", "w")
    for r in dec_rate:
        maxIdx = int(round(r * nsample, 0))
        rergo = res_ergo[:maxIdx]
        rorigo = res_origo[:maxIdx]
        rdc = res_dc[:maxIdx]

        b1, b2 = binom(maxIdx, 0.5).interval(0.95)
        b1 = b1 / maxIdx
        b2 = b2 / maxIdx
        accs_ergo.append(sum(rergo) / len(rergo))
        accs_origo.append(sum(rorigo) / len(rorigo))
        accs_dc.append(sum(rdc) / len(rdc))
        fp.write("%.2f %.2f %.2f %.2f %.2f %.2f\n" %
                 (r, sum(rergo) / len(rergo), sum(rorigo) /
                  len(rorigo), sum(rdc) / len(rdc), b1, b2))
    fp.close()


def run_on_pair(pair_num):
    assert 1 <= pair_num <= 100
    print("running on pair %d" % pair_num)
    tubingen_pairs = get_all_tubingen_pairs()
    pair = tubingen_pairs[pair_num - 1]
    print("pair path is %s" % pair)
    X, Y = load_tubingen_pair(pair)

    X = normalise(X)
    Y = normalise(Y)

    epsilon_X = float(input("epsilon_X = "))
    delta_X = float(input("delta_X = "))

    epsilon_Y = float(input("epsilon_Y = "))
    delta_Y = float(input("delta_Y = "))

    discretizer = UnivariateNMLDiscretizer(
        X, Y, max_bins_X=5, max_bins_Y=5, epsilon_X=epsilon_X, epsilon_Y=epsilon_Y, delta_X=delta_X, delta_Y=delta_Y)
    aX, Xd, aY, Yd = discretizer.discretize()
    rows = np.column_stack((Xd, Yd)).tolist()
    store = Store(None, rows, aX + aY)
    ergo, origo, dc = score_ergo_origo_dc(aX, aY, store)

    truths = get_ground_truths_of_tubingen_pairs()
    truth = truths[pair_num - 1]
    if truth[0] == "X":
        correct_ergo = ergo.X_to_Y < ergo.Y_to_X
        correct_origo = origo.X_to_Y < origo.Y_to_X
    else:
        correct_ergo = ergo.X_to_Y > ergo.Y_to_X
        correct_origo = origo.X_to_Y > origo.Y_to_X

    inferred = lambda x: "✓" if x else "✗"
    print("ground truth : %s → %s" % (truth[0], truth[1]))
    print("ergo         : %s" % inferred(correct_ergo))
    print("ergo+        : %s" % inferred(correct_origo))
    print()


if __name__ == "__main__":
    test_tubingen_pairs()
    # run_on_pair(1)
