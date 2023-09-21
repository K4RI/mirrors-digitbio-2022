#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Mines correlated pairs from fimi dataset. It expects the <dataset> directory inside the `data` directory. In the <dataset> directory, there should be a <dataset.dat> data file and <dataset.labels> labels file.
"""
import os
from operator import itemgetter
from subprocess import check_output

from score import score_ergo_origo_dc
from store import Store
from utils import partition


def run_opus_miner(fpath):
    opus_miner_dir = os.path.join(os.path.dirname(__file__), "opus")
    opus_miner_exec = os.path.join(opus_miner_dir, "opus_miner")
    output_path = os.path.join(opus_miner_dir, "out.txt")
    cmd = [opus_miner_exec, fpath, output_path]
    print(" ".join(c for c in cmd))
    check_output(cmd)

    with open(output_path, "r") as reader:
        rows = reader.readlines()

    is_next_interesting = False
    self_sufficient_itemsets = []
    for row in rows:
        row = row.rstrip()
        if row.startswith("SELF-SUFFICIENT ITEMSETS:"):
            is_next_interesting = True
            continue

        if is_next_interesting:
            if not row:
                break
            itemset = row.split('[')[0].rstrip()
            itemset = list(map(int, itemset.split(',')))
            self_sufficient_itemsets.append(itemset)
    return self_sufficient_itemsets


def mine_using_opus_miner(dataset):
    justify = 100
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    dat_path = os.path.join(data_dir, dataset, "%s.dat" % dataset)
    labels_path = os.path.join(data_dir, dataset, "%s.labels" % dataset)

    with open(labels_path, "r") as reader:
        labels = reader.readlines()
        labels = [label.rstrip() for label in labels]

    print("running opus miner", end=' ')
    self_sufficient_itemsets = run_opus_miner(dat_path)
    print(" [DONE]")
    store = Store(dat_path)

    print("searching for causal pairs")
    causal_pairs = []
    for itemset in self_sufficient_itemsets:
        # create all the partitions of itemset of cardinality two
        for pair in partition(itemset):
            if len(pair) != 2:
                continue

            colnames_x, colnames_y = pair
            new_colnames_x = list(range(len(colnames_x)))
            new_colnames_y = list(range(len(colnames_x), len(
                colnames_x) + len(colnames_y)))

            new_colnames = dict(
                list(zip(colnames_x + colnames_y, new_colnames_x + new_colnames_y)))

            rows = store.project(colnames_x + colnames_y)
            for row in rows:
                for i, old_colname in enumerate(row):
                    if old_colname != None:
                        row[i] = new_colnames[old_colname]

            pair_store = Store(None, rows)
            _, origo, _ = score_ergo_origo_dc(
                new_colnames_x, new_colnames_y, pair_store)
            delta = abs(origo.X_to_Y - origo.Y_to_X)

            words_x = " ".join(labels[i] for i in colnames_x)
            words_y = " ".join(labels[i] for i in colnames_y)
            if origo.X_to_Y < origo.Y_to_X:
                if dataset == "adult" and len(words_y) > 5:
                    continue
                to_print = "%s → %s (%f)" % (words_x, words_y, delta)
                causal_pairs.append((to_print, delta))
                print(to_print)
            elif origo.X_to_Y > origo.Y_to_X:
                if dataset == "adult" and len(words_x) > 5:  # for income only
                    continue
                to_print = "%s → %s (%f)" % (words_y, words_x, delta)
                causal_pairs.append((to_print, delta))
                print(to_print)

    sorted_causal_pairs = sorted(causal_pairs, key=itemgetter(1), reverse=True)
    out_path = os.path.join(os.path.dirname(__file__),
                            "results", "%s_out.txt" % dataset)
    with open(out_path, "w") as writer:
        for i, causal_pair in enumerate(causal_pairs):
            causal_pair = causal_pair[0].rjust(justify)
            causal_pair_sorted = sorted_causal_pairs[i][0].rjust(justify)
            row = "%s %s" % (causal_pair, causal_pair_sorted)
            writer.write(row + "\n")


if __name__ == "__main__":
    mine_using_opus_miner("icdm")
    mine_using_opus_miner("adult")
