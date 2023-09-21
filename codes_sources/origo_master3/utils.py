#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys


def concatenate(X, Y):
    assert len(X[0]) == len(Y[0])
    num_rows = len(X[0])
    rows = []
    for i in range(num_rows):
        row = []
        for j in range(len(X)):
            row.append(X[j][i])

        for j in range(len(Y)):
            row.append(Y[j][i])
        rows.append(row)
    return rows


def normalise(data):
    max_val = max(data)
    norm = lambda x: x * 1.0 / max_val
    return list(map(norm, data))


def pad_sequence(sequence, size, padding):
    extension = [padding] * abs(len(sequence) - size)
    padded_sequence = sequence + extension
    return padded_sequence


def partition(collection):
    if len(collection) == 1:
        yield [collection]
        return

    first = collection[0]
    for smaller in partition(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[first] + subset] + smaller[n + 1:]
        # put `first` in its own subset
        yield [[first]] + smaller


def progress(count, total, suffix=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))
    sys.stdout.flush()  # As suggested by Rom Ruben


def reverse_argsort(X):
    indices = list(range(len(X)))
    indices.sort(key=X.__getitem__, reverse=True)
    return indices
