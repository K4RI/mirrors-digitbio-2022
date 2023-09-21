#!/usr/bin/env python
# -*- coding: utf-8 -*-
import random

if __name__ == "__main__":
    from pack import pack

else:
    from .pack import pack


def bootstrap(rows):
    return [random.choice(rows) for i in range(len(rows))]


def compress(rows, lcols, rcols, at_once):
    if at_once:
        size = pack.compress_together(rows, lcols, rcols, False)
    else:
        size = pack.compress(rows, lcols, rcols)
    size = max(size, 1.0)
    return size


def bagging(rows, lcols, rcols, at_once=False):
    B = 50
    return [compress(bootstrap(rows), lcols, rcols, at_once)
            for i in range(B)]


if __name__ == "__main__":
    rows = [
        [0],
        [0, 1],
        [1, 2, 3],
        [1, 2]
    ]

    print(bootstrap(rows))
