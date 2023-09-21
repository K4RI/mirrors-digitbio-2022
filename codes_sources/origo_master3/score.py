#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Computes different scores.
"""
if __name__ == "__main__":
    from dc import dc
    from pack import pack
    from store import Store
    from bagging import bagging

else:
    from .dc import dc
    from .pack import pack
    from .store import Store
    from .bagging import bagging



class Score(object):

    def __init__(self, X_to_Y, Y_to_X):
        self.X_to_Y = X_to_Y
        self.Y_to_X = Y_to_X

    def __str__(self):
        return "X->Y = %.8f\tY->X = %.8f" % (self.X_to_Y, self.Y_to_X)


def score_ergo_origo_dc(X, Y, store):
    X_rows = store.project(X)
    Y_rows = store.project(Y)
    X_Y_rows = Store.combine(X_rows, Y_rows)

    pack_X = pack.compress(X_rows, X)
    pack_Y = pack.compress(Y_rows, Y)

    pack_X_and_Y = pack.compress_together(X_Y_rows, X, Y, False)
    pack_X_to_Y = pack.compress(X_Y_rows, X + Y, Y)
    pack_Y_to_X = pack.compress(X_Y_rows, Y + X, X)

    pack_Y_given_X = pack_X_to_Y - pack_X
    pack_X_given_Y = pack_Y_to_X - pack_Y

    # to make it fair for encoding cost when computing K(Y|X) and K(Y)
    ergo_X_to_Y = pack_Y_given_X / max(pack_X_and_Y - pack_X, 1.0)
    ergo_Y_to_X = pack_X_given_Y / max(pack_X_and_Y - pack_Y, 1.0)

    dc_result = dc(X_rows, Y_rows)
    dc_score = Score(dc_result[0], dc_result[1])
    ergo_score = Score(ergo_X_to_Y, ergo_Y_to_X)
    origo_score = Score(pack_X_to_Y, pack_Y_to_X)

    return ergo_score, origo_score, dc_score


def score_ergo_origo_dc_old(X, Y, store):
    X_rows = store.project(X)
    Y_rows = store.project(Y)
    X_Y_rows = Store.combine(X_rows, Y_rows)

    pack_X = pack.compress(X_rows, X)
    pack_Y = pack.compress(Y_rows, Y)

    pack_X_and_Y = pack.compress_together(X_Y_rows, X, Y, False)
    pack_X_to_Y = pack.compress(X_Y_rows, X + Y, Y)
    pack_Y_to_X = pack.compress(X_Y_rows, Y + X, X)

    pack_Y_given_X = pack_X_to_Y - pack_X
    pack_X_given_Y = pack_Y_to_X - pack_Y

    # to make it fair for encoding cost when computing K(Y|X) and K(Y)
    ergo_X_to_Y = pack_Y_given_X / max(pack_X_and_Y - pack_X, 1.0)
    ergo_Y_to_X = pack_X_given_Y / max(pack_X_and_Y - pack_Y, 1.0)

    origo_X_to_Y = pack_X_to_Y / pack_X_and_Y
    origo_Y_to_X = pack_Y_to_X / pack_X_and_Y

    dc_result = dc(X_rows, Y_rows)
    dc_score = Score(dc_result[0], dc_result[1])
    ergo_score = Score(ergo_X_to_Y, ergo_Y_to_X)
    origo_score = Score(origo_X_to_Y, origo_Y_to_X)

    return ergo_score, origo_score, dc_score


def score_origo_bagging(X, Y, store):
    X_rows = store.project(X)
    Y_rows = store.project(Y)
    X_Y_rows = Store.combine(X_rows, Y_rows)

    pack_X_to_Y = pack.compress(X_Y_rows, X + Y, Y)
    pack_Y_to_X = pack.compress(X_Y_rows, Y + X, X)

    baggs_XY = bagging(X_Y_rows, X + Y, Y)
    baggs_YX = bagging(X_Y_rows, X + Y, X)

    origo = Score(pack_X_to_Y, pack_Y_to_X)
    origo_bagging = Score(baggs_XY, baggs_YX)

    return origo, origo_bagging
