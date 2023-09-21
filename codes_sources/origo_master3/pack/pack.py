#!/usr/bin/env python
import sys
import os
from math import log
if __name__ == "__main__":
    from mdl import binomial_mdl
    from out import out_dag, out_scheme, out_sets, out_model
else:
    from .mdl import binomial_mdl
    from .out import out_dag, out_scheme, out_sets, out_model

import profile


bmdl = None


def binmdl(n):
    global bmdl
    if bmdl[n] == -1:
        bmdl[n] = binomial_mdl(n)
    return bmdl[n]


def isin(v, t):
    s = 0
    e = len(t) - 1
    while (s <= e):
        m = (s + e) // 2
        if t[m] == v:
            return True
        elif t[m] < v:
            s = m + 1
        elif t[m] > v:
            e = m - 1
    return False


def split(v, data):
    p = []
    n = []
    for t in data:
        if isin(v, t):
            p.append(t)
        else:
            n.append(t)
    return [p, n]


class Tree:

    def __init__(s):
        s.leaf = True
        s.neg = None
        s.pos = None
        s.attr = None
        s.score = 0
        s.saved = 0
        s.total = 0
        s.ncount = 0
        s.pcount = 0

    def split(s, neg, pos, attr):
        s.leaf = False
        s.neg = neg
        s.pos = pos
        s.attr = attr

    def update_total(s):
        if s.leaf == False:
            s.total = s.score + s.neg.update_total() + s.pos.update_total()
        else:
            s.total = s.score
        return s.total


class Scheme:

    def __init__(s, data, target, cands, K):
        """s = self
        target = attribute for tree is being generated
        cands = candidates
        K = maximum integer in the data
        """
        s.target = target
        s.tree = Tree()
        s.fragments = [Fragment([[]], split(target, data), cands, K, s.tree)]
        s.update_scores()

    def update_scores(s):
        s.scores = [s.fragments[i].best_score() + [i]
                    for i in range(len(s.fragments))]
        s.best_score = min(s.scores)

    def remove_cands(s, black):
        for f in s.fragments:
            f.remove_cands(black)
        s.update_scores()

    def split(s):
        i = s.best_score[2]
        s.fragments[i].tree.saved = -s.best_score[0]
        s.fragments = s.fragments[
            0:i] + s.fragments[i].split(s.best_score[1]) + s.fragments[(i + 1):]
        s.update_scores()

    def out(s):
        print(s.target, ':')
        for f in s.fragments:
            f.out()


class Fragment:

    def __init__(s, fr, data, cands, K, tree):
        """
        fr: empty list of list
        data: [positive rows, negative rows]
        tree: each fragment gets an instance of Tree
        """
        s.data = data
        s.fr = fr
        s.K = K
        s.cands = cands[:]
        s.tree = tree

        s.cm = [-1] * s.K
        for i in range(len(s.cands)):
            s.cm[s.cands[i]] = i

        s.tree.score = s.bits()
        s.tree.pcount = len(data[0])  # rows containing this fragment (attr)
        s.tree.ncount = len(data[1])  # rows that do not contain this fragment
        # score of splitting for all candidates
        s.score = s.calculate_score(s.calculate_freqs())

    def calculate_freqs(s):
        """computes the frequencies of items in positive and negative rows for
        this fragment
        """
        freqs = [[0] * len(s.cands), [0] * len(s.cands)]

        for t in s.data[0]:
            for i in t:
                if s.cm[i] >= 0:
                    freqs[0][s.cm[i]] += 1

        for t in s.data[1]:
            for i in t:
                if s.cm[i] >= 0:
                    freqs[1][s.cm[i]] += 1
        return freqs

    def remove_cands(s, black):
        for k in black:
            if s.cm[k] >= 0:
                s.cands[s.cm[k]] = -1
                s.score[s.cm[k]] = 1
        return

    def bits(s):
        entr = lambda x: 0 if x == 0 else -x * log(x, 2)
        pen = lambda x: binmdl(int(x)) + 1.0
        f = [len(s.data[0]), len(s.data[1])]
        n = sum(f)
        # A dirty trick to care of zero values.
        f = [x + int(x == 0) for x in f]
        return sum(map(entr, f)) - entr(n) + pen(n)

    def split_bits(s, f):
        # pos: rows where this fragment is present
        # f = [items_pres_freq_pos, items_abs_freq_pos, items_pres_freq_neg, items_abs_freq_neg]
        n1 = [f[0] + f[2], f[1] + f[3]]  # items_pres_freq, items_abs_freq
        n2 = [f[0] + f[1], f[2] + f[3]]  # items_freq_pos, items_freq_neg
        n = [n1[0] + n1[1]]
        # print f
        # print n1
        # print n2
        # print n

        if (min(n1) == 0 or min(n2) == 0):
            return 1

        # A dirty trick to care of zero values.
        f = [x + int(x == 0) for x in f]

        pen2 = 1.0 + log(s.K, 2)  # for intermediate node (data cost)
        # all the penalties are part of leaves node cost (data cost)
        # print "model cost = %.2f" % (entropy(f) - entropy(n1) - entropy(n2) + entropy(n))
        # print
        return entropy(f) - entropy(n1) - entropy(n2) + entropy(n) + penalty(n1) - penalty(n) + pen2

    def calculate_score(s, freqs):
        score = [0] * len(s.cands)

        for i in range(len(s.cands)):
            f = [0] * 4
            # frquencies of items present in positive data
            f[0] = float(freqs[0][i])
            # frquencies of items absent in positive data
            f[1] = float(len(s.data[0]) - freqs[0][i])
            # frquencies of items in negative data
            f[2] = float(freqs[1][i])
            # frquencies of items absent in negative data
            f[3] = float(len(s.data[1]) - freqs[1][i])
            score[i] = s.split_bits(f)

        return score

    def best_score(s):
        if len(s.score) == 0:
            return [1, -1]
        i = s.score.index(min(s.score))
        return [s.score[i], s.cands[i]]

    def split(s, v):
        datap = split(v, s.data[0])
        datan = split(v, s.data[1])
        cands = [x for x in s.cands if x >= 0 and x != v]

        neg = Tree()
        pos = Tree()
        s.tree.split(neg, pos, v)
        s.tree.score = 1.0 + log(s.K, 2)
        frsets = [x + [v] for x in s.fr]
        fr1 = Fragment(frsets, [datap[0], datan[0]], cands, s.K, pos)
        fr2 = Fragment(s.fr + frsets, [datap[1], datan[1]], cands, s.K, neg)
        return [fr1, fr2]

    def out(s):
        print(s.fr)
        # sys.stdout.write('[')
        # for f in s.fr:
        # if f < 0:
        # sys.stdout.write('-')
        # else:
        # sys.stdout.write('+')
        # sys.stdout.write(str(abs(f)-1))
        #sys.stdout.write(', ')
        # print ']'


def entropy(x):
    res = 0
    for y in x:
        res -= y * log(y, 2)
    return res


def penalty(x):
    res = 0
    for y in x:
        res += binmdl(int(y)) + 1.0
    return res


class Dag:

    def __init__(s, K):
        s.edges = [set() for i in range(K)]
        s.children = [set() for i in range(K)]

    def add(s, par, child):
        s.edges[child] |= set([par])
        s.update(par, s.children[child] | set([child]))

    def update(s, node, children):
        if not children <= s.children[node]:
            s.children[node] |= children
            for c in s.edges[node]:
                s.update(c, children)


def build_trees(data, K, target=[]):
    # When computing K(X) + K(Y|X), we do not want any attribute in X to
    # depend on Y
    if target:
        schemes = [None] * K
        cause = sorted(list(set(range(K)).difference(target)))
        # print
        # print "cause=", cause, " target=", target

        for citem in cause:
            temp = list(cause)
            temp.remove(citem)
            scheme = Scheme(data, citem, temp, K)
            # print citem, " __ccands__", temp
            schemes[citem] = scheme

        for titem in target:
            temp = list(target)
            temp.remove(titem)
            temp = temp + cause
            scheme = Scheme(data, titem, temp, K)
            # print titem, " __tcands__", temp
            schemes[titem] = scheme
    else:
        schemes = [Scheme(data, i, list(range(i)) + list(range(i + 1, K)), K)
                   for i in range(K)]
    deps = Dag(K)
    orig_score = sum([s.tree.score for s in schemes])
    while True:
        # best_score = [score, split_attr, fragment, leaf]
        scores = [schemes[i].best_score + [i] for i in range(K)]
        bs = min(scores)
        # print '%.3f at %d with %d' % (bs[0], bs[3], bs[1])
        if bs[0] > 0:
            break
        # split the scheme corresponding to the fragment
        schemes[bs[3]].split()
        # print 'updating graph'
        # todo(kailash): here is the place to stop the edges if any
        # bs[1] is parent and bs[3] the child
        # print "splitting parent=%d child=%d" % (bs[1], bs[3])
        deps.add(bs[1], bs[3])
        # print 'removing cands'
        for i in range(len(schemes)):
            schemes[i].remove_cands(deps.children[i])
    # print deps.edges
    return (schemes, orig_score, deps)


def build_trees_together(data, K, left_columns, right_columns, conditional):
    # When computing K(X) + K(Y|X), we do not want any attribute in X to
    # depend on Y
    schemes = [None] * K
    for col in left_columns:
        temp = list(left_columns)
        temp.remove(col)
        # print col, '->', temp
        scheme = Scheme(data, col, temp, K)
        schemes[col] = scheme

    for col in right_columns:
        temp = []
        if conditional:
            temp += left_columns
        temp += right_columns
        temp.remove(col)
        # print col, '->', temp
        scheme = Scheme(data, col, temp, K)
        schemes[col] = scheme

    deps = Dag(K)
    orig_score = sum([s.tree.score for s in schemes])
    num_causal_edges = 0
    sizes = []
    while True:
        scores = [schemes[i].best_score + [i] for i in range(K)]
        bs = min(scores)
        # print '%.3f at %d with %d' % (bs[0], bs[3], bs[1])
        if bs[0] > 0:
            break
        schemes[bs[3]].split()
        # print 'updating graph'
        if bs[3] in right_columns and bs[1] in left_columns:
            num_causal_edges += 1
            for scheme in schemes:
                scheme.tree.update_total()
            size_total = sum([s.tree.total for s in schemes])
            # print "%d x %.2f" % (num_causal_edges, size_total)
            sizes.append(size_total)
        deps.add(bs[1], bs[3])
        # print 'removing cands'
        for i in range(len(schemes)):
            schemes[i].remove_cands(deps.children[i])
    # print sizes
    return (schemes, orig_score, deps)


def compress_together(rows, left_columns, right_columns, conditional=False):
    global bmdl
    K = 1 + max(left_columns + right_columns)
    bmdl = [-1] * (1 + len(rows))
    (schemes, orig_score, deps) = build_trees_together(
        rows, K, left_columns, right_columns, conditional)
    # print column_names, ' -> ', deps.edges
    for scheme in schemes:
        scheme.tree.update_total()
    size_total = sum([s.tree.total for s in schemes])
    return size_total or 1.0


def compress(rows, column_names, target=[]):
    """ rows = list(list(int))
        column_names = list(int)
    """
    global bmdl
    K = 1 + max(column_names)
    bmdl = [-1] * (1 + len(rows))
    (schemes, orig_score, deps) = build_trees(rows, K, target)
    # print column_names, ' -> ', deps.edges
    for scheme in schemes:
        scheme.tree.update_total()
    size_total = sum([s.tree.total for s in schemes])
    return size_total or 1.0


def discover_patterns(rows, column_names):
    """ rows = list(list(int))
        column_names = list(int)
    """
    global bmdl
    K = 1 + max(column_names)
    bmdl = [-1] * (1 + len(rows))
    (schemes, orig_score, deps) = build_trees(rows, K)
    for scheme in schemes:
        scheme.tree.update_total()
    patterns = []
    for attr, edge in enumerate(deps.edges):
        if edge:
            patterns.append((edge, attr))
    return patterns


if __name__ == '__main__':
    print(compress_together([
        [1, 2],
        [1, 2],
        [1, 2],
        [1, 2],
        [1, 2],
        [1, 2],
        [1, 2],
        [1, 2],
        [1, 2],
        [1, 2],
        [1, 2],
        [1, 2],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [3, 4],
        [3, 4],
        [3, 4],
        [3, 4],
        [3, 4],
        [3, 4],
        [3, 4],
        [3, 4],
        [3, 4],
        [3, 4],
        [3, 4],
        [3, 4],
        [3],
        [3],
        [3],
        [3],
        [3],
        [3],
        [3],
        [3],
        [3],
        [3]
        ], [0, 1, 2], [3, 4], True))
