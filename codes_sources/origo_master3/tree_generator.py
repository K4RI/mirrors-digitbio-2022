#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Generates tree patterns, which are planted in data. Left side of the decision tree represents attribute=1 decision and right side attribute=0. We start with uniform distribution of leafs.
"""
import collections
import math
import numpy
import random
import time

from sequence import Queue, Stack


MAX_CONFIG_ATTEMPTS = 10000


def compositions(t, s):
    """credit: https://gist.github.com/jasonmc/989158/682cfe1c25d39d5526acaeacc1426d287ef4f5de
    """
    q = [0] * (t + 1)
    r = None
    q[0] = s
    while True:
        if t == 0 and s == 1:
            q[0] = 1
            yield q
            break

        yield list(q)
        if s == 0:
            break
        if q[0] == 0:
            if r == t:
                break
            else:
                q[0] = q[r] - 1
                q[r] = 0
                r = r + 1
        else:
            q[0] = q[0] - 1
            r = 1
        q[r] = q[r] + 1


def get_random():
    random.seed(time.time())
    return random.random()


def pick_random(seq):
    random.seed(time.time())
    return random.choice(seq)


def gen_rand_prob_pair():
    p = round(pick_random(numpy.arange(0.01, 0.49, 0.01)), 2)
    q = 1 - p
    return [p, q]


def gen_rand_prob():
    return pick_random(gen_rand_prob_pair())


def log2(x):
    return math.log(x, 2)


def safely_remove(container, item):
    # container is a list type
    try:
        container.remove(item)
    except ValueError:
        pass


def remove_container(src_container, garbage):
    return [x for x in src_container if x not in garbage]


def classify(tree, sample):
    # sample is a dictionary, attribute -> value
    # or an implicit array, where attribute is its index

    # handle empty tree
    if tree.attribute == None:
        label = pick_random([0, 1])
        return label

    label = pick_random([0, 1])
    current_node = tree
    while True:
        current_attr = current_node.attribute
        value = sample[current_attr]
        next_node = current_node.left if value == 1 else current_node.right

        if next_node.is_leaf():
            prob_leaf_presence = current_node.left_prob_leaf_presence if value == 1 else current_node.right_prob_leaf_presence
            label = 1 if get_random() <= prob_leaf_presence else 0
            break
        else:
            current_node = next_node
    return label


class Node(object):

    def __init__(self, attribute, parent, level):
        self.left = None
        self.right = None
        self.level = level      # for pretty print
        self.parent = parent
        self.attribute = attribute
        self.left_prob_instances = 0.0
        self.right_prob_instances = 0.0
        self.left_prob_leaf_presence = 0.0
        self.right_prob_leaf_presence = 0.0

    def is_leaf(self):
        return self.attribute == None

    def make_leaf(self):
        self.__init__(None, self.parent, self.level)

    def has_one_of(self, attrs):
        found = False
        sequence = Stack()
        sequence.add(self.left)
        sequence.add(self.right)
        while not sequence.empty():
            child = sequence.get()
            if child.attribute in attrs:
                found = True
                break

            if not child.is_leaf():
                sequence.add(child.left)
                sequence.add(child.right)
        return found


class TreePatternGenerator(object):

    def __init__(self, causes, effects, prob_causal_dependency, prob_split, max_height=0, depth_first=False, num_causal_edges=None):
        """We choose one attribute either from causes or effects for every split decision. The probability of causal dependency indicates how often do we choose causes as the parents of effects. Obviously the higher the causal dependency, the likelier we are to infer the correct causal relationship. We use the given ordering to find the parents for a node. We expect fairly high (but not too high) probability of dependency and split to match the the maximum number of causal edges.
        """
        self.causes = causes
        self.effects = effects
        self.max_height = max_height
        self.prob_split = prob_split
        self.depth_first = depth_first
        self.num_causal_edges = num_causal_edges
        self.causal_edges_per_effect = None
        self.causal_edges_per_effect_added = None
        self.prob_causal_dependency = prob_causal_dependency

    @staticmethod
    def calculate_information_gain(parent, child):
        parent_entropy = - parent.left_prob_instances * \
            log2(parent.left_prob_instances) - \
            parent.right_prob_instances * log2(parent.right_prob_instances)
        child_entropy_left = - child.left_prob_leaf_presence * log2(child.left_prob_leaf_presence) - (
            1 - child.left_prob_leaf_presence) * log2(1 - child.left_prob_leaf_presence)
        child_entropy_right = - child.right_prob_leaf_presence * log2(child.right_prob_leaf_presence) - (
            1 - child.right_prob_leaf_presence) * log2(1 - child.right_prob_leaf_presence)
        information_gain = parent_entropy - child.left_prob_instances * \
            child_entropy_left - child.right_prob_instances * child_entropy_right
        return information_gain

    def _count_causal_edges(self, tree):
        # number of causal variables in the tree
        used_causes = set()
        queue = Queue()
        queue.add(tree)
        while not queue.empty():
            node = queue.get()
            attr = node.attribute

            if attr in self.causes:
                used_causes.add(attr)

            if node.left and not node.left.is_leaf():
                queue.add(node.left)
            if node.right and not node.right.is_leaf():
                queue.add(node.right)
        num_causal_edges = len(used_causes)
        return num_causal_edges

    def count_causal_edges(self, trees):
        return sum(self._count_causal_edges(tree) for tree in trees[len(self.causes):])

    def configure_child(self, child):
        prob_pair = gen_rand_prob_pair()
        child.left_prob_leaf_presence = gen_rand_prob()
        child.right_prob_leaf_presence = gen_rand_prob()
        child.left_prob_instances = pick_random(prob_pair)
        child.right_prob_instances = 1 - child.left_prob_instances

    def generate_edge_constrained_tree(self, leaf, cand_causes, cand_effects):
        # until we get the fixed number of causal edges, we do not look into
        # other parameters
        leaf_idx = len(cand_effects)
        root_attr = None
        if self.causal_edges_per_effect[leaf_idx] != 0:
            root_attr = pick_random(cand_causes)
            self.causal_edges_per_effect_added[leaf_idx] += 1
            # print "\n\nadding ", root_attr
        else:
            if get_random() <= self.prob_split and get_random() <= 1 - self.prob_causal_dependency and cand_effects:
                root_attr = pick_random(cand_effects)

        root = Node(root_attr, None, 0)

        if root_attr != None:
            prob_pair = gen_rand_prob_pair()
            root.left = Node(None, root, 1)
            root.right = Node(None, root, 1)
            root.left_prob_leaf_presence = gen_rand_prob()
            root.right_prob_leaf_presence = gen_rand_prob()
            root.left_prob_instances = pick_random(prob_pair)
            root.right_prob_instances = 1 - root.left_prob_instances
            self.split_edge_constrained(root, cand_causes, cand_effects)

        # print "LEAF =", str(leaf), " TREE_HEIGHT =", str(self.get_height(root))
        # TreePatternGenerator.pretty_print(root)
        # print "\n"

        return root

    def generate_edge_constrained_trees(self):
        max_num_causal_edges = len(self.causes) * len(self.effects)
        assert 0 <= self.num_causal_edges <= max_num_causal_edges

        edge_dist = list(compositions(
            len(self.effects) - 1, self.num_causal_edges))
        valid_edge_dist = []
        for dist in edge_dist:
            valid = True
            for num_causal_edges in dist:
                if num_causal_edges > len(self.causes):
                    valid = False
                    break
            if valid:
                valid_edge_dist.append(dist)

        random.seed(time.time())
        rand_idx = random.randint(0, len(valid_edge_dist) - 1)

        self.causal_edges_per_effect = valid_edge_dist[rand_idx]
        self.causal_edges_per_effect_added = [0] * len(self.effects)

        trees_for_causes = [self.generate_normal_tree(
            cause, self.causes[:pos], [], True) for pos, cause in enumerate(self.causes)]
        trees_for_effects = [self.generate_edge_constrained_tree(effect, list(self.causes), self.effects[
            :pos]) for pos, effect in enumerate(self.effects)]
        return trees_for_causes + trees_for_effects

    def generate_normal_tree(self, leaf, cand_causes, cand_effects, leaf_is_cause):
        root_attr = None
        if get_random() <= self.prob_split:
            root_attr = self.get_split_attr(
                cand_causes, cand_effects, [], leaf_is_cause)
        root = Node(root_attr, None, 0)

        if root_attr != None:
            prob_pair = gen_rand_prob_pair()
            root.left = Node(None, root, 1)
            root.right = Node(None, root, 1)
            root.left_prob_leaf_presence = gen_rand_prob()
            root.right_prob_leaf_presence = gen_rand_prob()
            root.left_prob_instances = pick_random(prob_pair)
            root.right_prob_instances = 1 - root.left_prob_instances

            if self.max_height > 1:
                self.split_normal(root, cand_causes,
                                  cand_effects, leaf_is_cause)
        # print "LEAF =", str(leaf), " TREE_HEIGHT =", str(self.get_height(root))
        # TreePatternGenerator.pretty_print(root)
        # print "\n"

        return root

    def generate_normal_trees(self):
        trees_for_causes = [self.generate_normal_tree(
            cause, self.causes[:pos], [], True) for pos, cause in enumerate(self.causes)]
        trees_for_effects = [self.generate_normal_tree(effect, list(self.causes), self.effects[
            :pos], False) for pos, effect in enumerate(self.effects)]
        return trees_for_causes + trees_for_effects

    def generate_trees(self):
        if self.num_causal_edges != None:
            return self.generate_edge_constrained_trees()
        return self.generate_normal_trees()

    def get_height(self, tree):
        # height = number of edges from the root to the leaf
        height = 0
        queue = Queue()
        queue.add(tree)
        while not queue.empty():
            node = queue.get()
            if node.level > height:
                height = node.level
            if node.left:
                queue.add(node.left)
            if node.right:
                queue.add(node.right)
        return height

    def get_split_attr(self, cand_causes, cand_effects, parents, leaf_is_cause=False):
        attr = None
        try:
            if leaf_is_cause:
                if get_random() <= 1 - self.prob_causal_dependency:
                    pool = list(set(cand_causes).difference(parents))
                    attr = pick_random(pool)
            else:
                # magic johnson here
                if get_random() <= self.prob_causal_dependency:
                    pool = list(set(cand_causes).difference(parents))
                    attr = pick_random(pool)
                else:
                    pool = list(set(cand_effects).difference(parents))
                    attr = pick_random(pool)
        except IndexError:
            pass
        return attr

    @staticmethod
    def _pretty_print(tree, indent):
        if tree.is_leaf():
            print(indent + "x")
            return

        if tree.left.is_leaf() and tree.right.is_leaf():
            print(indent + str(tree.attribute) + "(%d%%,%d%%)" % (int(tree.left_prob_leaf_presence * 100), int(tree.right_prob_leaf_presence * 100)) + "[%d%%,%d%%]" % (int(tree.left_prob_instances * 100), int(tree.right_prob_instances * 100)))
        elif tree.left.is_leaf():
            print(indent + str(tree.attribute) + "(%d%%,)" % int(tree.left_prob_leaf_presence * 100) + "[%d%%,%d%%]" % (int(tree.left_prob_instances * 100), int(tree.right_prob_instances * 100)))
        elif tree.right.is_leaf():
            print(indent + str(tree.attribute) + "(,%d%%)" % int(tree.right_prob_leaf_presence * 100) + "[%d%%,%d%%]" % (int(tree.left_prob_instances * 100), int(tree.right_prob_instances * 100)))
        else:
            print(indent + str(tree.attribute) + "[%d%%,%d%%]" % (int(tree.left_prob_instances * 100), int(tree.right_prob_instances * 100)))
        TreePatternGenerator._pretty_print(
            tree.left, indent + "|___________________")
        TreePatternGenerator._pretty_print(
            tree.right, indent + "|___________________")

    @staticmethod
    def pretty_print(tree):
        TreePatternGenerator._pretty_print(tree, "")

    def split_edge_constrained(self, root, cand_causes, cand_effects):
        leaf_idx = len(cand_effects)
        sequence = Stack() if self.depth_first else Queue()
        sequence.add((root.left, [root.attribute]))
        sequence.add((root.right, [root.attribute]))

        used_causes = set([root.attribute])
        while not sequence.empty():
            child, parents = sequence.get()
            child_attr = None
            if self.causal_edges_per_effect_added[leaf_idx] < self.causal_edges_per_effect[leaf_idx]:
                pool = list(set(cand_causes).difference(used_causes))
                if pool:
                    child_attr = pick_random(pool)
                    used_causes.add(child_attr)
                    # print "adding ", child_attr
                    self.causal_edges_per_effect_added[leaf_idx] += 1
            else:
                # normal routine here
                if get_random() > self.prob_split or len(parents) >= self.max_height:
                    continue

                if get_random() <= (1 - self.prob_causal_dependency):
                    pool = list(set(cand_effects).difference(parents))
                    if pool:
                        child_attr = pick_random(pool)

            if child_attr == None:
                continue

            num_config_attempts = 0
            while True:
                self.configure_child(child)
                information_gain = TreePatternGenerator.calculate_information_gain(
                    child.parent, child)

                if information_gain > 0.0:
                    child.attribute = child_attr
                    child.left = Node(None, child, child.level + 1)
                    child.right = Node(None, child, child.level + 1)

                    sequence.add((child.left, parents + [child_attr]))
                    sequence.add((child.right, parents + [child_attr]))
                    break
                else:
                    num_config_attempts += 1

                if num_config_attempts == MAX_CONFIG_ATTEMPTS:
                    child.make_leaf()
                    break

    def split_normal(self, root, cand_causes, cand_effects, leaf_is_cause):
        """We do a breadth-first traversal of the tree. We decide whether to split at each step and randomly find a configuration that yields information gain. The configuration comprises of the following.
        1. number of instances on left and right branches
        2. prob(leaf=1) on left and right branches
        """
        sequence = Stack() if self.depth_first else Queue()
        sequence.add((root.left, [root.attribute]))
        sequence.add((root.right, [root.attribute]))

        while not sequence.empty():
            child, parents = sequence.get()

            # at this level, we reach max height
            if len(parents) >= self.max_height:
                continue

            # should we split here
            if get_random() > self.prob_split:
                continue

            # do something with the parents
            child_attr = self.get_split_attr(
                cand_causes, cand_effects, parents, leaf_is_cause)

            if child_attr == None:
                continue

            num_config_attempts = 0
            while True:
                self.configure_child(child)
                information_gain = TreePatternGenerator.calculate_information_gain(
                    child.parent, child)

                if information_gain > 0.0:
                    child.attribute = child_attr
                    child.left = Node(None, child, child.level + 1)
                    child.right = Node(None, child, child.level + 1)

                    sequence.add((child.left, parents + [child_attr]))
                    sequence.add((child.right, parents + [child_attr]))
                    break
                else:
                    num_config_attempts += 1

                if num_config_attempts == MAX_CONFIG_ATTEMPTS:
                    child.make_leaf()
                    break


def test_tree():
    """generate the following tree and test the features.
            1
          /   \
         2     3
        / \   / \
    """
    print("generating a synthetic tree")
    root = Node(1, None, 0)
    root.left_prob_leaf_presence = 0.4
    root.right_prob_leaf_presence = 0.7
    root.left_prob_instances = 0.6
    root.right_prob_instances = 0.4

    left = Node(2, root, 1)
    left.left = Node(None, left, 1)
    left.right = Node(None, left, 1)
    left.left_prob_leaf_presence = 0.8
    left.right_prob_leaf_presence = 0.2
    left.left_prob_instances = 0.7
    left.right_prob_instances = 0.3

    right = Node(3, root, 1)
    right.left = Node(None, right, 1)
    right.right = Node(None, right, 1)
    right.left_prob_leaf_presence = 0.8
    right.right_prob_leaf_presence = 0.9
    right.left_prob_instances = 0.3
    right.right_prob_instances = 0.7

    root.left = left
    root.right = right

    TreePatternGenerator.pretty_print(root)

    print("labelling the data")
    # 80 % probability of getting label 1
    num_ones = 0
    for i in range(100000):
        val = classify(root, {1: 0, 2: 1, 3: 0})
        num_ones += val

    assert num_ones >= 0.75  # just to be safe with randomisation
    print("alright, alright, alright")


if __name__ == '__main__':
    num_causes = 5
    num_effects = 5
    causes = list(range(num_causes))
    effects = list(range(num_causes, num_causes + num_effects))
    print(causes, effects)
    pattern_generator = TreePatternGenerator(
        causes, effects, 0.5, 0.6, 9, num_causal_edges=3)
    print(pattern_generator.count_causal_edges(pattern_generator.generate_trees()))
    print()
