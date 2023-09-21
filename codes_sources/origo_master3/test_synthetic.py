#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""We generate synthetic data on a per transaction basis. For each transaction, the synthetic tree pattern dictates the presence or absence of an item. All the experiments reported on synthetic data are in this module.
"""

from collections import Counter
from copy import deepcopy
import random

import numpy as np

from score import score_ergo_origo_dc, score_origo_bagging
from store import Store
from tree_generator import classify, TreePatternGenerator
from utils import progress


def generate_data(nrows, causes, effects, dependency, prob_split, max_height):
    num_causes = len(causes)
    num_effects = len(effects)
    num_variables = num_causes + num_effects

    pattern_generator = TreePatternGenerator(
        causes, effects, dependency, prob_split, max_height)
    trees = pattern_generator.generate_trees()
    ncausal_edges = pattern_generator.count_causal_edges(trees)

    rows = []
    row_count = 0
    while True:
        if row_count == nrows:
            break

        row_vector = [0] * num_variables
        if num_variables == 2:
            tree_effect = trees[num_causes]
            if tree_effect.attribute != None:
                # prob_presence = tree_effect.left_prob_instances
                row_vector[0] = 1 if random.random() <= 0.9 else 0
                row_vector[1] = classify(tree_effect, row_vector)
            else:
                row_vector[0] = random.choice([0, 1])
                row_vector[1] = random.choice([0, 1])
        else:
            for i in range(num_causes):
                tree = trees[i]
                cause = causes[i]
                row_vector[cause] = classify(tree, row_vector)

            for i in range(num_effects):
                tree = trees[num_causes + i]
                effect = effects[i]
                row_vector[effect] = classify(tree, row_vector)

        row = [i for i in range(num_variables) if row_vector[i] == 1]
        rows.append(row)
        row_count += 1
    return rows, ncausal_edges


def test_scalability_asymmetric():
    num_rows = 5000
    num_simulations = 500
    pairs = [(5, 1), (5, 3), (5, 5), (5, 7), (5, 10)]
    prob_split = 1.0
    prob_dep = 0.7
    max_height = 5  # to control the tree generation time

    accuracy_pairs = []
    progress(0, len(pairs))

    for i, pair in enumerate(pairs):
        causes = list(range(pair[0]))
        effects = list(range(pair[0], sum(pair)))

        num_correct_dc = 0
        num_correct_ergo = 0
        num_correct_origo = 0

        num_pairs_type_one = int(0.5 * num_simulations)
        for _ in range(num_pairs_type_one):
            rows, _ = generate_data(num_rows, causes, effects, prob_dep,
                                    prob_split, max_height)
            store = Store(None, rows, causes + effects)
            ergo, origo, dc = score_ergo_origo_dc(causes, effects, store)

            num_correct_dc += int(dc.X_to_Y < dc.Y_to_X)
            num_correct_ergo += int(ergo.X_to_Y < ergo.Y_to_X)
            num_correct_origo += int(origo.X_to_Y < origo.Y_to_X)

        causes = list(range(pair[1]))
        effects = list(range(pair[1], sum(pair)))

        num_pairs_type_two = num_simulations - num_pairs_type_one
        for _ in range(num_pairs_type_two):
            rows, _ = generate_data(num_rows, causes, effects, prob_dep,
                                    prob_split, max_height)
            store = Store(None, rows, causes + effects)
            ergo, origo, dc = score_ergo_origo_dc(causes, effects, store)

            num_correct_dc += int(dc.X_to_Y < dc.Y_to_X)
            num_correct_ergo += int(ergo.X_to_Y < ergo.Y_to_X)
            num_correct_origo += int(origo.X_to_Y < origo.Y_to_X)

        accuracy_dc = num_correct_dc / num_simulations
        accuracy_ergo = num_correct_ergo / num_simulations
        accuracy_origo = num_correct_origo / num_simulations

        print((accuracy_ergo, accuracy_origo, accuracy_dc))
        accuracy_pairs.append((accuracy_ergo, accuracy_origo, accuracy_dc))

        progress(i + 1, len(pairs))

    print()
    print(accuracy_pairs)
    print("writing results in `results/scalability_asymmetric.dat`")
    print("FORMAT: <cardinality> <accuracy_ergo> <accuracy_origo> <accuracy_dc>")
    with open("results/scalability_asymmetric.dat", "w") as writer:
        for i, pair in enumerate(pairs):
            row = [pair[1], accuracy_pairs[i][0],
                   accuracy_pairs[i][1], accuracy_pairs[i][2]]
            row = " ".join(str(item) for item in row)
            writer.write(row + "\n")


def test_scalability_symmetric():
    num_rows = 5000
    num_simulations = 500
    pairs = [(1, 1), (3, 3), (5, 5), (7, 7), (10, 10)]
    prob_split = 1.0
    prob_dep = 0.7
    max_height = 5  # to control the tree generation time

    accuracy_pairs = []
    progress(0, len(pairs))

    for i, pair in enumerate(pairs):
        causes = list(range(pair[0]))
        effects = list(range(pair[0], sum(pair)))

        num_correct_dc = 0
        num_correct_ergo = 0
        num_correct_origo = 0
        for _ in range(num_simulations):

            rows, _ = generate_data(num_rows, causes, effects, prob_dep,
                                    prob_split, max_height)
            store = Store(None, rows, causes + effects)
            ergo, origo, dc = score_ergo_origo_dc(causes, effects, store)

            num_correct_dc += int(dc.X_to_Y < dc.Y_to_X)
            num_correct_ergo += int(ergo.X_to_Y < ergo.Y_to_X)
            num_correct_origo += int(origo.X_to_Y < origo.Y_to_X)

        accuracy_dc = num_correct_dc / num_simulations
        accuracy_ergo = num_correct_ergo / num_simulations
        accuracy_origo = num_correct_origo / num_simulations

        print((accuracy_ergo, accuracy_origo, accuracy_dc))
        accuracy_pairs.append((accuracy_ergo, accuracy_origo, accuracy_dc))

        progress(i + 1, len(pairs))

    print()
    print(accuracy_pairs)
    print("writing the results in `results/scalability_symmetric.dat`")
    print("FORMAT: <num_variables> <accuracy_ergo> <accuracy_origo> <accuracy_dc>")
    with open("results/scalability_symmetric.dat", "w") as writer:
        for i, pair in enumerate(pairs):
            row = [pair[0], accuracy_pairs[i][0],
                   accuracy_pairs[i][1], accuracy_pairs[i][2]]
            row = " ".join(str(item) for item in row)
            writer.write(row + "\n")


def test_type2_error():
    # false negatives: telling that a pregnant woman is not pregnant
    # analogy: telling a causal data is not causal
    num_causes = 3
    num_effects = 3
    num_rows = 5000
    prob_split = 1.0

    causes = list(range(num_causes))
    effects = list(range(num_causes, num_causes + num_effects))
    max_height = 5  # for speeding tree generation

    num_sim_cutoff = 250
    num_sim_power = 250
    dependencies = np.arange(0.0, 1.1, 0.1)

    cutoffs_ergo = []
    cutoffs_origo = []
    cutoffs_dc = []

    powers_ergo = []
    powers_origo = []
    powers_dc = []

    progress(0, len(dependencies))

    for prog, dependency in enumerate(dependencies):
        deltas_cutoff_ergo = []
        deltas_cutoff_origo = []
        deltas_cutoff_dc = []

        deltas_power_ergo = []
        deltas_power_origo = []
        deltas_power_dc = []

        for i in range(num_sim_cutoff):
            rows, _ = generate_data(num_rows, causes, effects,
                                    dependency, prob_split, max_height)
            # randomise the effects
            randomised_rows = []
            for row in rows:
                new_row = [item for item in row if item in causes]
                for effect in effects:
                    if random.random() <= 0.5:
                        new_row.append(effect)
                randomised_rows.append(new_row)

            store = Store(None, randomised_rows, causes + effects)
            ergo, origo, dc = score_ergo_origo_dc(causes, effects, store)

            delta_ergo = abs(ergo.Y_to_X - ergo.X_to_Y)
            delta_origo = abs(origo.Y_to_X - origo.X_to_Y)
            delta_dc = abs(dc.Y_to_X - dc.X_to_Y)

            deltas_cutoff_ergo.append(delta_ergo)
            deltas_cutoff_origo.append(delta_origo)
            deltas_cutoff_dc.append(delta_dc)

        # 95% of datasets have delta more than cutoff_score
        cutoff_ergo = np.percentile(deltas_cutoff_ergo, 5)
        cutoff_origo = np.percentile(deltas_cutoff_origo, 5)
        cutoff_dc = np.percentile(deltas_cutoff_dc, 5)

        cutoffs_ergo.append(cutoff_ergo)
        cutoffs_origo.append(cutoff_origo)
        cutoffs_dc.append(cutoff_dc)

        for i in range(num_sim_power):
            rows, _ = generate_data(num_rows, causes, effects,
                                    dependency, prob_split, max_height)
            store = Store(None, rows, causes + effects)
            ergo, origo, dc = score_ergo_origo_dc(causes, effects, store)

            delta_ergo = abs(ergo.Y_to_X - ergo.X_to_Y)
            delta_origo = abs(origo.Y_to_X - origo.X_to_Y)
            delta_dc = abs(dc.Y_to_X - dc.X_to_Y)

            deltas_power_ergo.append(delta_ergo)
            deltas_power_origo.append(delta_origo)
            deltas_power_dc.append(delta_dc)

        # % of cause-effect pairs having delta greater than cutoff
        power_ergo = sum(deltas_power_ergo > cutoff_ergo) / num_sim_power
        power_origo = sum(deltas_power_origo > cutoff_origo) / num_sim_power
        power_dc = sum(deltas_power_dc > cutoff_dc) / num_sim_power

        powers_ergo.append(power_ergo)
        powers_origo.append(power_origo)
        powers_dc.append(power_dc)

        progress(prog + 1, len(dependencies))

    print("\nwriting the results in `results/statistical_power.dat`")
    print("FORMAT: <dependency> <power_ergo> <power_origo> <power_dc>")
    with open("results/statistical_power.dat", "w") as writer:
        for i, dependency in enumerate(dependencies):
            writer.write("%f %f %f %f\n" % (dependency, powers_ergo[
                         i], powers_origo[i], powers_dc[i]))


def test_type1_error():
    # false positives: telling a man that he is pregnant
    # analogy: telling a random data it is causal
    # Pr(reject H0 when H0 is true)
    # H0 = delta_randomised_data < delta_original_data
    #    = random data cannot have greater than same delta as the causal data
    # Type 1 error = in what percentage of data do we observe H0 to be False
    num_causes = 3
    num_effects = 3
    prob_split = 1.0
    num_rows = 5000
    num_rand_data = 500
    dependency = 1.0

    causes = list(range(num_causes))
    effects = list(range(num_causes, num_causes + num_effects))
    max_height = num_causes + num_effects - 1

    progress(0, num_rand_data + 1)

    causal_data, _ = generate_data(
        num_rows, causes, effects, dependency, prob_split, max_height)
    store = Store(None, causal_data, causes + effects)
    _, origo, _ = score_ergo_origo_dc(causes, effects, store)
    delta = abs(origo.Y_to_X - origo.X_to_Y)
    deltas = [delta]

    progress(1, num_rand_data + 1)

    for i in range(num_rand_data):
        new_store = deepcopy(store)
        new_store.swap_randomise()
        _, origo, _ = score_ergo_origo_dc(causes, effects, new_store)
        delta = abs(origo.Y_to_X - origo.X_to_Y)
        deltas.append(delta)
        progress(i + 2, num_rand_data + 1)

    print()
    print(deltas)
    num_as_good_deltas = 0
    delta_causal_data = deltas[0]
    for delta_random_data in deltas[1:]:
        if delta_random_data >= delta_causal_data:
            num_as_good_deltas += 1

    counter = Counter(deltas)
    print("writing the results in `results/swap_randomisation.dat`")
    print("ALGORITHM    : Origo")
    print("ORIG_DELTA   : %.f" % delta_causal_data)
    print("FORMAT       : <delta> <frequency>")
    with open("results/swap_randomisation.dat", "w") as writer:
        for delta, freq in counter.items():
            writer.write("%f %d\n" % (delta, freq))

    p_value = num_as_good_deltas / num_rand_data
    print("P-VALUE      : %f" % p_value)


def test_compare_accuracy_symmetric(cardinality):
    prob_split = 1.0
    num_rows = 5000
    num_simulations = 500
    causes = list(range(cardinality))
    effects = list(range(cardinality, 2 * cardinality))
    max_height = cardinality + cardinality - 1
    dependencies = np.arange(0, 1.05, 0.1)

    ergo_accuracies = []
    origo_accuracies = []
    dc_accuracies = []

    progress(0, len(dependencies))

    for i, dependency in enumerate(dependencies):
        num_correct_ergo = 0
        num_correct_origo = 0
        num_correct_dc = 0
        num_samples = 0

        for _ in range(num_simulations):
            rows, ndeps = generate_data(num_rows, causes, effects,
                                        dependency, prob_split, max_height)
            if not ndeps:
                continue

            num_samples += 1
            store = Store(None, rows, causes + effects)
            ergo, origo, dc = score_ergo_origo_dc(causes, effects, store)

            num_correct_ergo += int(ergo.X_to_Y < ergo.Y_to_X)
            num_correct_origo += int(origo.X_to_Y < origo.Y_to_X)
            num_correct_dc += int(dc.X_to_Y < dc.Y_to_X)

        num_samples = num_samples or 1
        num_simulations = num_samples

        ergo_accuracies.append(num_correct_ergo / num_simulations)
        origo_accuracies.append(num_correct_origo / num_simulations)
        dc_accuracies.append(num_correct_dc / num_simulations)

        num_simulations = 500
        progress(i + 1, len(dependencies))

    print("\nwriting the results in `results/accuracy_symmetric.dat`")
    print("FORMAT: <dependency> <ergo_accuracy> <origo_accuracy> <dc_accuracy>")
    with open("results/accuracy_symmetric.dat", "w") as writer:
        for i, dependency in enumerate(dependencies):
            row = "%f %f %f %f" % (dependency, ergo_accuracies[i], origo_accuracies[
                                   i], dc_accuracies[i])
            writer.write(row + "\n")


def test_compare_accuracy_asymmetric(num_causes, num_effects):
    prob_split = 1.0
    num_rows = 5000
    num_simulations = 500
    max_height = num_causes + num_effects - 1
    dependencies = np.arange(0, 1.05, 0.1)

    ergo_accuracies = []
    origo_accuracies = []
    dc_accuracies = []

    progress(0, len(dependencies))

    for i, dependency in enumerate(dependencies):
        num_correct_ergo = 0
        num_correct_origo = 0
        num_correct_dc = 0
        num_samples = 0

        causes = list(range(num_causes))
        effects = list(range(num_causes, num_causes + num_effects))
        num_pairs_type_one = int(num_simulations / 2)
        for _ in range(num_pairs_type_one):
            rows, ndeps = generate_data(num_rows, causes, effects,
                                        dependency, prob_split, max_height)
            if not ndeps:
                continue

            num_samples += 1
            store = Store(None, rows, causes + effects)
            ergo, origo, dc = score_ergo_origo_dc(causes, effects, store)

            num_correct_ergo += int(ergo.X_to_Y < ergo.Y_to_X)
            num_correct_origo += int(origo.X_to_Y < origo.Y_to_X)
            num_correct_dc += int(dc.X_to_Y < dc.Y_to_X)

        causes = list(range(num_effects))
        effects = list(range(num_effects, num_causes + num_effects))
        for _ in range(num_simulations - num_pairs_type_one):
            rows, ndeps = generate_data(num_rows, causes, effects,
                                        dependency, prob_split, max_height)
            if not ndeps:
                continue

            num_samples += 1
            store = Store(None, rows, causes + effects)
            ergo, origo, dc = score_ergo_origo_dc(causes, effects, store)

            num_correct_ergo += int(ergo.X_to_Y < ergo.Y_to_X)
            num_correct_origo += int(origo.X_to_Y < origo.Y_to_X)
            num_correct_dc += int(dc.X_to_Y < dc.Y_to_X)

        num_samples = num_samples or 1
        num_simulations = num_samples

        ergo_accuracies.append(num_correct_ergo / num_simulations)
        origo_accuracies.append(num_correct_origo / num_simulations)
        dc_accuracies.append(num_correct_dc / num_simulations)

        num_simulations = 500
        progress(i + 1, len(dependencies))

    print("\nwriting the results in `results/accuracy_asymmetric.dat`")
    print("FORMAT: <dependency> <ergo_accuracy> <origo_accuracy> <dc_accuracy>")
    with open("results/accuracy_asymmetric.dat", "w") as writer:
        for i, dependency in enumerate(dependencies):
            row = "%f %f %f %f" % (dependency, ergo_accuracies[
                                   i], origo_accuracies[i], dc_accuracies[i])
            writer.write(row + "\n")


def test_split_probabilty_vs_accuracy():
    num_causes = 3
    num_effects = 3
    num_rows = 5000
    num_simulations = 500

    accuracies = []
    causes = list(range(num_causes))
    effects = list(range(num_causes, num_causes + num_effects))
    max_height = num_causes + num_effects - 1
    dependency = 1.0
    prob_splits = np.arange(0.0, 1.05, 0.1)

    progress(0, len(prob_splits))
    for i, prob_split in enumerate(prob_splits):
        num_correct = 0

        for _ in range(num_simulations):
            rows, ndeps = generate_data(num_rows, causes, effects,
                                        dependency, prob_split, max_height)
            store = Store(None, rows, causes + effects)
            _, origo, _ = score_ergo_origo_dc(causes, effects, store)
            num_correct += int(origo.X_to_Y < origo.Y_to_X)

        accuracies.append(num_correct / num_simulations)
        progress(i + 1, len(prob_splits))

    print()
    print(accuracies)
    print("writing the results in `results/accuracy_vs_prob_split.dat`")
    print("ALGORITHM    : Origo")
    print("FORMAT       : <dependency> <accuracy>")
    with open("results/accuracy_at_split_probability.dat", "w") as writer:
        for i, prob_split in enumerate(prob_splits):
            writer.write("%f %f\n" % (prob_split, accuracies[i]))


def test_dependency_vs_accuracy_at_heights():
    num_causes = 3
    num_effects = 3
    num_simulations = 500
    prob_split = 1.0
    num_rows = 2500

    accuracies = []
    causes = list(range(num_causes))
    effects = list(range(num_causes, num_causes + num_effects))
    dependencies = np.arange(0, 1.05, 0.1)

    for max_height in range(1, num_causes + num_effects):
        print("max_height =", max_height)
        progress(0, len(dependencies))
        accuracies_at_maxh = []

        for i, dependency in enumerate(dependencies):
            num_correct = 0
            for _ in range(num_simulations):
                rows, _ = generate_data(
                    num_rows, causes, effects, dependency, prob_split, max_height)
                store = Store(None, rows, causes + effects)
                _, origo, _ = score_ergo_origo_dc(causes, effects, store)
                num_correct += int(origo.X_to_Y < origo.Y_to_X)
            accuracy = num_correct / num_simulations
            accuracies_at_maxh.append(accuracy)
            progress(i + 1, len(dependencies))

        print()
        accuracies.append(accuracies_at_maxh)

    print()
    print(accuracies)

    print("writing the results in `results/accuracy_vs_dependency_at_heights.dat`")
    print("ALGORITHM    : Origo")
    print("FORMAT       : <dependency> <h0_accuracy> <h1_accuracy> ... <hn_accuracy>")
    with open("results/accuracy_at_heights.dat", "w") as writer:
        for i, dependency in enumerate(dependencies):
            row = [dependency]
            for accuracies_at_maxh in accuracies:
                row.append(accuracies_at_maxh[i])
            row = " ".join(str(item) for item in row)
            writer.write(row + "\n")


def test_dependency_vs_various_metrics():
    num_causes = 3
    num_effects = 3
    num_rows = 5000
    prob_split = 1.0
    num_simulations = 500
    dependencies = np.arange(0.0, 1.01, 0.1)
    max_height = num_causes + num_effects - 1
    causes = list(range(num_causes))
    effects = list(range(num_causes, num_causes + num_effects))

    stats = []
    progress(0, len(dependencies))
    for i, prob_dep in enumerate(dependencies):
        num_correct = 0
        num_indecisive = 0
        num_incorrect = 0

        for _ in range(num_simulations):
            rows, _ = generate_data(num_rows, causes, effects,
                                    prob_dep, prob_split, max_height)
            store = Store(None, rows, causes + effects)
            _, origo, _ = score_ergo_origo_dc(causes, effects, store)

            num_correct += int(origo.X_to_Y < origo.Y_to_X)
            num_indecisive += int(origo.X_to_Y == origo.Y_to_X)
            num_incorrect += int(origo.X_to_Y > origo.Y_to_X)

        percentage_correct = num_correct / num_simulations
        percentage_indecisive = num_indecisive / num_simulations
        percentage_incorrect = num_incorrect / num_simulations

        stats.append(
            (percentage_correct, percentage_indecisive, percentage_incorrect))
        progress(i + 1, len(dependencies))

    print()
    print(stats)
    print("writing the results in `results/dependency_vs_metrics.dat`")
    print("ALGORITHM    : Origo")
    print("FORMAT       : <dependency> <correct> <indecisive> <incorrect>")
    with open("results/dependency_vs_metrics.dat", "w") as writer:
        for i, stat in enumerate(stats):
            row = [dependencies[i], stat[0], stat[1], stat[2]]
            row = " ".join(str(item) for item in row)
            writer.write(row + "\n")


def test_scalability_symmetric_bagging():
    num_rows = 5000
    num_simulations = 500
    pairs = [(1, 1), (3, 3), (5, 5), (7, 7), (10, 10)]
    prob_split = 1.0
    dependency = 0.7
    epsilon = 0.0
    max_height = 5  # to control the tree generation time

    accuracy_pairs = []
    progress(0, len(pairs))

    for i, pair in enumerate(pairs):
        X = list(range(pair[0]))
        Y = list(range(pair[0], sum(pair)))

        num_correct_origo = 0
        num_correct_origob = 0
        for _ in range(num_simulations):
            rows, _ = generate_data(num_rows, X, Y, dependency,
                                    prob_split, max_height)
            store = Store(None, rows, X + Y)
            origo, origob = score_origo_bagging(X, Y, store)
            origob_XY = np.array(origob.X_to_Y)
            origob_YX = np.array(origob.Y_to_X)

            num_correct_origo += int(origo.Y_to_X - origo.X_to_Y > epsilon)
            num_correct_origob += int(sum(origob_XY < origob_YX)
                                      > sum(origob_XY > origob_YX))

        accuracy_origo = num_correct_origo / num_simulations
        accuracy_origob = num_correct_origob / num_simulations

        print((accuracy_origo, accuracy_origob))
        accuracy_pairs.append((accuracy_origo, accuracy_origob))

        progress(i + 1, len(pairs))

    print()
    print(accuracy_pairs)
    # print "writing the results in `results/scalability_symmetric.dat`"
    # print "FORMAT: <num_variables> <accuracy_dc> <accuracy_origi> <accuracy_origo>"
    # with open("results/scalability_symmetric.dat", "w") as writer:
    #     for i, pair in enumerate(pairs):
    #         row = [pair[0], accuracy_pairs[i][0],
    #                accuracy_pairs[i][1], accuracy_pairs[i][2]]
    #         row = " ".join(str(item) for item in row)
    #         writer.write(row + "\n")


def run_tests():
    print("1  test dependency vs. various metrics (origo)")
    print("2  test dependency vs. accuracy at various heights (origo)")
    print("3  test split probablity vs accuracy (origo)")
    print("4  test dependency vs. accuracy (asymmetric, compare all)")
    print("5  test dpeendency vs. accuracy (symmetric, compare all)")
    print("6  test swap randomisation (origo)")
    print("7  test statistical power (compare all)")
    print("8  test scalability (symmetric)")
    print("9  test scalability (asymmetric)")
    print("10 test scalability bagging (symmetric)")

    command = input("\nWhat would you wish of me? : ")
    command = command.strip()
    if command == "1":
        test_dependency_vs_various_metrics()
    elif command == "2":
        test_dependency_vs_accuracy_at_heights()
    elif command == "3":
        test_split_probabilty_vs_accuracy()
    elif command == "4":
        test_compare_accuracy_asymmetric(1, 3)
    elif command == "5":
        test_compare_accuracy_symmetric(3)
    elif command == "6":
        test_type1_error()
    elif command == "7":
        test_type2_error()
    elif command == "8":
        test_scalability_symmetric()
    elif command == "9":
        test_scalability_asymmetric()
    elif command == "10":
        test_scalability_symmetric_bagging()
    else:
        print("Invalid command.")


if __name__ == "__main__":
    run_tests()
