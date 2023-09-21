#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""discretizes the continuous real-valued data using one of the discretization strategies.
"""
import os
import random
import re
import subprocess

import numpy as np


PATTERN_CUT_POINTS = re.compile(r"Cut points:\n(.*)")


class MvDiscretizer(object):

    def __init__(self, X, Y, max_bins=None, epsilon=None, delta=None, num_bins=None):
        self.X = X  # list of random variables
        self.Y = Y
        self.max_bins = max_bins
        self.epsilon = epsilon
        self.delta = delta
        self.num_bins = num_bins

    def map_alphabets(self, data, cut_points, alphabets):
        # first bin corresponds to the first alphabet, 2nd to the second, etc.
        # each interval is close-open interval [x,y)
        cut_points[-1] += 1  # to close the last bin with max value
        intervals = [(cut_points[i], cut_points[i + 1])
                     for i in range(len(cut_points) - 1)]
        discretized_data = []
        for point in data:
            for i, interval in enumerate(intervals):
                if interval[0] <= point < interval[1]:
                    discretized_data.append(alphabets[i])
                    break
        assert len(discretized_data) == len(data)
        return discretized_data

    def discretize(self):
        raise NotImplementedError("Subclass must implement abstract method.")


class NMLMvDiscretizer(MvDiscretizer):

    def get_cut_points(self, data, max_bins, epsilon, delta):
        # write the data into the file
        unrolled_data = " ".join(str(observation) for observation in data)
        nml_data_path = os.path.join("nml", "data.dat")
        with open(nml_data_path, "w") as fp:
            fp.write(unrolled_data)

        # run NML histogram tool
        nml_hist_exec = os.path.join("nml", "NML_histogram.exe")
        # ./NML_histogram iris2.dat 10 0.1 0.2
        nml_hist_cmd = [nml_hist_exec, nml_data_path,
                        str(max_bins), str(epsilon), str(delta)]
        output = subprocess.check_output(nml_hist_cmd)
        output = output.decode("utf-8").replace('\r\n', '\n')

        # extract the cut points using regex from the output byte string
        cut_points = PATTERN_CUT_POINTS.findall(output)[0].split(",")
        cut_points = list(map(float, cut_points))

        return cut_points

    def discretize(self):
        cut_points_X = [self.get_cut_points(
            data, self.max_bins, self.epsilon, self.delta) for data in self.X]
        cut_points_Y = [self.get_cut_points(
            data, self.max_bins, self.epsilon, self.delta) for data in self.Y]

        num_bins_X = [len(cut_points) - 1 for cut_points in cut_points_X]
        num_bins_Y = [len(cut_points) - 1 for cut_points in cut_points_Y]

        pool_X = list(range(sum(num_bins_X)))
        pool_Y = list(range(sum(num_bins_X), sum(num_bins_X) + sum(num_bins_Y)))

        alphabets_X = [pool_X[sum(num_bins_X[:i]):sum(
            num_bins_X[:i]) + n] for i, n in enumerate(num_bins_X)]
        alphabets_Y = [pool_Y[sum(num_bins_Y[:i]):sum(
            num_bins_Y[:i]) + n] for i, n in enumerate(num_bins_Y)]

        discretized_X = [self.map_alphabets(data, cut_points_X[i], alphabets_X[
                                            i]) for i, data in enumerate(self.X)]
        discretized_Y = [self.map_alphabets(data, cut_points_Y[i], alphabets_Y[
                                            i]) for i, data in enumerate(self.Y)]

        return alphabets_X, discretized_X, alphabets_Y, discretized_Y


class EquiWidthMvDiscretizer(MvDiscretizer):

    def discretize_individual(self, data, data_num):
        min_val, max_val = min(data), max(data)
        interval = (max_val - min_val) * 1.0 / self.num_bins
        cut_points = np.append(
            np.arange(min_val, max_val, interval), [max_val])
        alphabets = list(range(data_num * self.num_bins,
                          data_num * self.num_bins + self.num_bins))
        discretized = self.map_alphabets(data, cut_points, alphabets)
        return alphabets, discretized

    def discretize(self):
        alphabets_X = []
        discretized_X = []
        data_num = 0
        for data in self.X:
            alphabets, discretized = self.discretize_individual(data, data_num)
            alphabets_X.append(alphabets)
            discretized_X.append(discretized)
            data_num += 1

        alphabets_Y = []
        discretized_Y = []
        for data in self.Y:
            alphabets, discretized = self.discretize_individual(data, data_num)
            alphabets_Y.append(alphabets)
            discretized_Y.append(discretized)
            data_num += 1

        return alphabets_X, discretized_X, alphabets_Y, discretized_Y


class EquiFrequencyMvDiscretizer(MvDiscretizer):

    def get_cut_points(self, sorted_data, num_bins):
        bin_freq = len(sorted_data) // num_bins
        cut_points = [sorted_data[0]]
        for i in range(1, num_bins):
            cut_point = sorted_data[i * bin_freq]
            cut_points.append(cut_point)
        cut_points.append(sorted_data[-1])
        return cut_points

    def discretize_individual(self, data, data_num):
        sorted_data = sorted(data)
        cut_points = self.get_cut_points(sorted_data, self.num_bins)
        alphabets = list(range(data_num * self.num_bins,
                          data_num * self.num_bins + self.num_bins))
        discretized = self.map_alphabets(data, cut_points, alphabets)
        return alphabets, discretized

    def discretize(self):
        alphabets_X = []
        discretized_X = []
        data_num = 0
        for data in self.X:
            alphabets, discretized = self.discretize_individual(data, data_num)
            alphabets_X.append(alphabets)
            discretized_X.append(discretized)
            data_num += 1

        alphabets_Y = []
        discretized_Y = []
        for data in self.Y:
            alphabets, discretized = self.discretize_individual(data, data_num)
            alphabets_Y.append(alphabets)
            discretized_Y.append(discretized)
            data_num += 1

        return alphabets_X, discretized_X, alphabets_Y, discretized_Y


class IPMvDiscretizer(MvDiscretizer):

    def run_ipd(self):
        max_val = max(max(data) for data in self.X + self.Y)
        self.output_path = os.path.join("ipd", "out.txt")
        runtime_path = os.path.join("ipd", "runtime.txt")
        cutpoint_path = os.path.join("ipd", "cuts.txt")
        exec_path = os.path.join("ipd", "jar", "ipd.jar")
        ipd_cmd = ["java", "-jar", exec_path, "-FILE_INPUT", self.data_path, "-FILE_CP_OUTPUT", cutpoint_path, "-FILE_RUNTIME_OUTPUT", runtime_path, "-FILE_DATA_OUTPUT",
                   self.output_path, "-NUM_ROWS", str(len(self.X[0])), "-NUM_MEASURE_COLS", str(len(self.X) + len(self.Y)), "-NUM_CAT_CONTEXT_COLS", "0", "-MAX_VAL", str(max_val), "-METHOD", "0"]
        output = subprocess.check_output(ipd_cmd)

    def read_ipd_output(self):
        with open(self.output_path, "r") as fp:
            raw_data = fp.read()

        minus_one = lambda x: int(x) - 1
        # read all alphabets
        alphabets_X = [list(map(minus_one, re.findall(
            "dim%d \{([\d+,?]+)\}" % i, raw_data)[0].split(","))) for i in range(len(self.X))]
        alphabets_Y = [list(map(minus_one, re.findall(
            "dim%d \{([\d+,?]+)\}" % i, raw_data)[0].split(","))) for i in range(len(self.X), len(self.X) + len(self.Y))]
        data = re.findall("[\d+,]+\"0\"", raw_data)

        discretized_X = [list() for x in self.X]
        discretized_Y = [list() for y in self.Y]
        for row in data:
            row = list(map(minus_one, row.split(",")[:-1]))
            j = 0
            for i in range(len(self.X)):
                discretized_X[i].append(row[j])
                j += 1

            for i in range(len(self.Y)):
                discretized_Y[i].append(row[j])
                j += 1

        return alphabets_X, discretized_X, alphabets_Y, discretized_Y

    def write_ipd_input(self):
        # IPD expects a class label in the end, we provide a dummy one
        input_rows = []
        num_rows = len(self.X[0])
        for i in range(num_rows):
            row = ";".join(str(self.X[j][i]) for j in range(len(self.X)))
            row += ";" + ";".join(str(self.Y[j][i])
                                  for j in range(len(self.Y))) + ";0"
            input_rows.append(row)
        input_raw = "\n".join(row for row in input_rows)
        self.data_path = os.path.join("ipd", "data.csv")
        with open(self.data_path, "w") as fp:
            fp.write(input_raw)

    def discretize(self):
        self.write_ipd_input()
        self.run_ipd()
        return self.read_ipd_output()


def test_discretizer():
    mu, sigma = 0, 0.1
    X = [list(map(abs, np.random.normal(mu, sigma, 100)))]
    Y = [list(map(abs, np.random.normal(mu, sigma, 100)))]
    X = [[random.randint(1, 20) for i in range(100)]]
    Y = [[random.randint(1, 20) for i in range(100)]]
    # X = [range(1, 13)]
    # Y = [range(1, 13)]
    print(X, Y)
    print()
    discretizer = NMLMvDiscretizer(
        X, Y, max_bins=3, epsilon=0.1, delta=0.2)

    aX, Xd, aY, Yd = discretizer.discretize()
    print(aX, aY)
    print(Xd)
    print(Yd)
    print("\n" * 2)

    discretizer = EquiWidthMvDiscretizer(
        X, Y, num_bins=3)
    aX, Xd, aY, Yd = discretizer.discretize()
    print(aX, aY)
    print(Xd)
    print(Yd)
    print("\n" * 2)

    discretizer = EquiFrequencyMvDiscretizer(
        X, Y, num_bins=4)
    aX, Xd, aY, Yd = discretizer.discretize()
    print(aX, aY)
    print(Xd)
    print(Yd)
    print("\n" * 2)


if __name__ == "__main__":
    mu, sigma = 0, 0.1
    X = [list(range(100)), list(range(100, 200))]
    Y = [list(range(200, 300)), list(range(300, 400))]
    discretizer = IPMvDiscretizer(X, Y)
    # discretizer = NMLMvDiscretizer(X, Y, max_bins=4, epsilon=1, delta=2)
    aX, Xd, aY, Yd = discretizer.discretize()
    print(aX, aY)
    print(Xd)
    print(Yd)
    print("\n" * 2)

    test_discretizer()
