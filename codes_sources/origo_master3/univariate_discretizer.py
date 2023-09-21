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


class UnivariateDiscretizer(object):

    def __init__(self, X, Y, max_bins_X=None, max_bins_Y=None, epsilon_X=None, delta_X=None, epsilon_Y=None, delta_Y=None, num_bins_X=None, num_bins_Y=None):
        self.X = X
        self.Y = Y
        self.max_bins_X = max_bins_X
        self.max_bins_Y = max_bins_Y
        self.num_bins_X = num_bins_X
        self.num_bins_Y = num_bins_Y
        self.epsilon_X = epsilon_X
        self.epsilon_Y = epsilon_Y
        self.delta_X = delta_X
        self.delta_Y = delta_Y

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


class UnivariateNMLDiscretizer(UnivariateDiscretizer):

    def get_cut_points(self, data, max_bins, epsilon, delta):
        # write the data into the file
        unrolled_data = " ".join(str(observation) for observation in data)
        nml_data_path = os.path.join("nml", "data.dat")
        with open(nml_data_path, "w") as fp:
            fp.write(unrolled_data)

        # run NML histogram tool
        nml_hist_exec = os.path.join("nml", "NML_histogram")
        # ./NML_histogram iris2.dat 10 0.1 0.2
        nml_hist_cmd = [nml_hist_exec, nml_data_path,
                        str(max_bins), str(epsilon), str(delta)]
        output = subprocess.check_output(nml_hist_cmd)

        # extract the cut points using regex from the output byte string
        cut_points = PATTERN_CUT_POINTS.findall(output)[0].split(",")
        cut_points = list(map(float, cut_points))

        return cut_points

    def discretize(self):
        cut_points_X = self.get_cut_points(
            self.X, self.max_bins_X, self.epsilon_X, self.delta_X)
        cut_points_Y = self.get_cut_points(
            self.Y, self.max_bins_Y, self.epsilon_Y, self.delta_Y)

        num_bins_X = len(cut_points_X) - 1
        num_bins_Y = len(cut_points_Y) - 1

        alphabets_X = list(range(num_bins_X))
        alphabets_Y = list(range(num_bins_X, num_bins_X + num_bins_Y))

        discretized_X = self.map_alphabets(self.X, cut_points_X, alphabets_X)
        discretized_Y = self.map_alphabets(self.Y, cut_points_Y, alphabets_Y)

        return alphabets_X, discretized_X, alphabets_Y, discretized_Y


class UnivariateEquiWidthDiscretizer(UnivariateDiscretizer):

    def discretize(self):
        min_X, max_X = min(self.X), max(self.X)
        min_Y, max_Y = min(self.Y), max(self.Y)

        interval_X = (max_X - min_X) * 1.0 / self.num_bins_X
        interval_Y = (max_Y - min_Y) * 1.0 / self.num_bins_Y

        cut_points_X = np.append(np.arange(min_X, max_X, interval_X), [max_X])
        cut_points_Y = np.append(np.arange(min_Y, max_Y, interval_Y), [max_Y])

        alphabets_X = list(range(self.num_bins_X))
        alphabets_Y = list(range(self.num_bins_X, self.num_bins_X + self.num_bins_Y))

        discretized_X = self.map_alphabets(self.X, cut_points_X, alphabets_X)
        discretized_Y = self.map_alphabets(self.Y, cut_points_Y, alphabets_Y)

        return alphabets_X, discretized_X, alphabets_Y, discretized_Y


class UnivariateEquiFrequencyDiscretizer(UnivariateDiscretizer):

    def get_cut_points(self, sorted_data, num_bins):
        bin_freq = len(sorted_data) / num_bins
        cut_points = [sorted_data[0]]
        for i in range(1, num_bins):
            cut_point = sorted_data[i * bin_freq]
            cut_points.append(cut_point)
        cut_points.append(sorted_data[-1])
        return cut_points

    def discretize(self):
        sorted_X = sorted(self.X)
        sorted_Y = sorted(self.Y)

        cut_points_X = self.get_cut_points(sorted_X, self.num_bins_X)
        cut_points_Y = self.get_cut_points(sorted_Y, self.num_bins_Y)

        alphabets_X = list(range(self.num_bins_X))
        alphabets_Y = list(range(self.num_bins_X, self.num_bins_X + self.num_bins_Y))

        discretized_X = self.map_alphabets(self.X, cut_points_X, alphabets_X)
        discretized_Y = self.map_alphabets(self.Y, cut_points_Y, alphabets_Y)

        return alphabets_X, discretized_X, alphabets_Y, discretized_Y


class UnivariateIPDiscretizer(UnivariateDiscretizer):

    def run_ipd(self):
        self.output_path = os.path.join("ipd", "out.txt")
        runtime_path = os.path.join("ipd", "runtime.txt")
        cutpoint_path = os.path.join("ipd", "cuts.txt")
        exec_path = os.path.join("ipd", "jar", "ipd.jar")
        ipd_cmd = ["java", "-jar", exec_path, "-FILE_INPUT", self.data_path, "-FILE_CP_OUTPUT", cutpoint_path, "-FILE_RUNTIME_OUTPUT", runtime_path, "-FILE_DATA_OUTPUT",
                   self.output_path, "-NUM_ROWS", str(len(self.X)), "-NUM_MEASURE_COLS", "2", "-NUM_CAT_CONTEXT_COLS", "0", "-MAX_VAL", str(max(max(self.X), max(self.Y))), "-METHOD", "0"]
        output = subprocess.check_output(ipd_cmd)

    def read_ipd_output(self):
        with open(self.output_path, "r") as fp:
            raw_data = fp.read()

        minus_one = lambda x: int(x) - 1
        alphabets_X = list(map(minus_one, re.findall(
            "dim0 \{([\d+,?]+)\}", raw_data)[0].split(",")))
        alphabets_Y = list(map(minus_one, re.findall(
            "dim1 \{([\d+,?]+)\}", raw_data)[0].split(",")))
        data = re.findall("\d+,\d+,\"0\"", raw_data)

        discretized_X = []
        discretized_Y = []
        for row in data:
            row = list(map(minus_one, row.split(",")[:-1]))
            discretized_X.append(row[0])
            discretized_Y.append(row[1])

        return alphabets_X, discretized_X, alphabets_Y, discretized_Y

    def write_ipd_input(self):
        # IPD expects a class label in the end, we provide a dummy one
        input_rows = []
        for i in range(len(self.X)):
            input_rows.append("%f;%f;0" % (self.X[i], self.Y[i]))
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
    X = list(map(abs, np.random.normal(mu, sigma, 100)))
    Y = list(map(abs, np.random.normal(mu, sigma, 100)))
    print(X, Y)
    # X = [random.randint(1, 20) for i in xrange(100)]
    # Y = [random.randint(1, 20) for i in xrange(100)]
    # X = range(1, 13)
    # Y = range(1, 13)
    print("X =", X)
    print("Y =", Y)
    print()
    discretizer = UnivariateNMLDiscretizer(
        X, Y, max_bins_X=3, max_bins_Y=3, epsilon_X=0.1, delta_X=0.2, epsilon_Y=0.1, delta_Y=0.2)
    aX, Xd, aY, Yd = discretizer.discretize()
    print(aX, aY)
    print(Xd)
    print(Yd)
    print("\n" * 2)

    discretizer = UnivariateEquiWidthDiscretizer(
        X, Y, num_bins_X=3, num_bins_Y=3)
    aX, Xd, aY, Yd = discretizer.discretize()
    print(aX, aY)
    print(Xd)
    print(Yd)
    print("\n" * 2)

    discretizer = UnivariateEquiFrequencyDiscretizer(
        X, Y, num_bins_X=4, num_bins_Y=4)
    aX, Xd, aY, Yd = discretizer.discretize()
    print(aX, aY)
    print(Xd)
    print(Yd)
    print("\n" * 2)


if __name__ == "__main__":
    mu, sigma = 0, 0.1
    X = list(range(100))
    Y = list(range(100, 200))
    discretizer = UnivariateIPDiscretizer(X, Y)
    aX, Xd, aY, Yd = discretizer.discretize()
    print(aX, aY)
    print(Xd)
    print(Yd)
    print("\n" * 2)

    # test_discretizer()
