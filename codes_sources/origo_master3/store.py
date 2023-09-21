#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Column-store data
"""
import random


class Store(object):

    def __init__(self, file_path, rows=None, column_names=None):
        # We read column names from the data by default. As a result, if there
        # is a zero column vector in the data, we will miss that column. That
        # is the reason why we take column_names as extra argument.
        self.rows = []
        self.columns = []
        self.supports = []
        self.column_names = []
        self.column_pos = {}
        self.num_rows = 0
        self.load(file_path, rows, column_names)

    @staticmethod
    def combine(left_rows, right_rows):
        assert(len(left_rows) == len(right_rows))
        new_rows = []
        for i in range(len(left_rows)):
            row = []
            row.extend(left_rows[i])
            row.extend(right_rows[i])
            new_rows.append(row)
        return new_rows

    def extract_columns(self):
        self.supports = [0] * len(self.column_names)
        for row in self.rows:
            for column_name in self.column_names:
                item = column_name if column_name in row else None
                pos = self.column_pos[column_name]
                self.columns[pos].append(item)
                self.supports[pos] += int(item != None)

    def extract_metadata(self, column_names=None):
        alphabet = set()
        for row in self.rows:
            alphabet.update(row)

        self.num_rows = len(self.rows)
        self.column_names = list(
            sorted(alphabet)) if not column_names else column_names
        for idx, column_name in enumerate(self.column_names):
            self.column_pos[column_name] = idx
            self.columns.append([])

    def flip_value(self, row, col, cur_val, col_name):
        flipped_val = col_name if cur_val is None else None
        self.set_value(row, col, flipped_val)

    def get_column(self, column_name):
        col_pos = self.column_pos[column_name]
        column = self.columns[col_pos]
        return column

    def get_column_binarized(self, column_name):
        column = self.get_column(column_name)
        return [int(col == column_name) for col in column]

    def get_column_names(self):
        return self.column_names

    def get_column_margins(self):
        col_margins = []
        for col in self.columns:
            margin = 0
            for item in col:
                if item != None:
                    margin += 1
            col_margins.append(margin)
        return col_margins

    def get_support(self, column_name):
        pos = self.column_pos[column_name]
        return self.supports[pos]

    def get_supports(self):
        return self.supports

    def get_num_ones(self):
        return sum(self.get_column_margins())

    def get_rows(self):
        return self.rows

    def get_row_margins(self):
        row_margins = []
        for row_num in range(self.num_rows):
            margin = 0
            for col in self.columns:
                if col[row_num] != None:
                    margin += 1
            row_margins.append(margin)
        return row_margins

    def get_value(self, row, col):
        return self.columns[col][row]

    def load(self, file_path, rows, column_names):
        self.rows = rows if rows else Store.read_rows(file_path)
        self.extract_metadata(column_names)
        self.extract_columns()

    def project(self, col_names):
        projection = []
        for i in range(self.num_rows):
            items = [self.get_column(col_name)[i]
                     for col_name in col_names]
            row = [item for item in items if item != None]
            projection.append(row)
        return projection

    @staticmethod
    def read_rows(file_path):
        rows = []
        with open(file_path, 'r') as fp:
            for line in fp:
                row = list(map(int, line.split()))
                rows.append(row)
        return rows

    def set_value(self, row, col, val):
        self.columns[col][row] = val

    def swap(self, row1, row2, col1, col2, col_name1, col_name2):
        r1c1 = self.get_value(row1, col1)
        r1c2 = self.get_value(row1, col2)
        r2c1 = self.get_value(row2, col1)
        r2c2 = self.get_value(row2, col2)

        r1 = int(r1c1 == None) + int(r1c2 == None)
        r2 = int(r2c1 == None) + int(r2c2 == None)
        c1 = int(r1c1 == None) + int(r2c1 == None)
        c2 = int(r1c2 == None) + int(r2c2 == None)

        if r1 == r2 == c1 == c2 == 1:
            self.flip_value(row1, col1, r1c1, col_name1)
            self.flip_value(row1, col2, r1c2, col_name2)
            self.flip_value(row2, col1, r2c1, col_name1)
            self.flip_value(row2, col2, r2c2, col_name2)
            return True
        return False

    def swap_randomise(self):
        # swap randomisation following Gionis et al. 2005
        rows = list(range(self.num_rows))
        num_ones = self.get_num_ones()
        num_swap_operation = 0
        while num_swap_operation < num_ones:
            row1, row2 = random.sample(rows, 2)
            col_name1, col_name2 = random.sample(self.column_names, 2)
            col1 = self.column_pos[col_name1]
            col2 = self.column_pos[col_name2]
            if self.swap(row1, row2, col1, col2, col_name1, col_name2):
                num_swap_operation += 1


if __name__ == "__main__":
    #store = Store(
    #    '/Users/kbudhath/Projects/dc/data/datasets/dat-files/wine.dat')
    store = Store(
        './data/adult/adult.dat')
    print(store.project([1, 3]))
