#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A Sequence is a container of type queue (FIFO) or stack (LIFO).
"""

class Sequence(object):

    def __init__(self):
        self.items = []

    def empty(self):
        return self.items == []

    def add(self, item):
        raise NotImplementedError()

    def get(self):
        return self.items.pop()

    def size(self):
        return len(self.items)


class Queue(Sequence):

    def __init__(self):
        super(Queue, self).__init__()

    def add(self, item):
        self.items.insert(0, item)


class Stack(Sequence):

    def __init__(self):
        super(Stack, self).__init__()

    def add(self, item):
        self.items.append(item)


if __name__ == "__main__":
    stack = Queue()
    stack.add(1)
    stack.add(2)
    stack.add(3)

    while not stack.empty():
        print(stack.get())
