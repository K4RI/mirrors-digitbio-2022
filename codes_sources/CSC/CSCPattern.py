import itertools


def pattern_parser(line, begin):
    signs = line.split(' ')
    element = set()
    i = begin
    for sign in signs:
        if sign == '-' or sign == '+':
            element.add(Item(i, sign))
        i += 1
    return element


class Item:
    def __init__(self, attribute, sign):
        self.attribute = str(attribute)
        self.sign = str(sign)

    def opposite_sign(self):
        if self.sign == '+':
            return '-'
        return '+'

    def __eq__(self, other):
        if self.attribute == other.attribute and self.sign == other.sign:
            return True
        return False

    def __hash__(self):
        return hash((self.attribute, self.sign))

    def __repr__(self):
        output = str(self.attribute) + str(self.sign)
        return output


class PartitionElement:  # partition element = partition component
    def __init__(self, items):
        self.set_of_items = frozenset(items)

    def __eq__(self, other):
        if self.set_of_items == other.set_of_items:
            return True
        # check if all signs are opposite
        equal = True
        self_switched = set()
        for item in self.set_of_items:
            self_switched.add(Item(item.attribute, item.opposite_sign()))
        if self_switched == other.set_of_items:
            return True
        return False

    def __hash__(self):  # an attribute can not appear in more than one PartitionElement
        set_without_sign = set()
        for item in self.set_of_items:
            set_without_sign.add(item.attribute)
        return hash(frozenset(set_without_sign))

    def __repr__(self):
        output = ''
        for item in self.set_of_items:
            output += str(item) + ' '
        return output[:-1]

    def size(self):
        return len(self.set_of_items)


class PatternConfig:
    def __init__(self, theta=1):
        self.theta = theta


class Pattern:
    def __init__(self, instance=None, dirty=True, config=None, object=-1):
        self.objects = []
        self.partition = set()
        if not config:
            self.cfg = PatternConfig(1)
        else:
            self.cfg = config

        if dirty:
            self.partition.add(PartitionElement(pattern_parser(instance, 1)))
            self.objects.append(object)
        #else:
            #self.partition = set()

    def intersect(self, other):
        pi = Pattern(instance=None, config=self.cfg, dirty=False)
        for element1 in self.partition:
            for element2 in other.partition:
                self.common(element1, element2, pi.partition)

        ignore = True
        for element in pi.partition:
            if element.size() >= self.cfg.theta:
                ignore = False
                break
        if ignore:
            pi.partition = set()
        return pi

    def common(self, el1, el2, partition):
        intersect_same = el1.set_of_items.intersection(el2.set_of_items)
        if len(intersect_same) > 0:
            partition.add(PartitionElement(intersect_same))

        remaining = el1.set_of_items.difference(intersect_same)
        remaining_switched = set()
        for item in remaining:
            remaining_switched.add(Item(item.attribute, item.opposite_sign()))
        intersect_opposite = remaining_switched.intersection(el2.set_of_items)
        if len(intersect_opposite) > 0:
            partition.add(PartitionElement(intersect_opposite))

    def __eq__(self, other):
        if self.partition == other.partition:
            return True
        return False

    def __le__(self, other):
        pi = self.intersect(other)
        if pi == self:
            return True
        return False

    def size(self):
        return len(self.objects)

    def __repr__(self):
        output = ''
        for element in self.partition:
            output += str(element) +  '|'
        return output[:-1]
