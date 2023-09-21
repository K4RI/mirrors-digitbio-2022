class PatternConfig:
    def __init__(self, min_col=1, min_row=1):
        self.min_col = min_col
        self.min_row = min_row


class Pattern:
    def __init__(self, instance=None, dirty=True, config=None, object=-1):
        self.objects = []
        self.partition = set()
        self.number_list = []
        if not config:
            self.cfg = PatternConfig(1)
        else:
            self.cfg = config

        if dirty:
            self.partition = self.pattern_parser(instance)
            self.objects.append(object)

    def pattern_parser(self, line):
        attribute_index = 0
        attribute_list = list()
        numbers = line.split()
        partition = set()
        for n in numbers:
            if n != '?':  # untuk ignore ?, jika ? dianggap missing value
                attribute_list.append((attribute_index, n))
                self.number_list.append(float(n))
            else:
                self.number_list.append(float(-1))
            attribute_index += 1
        attribute_list.sort(key=lambda tup: tup[1])

        partition_element = set()
        before = '-'
        for attribute in attribute_list:
            if attribute[1] != before and before != '-':
                if len(partition_element) >= self.cfg.min_col:
                    partition.add(frozenset(partition_element))
                partition_element = set()
            partition_element.add(attribute[0])
            before = attribute[1]
        if len(partition_element) >= self.cfg.min_col:
            partition.add(frozenset(partition_element))
        return partition

    def intersect(self, other):
        pi = Pattern(instance=None, config=self.cfg, dirty=False)
        pi.objects = self.objects + other.objects
        for element1 in self.partition:
            for element2 in other.partition:
                common = element1.intersection(element2)
                if len(common) >= self.cfg.min_col:
                    pi.partition.add(frozenset(common))
        return pi

    def partition_size_sum(self):
        size = 0
        for element in self.partition:
            size += len(element)
        return size

    def set_of_filled_columns(self):
        return_set = set()
        for n in range(len(self.number_list)):
            if self.number_list[n] > -1:
                return_set.add(n)
        return return_set

    def __eq__(self, other):
        if self.partition == other.partition:
            return True
        return False

    def __le__(self, other):
        '''pi = self.intersect(other)
        if pi == self:
            return True
        return False'''
        for i in self.partition:
            component_is_subsumed = False
            for j in other.partition:
                if i.issubset(j):
                    component_is_subsumed = True
                    break
            if not component_is_subsumed:
                return False
        return True

    def size(self):
        return len(self.objects)

    def __repr__(self):
        output = ''
        for partition_element in self.partition:
            for item in partition_element:
                output += str(item) + ' '
            output += '- '
        return output[:-2]
