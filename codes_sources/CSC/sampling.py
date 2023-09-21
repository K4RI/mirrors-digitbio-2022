import random
import sys
from time import time
# from CSCPattern import Pattern, PatternConfig
from CRPattern import Pattern, PatternConfig
from diagram import init_diagram, add_intent, add_object, clean_flags, print_lattice


def get_objects_by_partition(pattern_set, wanted_pattern):
    objects = list()
    for pattern in pattern_set:
        if wanted_pattern <= pattern:
            objects.append(pattern.objects[0])
    return objects


def get_partition_by_objects(pattern_set, objects):
    result = pattern_set[objects[0]]
    for j in range(1, len(objects)):
        result = result.intersect(pattern_set[objects[j]])
    return result


if __name__ == "__main__":
    start_time = time()
    path = sys.argv[1]
    number_of_iterations = int(sys.argv[2])

    cfg = PatternConfig(1)
    patterns = list()
    with open(path, 'r') as f:
        for object_id, line in enumerate(f):
            raw_entry = line.replace('\n', '').replace('\r', '')
            partition = Pattern(instance=raw_entry, config=cfg, object=object_id)
            patterns.append(partition)
            if object_id == 0:
                number_of_attributes = len(raw_entry.split())
            # print 'line', object_id, '=', pattern

    all_objects = [x for x in range(0, len(patterns))]
    top_intent = get_partition_by_objects(patterns, all_objects)
    dummy_object = '1 '*number_of_attributes
    bottom_intent = Pattern(instance=dummy_object[:-1], dirty=True)

    top_or_bottom = random.choice(['T', 'B'])
    top_or_bottom = 'T'
    current_extent = []
    if top_or_bottom == 'T':
        current_intent = top_intent
        current_extent = list(all_objects)
    else:
        current_intent = bottom_intent
        current_extent = get_objects_by_partition(patterns, current_intent)

    for i in range(0, number_of_iterations):
        print ('current_extent', current_extent)
        print ('current_intent', current_intent)
        if current_intent == top_intent:
            direction = 'down'
        elif current_intent == bottom_intent:
            direction = 'up'
        else:
            direction = random.choice(['up', 'down'])
            direction = 'down'

        if direction == 'up':
            o = random.choice(list(set(all_objects) - set(current_extent)))
            temp_extent = list(current_extent)
            temp_extent.append(o)
            candidate_intent = get_partition_by_objects(patterns, temp_extent)
            candidate_extent = get_objects_by_partition(patterns, candidate_intent)
            alpha = float((len(candidate_extent) - len(current_extent)) * number_of_attributes) / float(((len(candidate_intent.partition) - len(current_intent.partition)) * len(all_objects)))
            print ('o', o)
            print ('candidate_extent', candidate_extent)
            print ('candidate_intent', candidate_intent)
            print ('alpha', alpha)
        elif direction == 'down':
            temp_intent = Pattern(dirty=False)
            temp_intent.partition = set(current_intent.partition)
            to_be_joined = random.sample(tuple(temp_intent.partition), 2)
            temp_intent.partition = temp_intent.partition.difference(set(to_be_joined))
            temp_intent.partition.add(to_be_joined[0].union(to_be_joined[1]))
            candidate_extent = get_objects_by_partition(patterns, temp_intent)
            candidate_intent = get_partition_by_objects(patterns, candidate_extent)
            alpha = float(((len(current_intent.partition) - len(candidate_intent.partition)) * len(all_objects))) / float((len(current_extent) - len(candidate_extent)) * number_of_attributes)
            print ('temp_intent', temp_intent)
            print ('candidate_extent', candidate_extent)
            print ('candidate_intent', candidate_intent)
            print ('alpha', alpha)

        x = random.choice([0, 1])
        if x < alpha:
            current_extent = candidate_extent
            current_intent = candidate_intent

    end_time = time()
    print(end_time - start_time)
