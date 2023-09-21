import sys
import argparse
from time import time
from CSCPattern import Pattern, PatternConfig
#from CRPattern import Pattern, PatternConfig
from diagram import init_diagram, add_intent, add_object, clean_flags, print_lattice

if __name__ == "__main__":
    start_time = time()
    __parser__ = argparse.ArgumentParser(description='AddIntent')
    __parser__.add_argument('context_path', metavar='context_path', type=str, help='path to the formal context')
    __parser__.add_argument('-c', '--min_columns', metavar='min_columns', type=int, help='minimum number of columns', default=1)  # di paper icfca = theta
    __parser__.add_argument('-r', '--min_rows', metavar='min_rows', type=int, help='minimum number of rows', default=1)
    __parser__.add_argument('-n', '--number_of_concepts', metavar='number_of_concepts', type=int, help='max number of concepts', default=10000000)
    __parser__.add_argument('-t', '--time_limit', metavar='time_limit', type=int, help='time limit in seconds', default=3600)
    __args__ = __parser__.parse_args()

    #cfg = PatternConfig(min_col=__args__.min_columns, min_row=__args__.min_rows) # for CRPattern
    cfg = PatternConfig(theta=__args__.min_columns)

    with open(__args__.context_path, 'r') as f:
        L = init_diagram()
        #print_lattice(L, path + '_lattice.txt')
        print (L.nodes[-1])
        for object_id, line in enumerate(f):
            # print ('line', object_id, 'lattice size', len(L.nodes(data=True)), 'time', time() - start_time) ################
            raw_entry = line.replace('\n', '').replace('\r', '')
            pattern = Pattern(instance=raw_entry, config=cfg, object=object_id)
            #print ('partition', pattern)
            object_concept_id = add_intent(pattern, -1, L, 0, max, start_time)
            add_object(object_concept_id, object_id, L)
            break_time = time() - start_time
            if len(L.nodes(data=True)) > __args__.number_of_concepts:  # or break_time > __args__.time_limit:
                break
            clean_flags(L, object_concept_id)
    L = L.reverse()
    end_time = time()
    print_lattice(L, 'addintent_' + __args__.context_path + '_c' + str(__args__.min_columns) + '_r' + str(__args__.min_rows) + '_lattice.txt')
    print ('lattice size', len(L.nodes(data=True)))
    print(end_time - start_time)
