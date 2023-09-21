import sys
from time import time
from CSCPattern import Pattern, PatternConfig
#from CRPattern import Pattern, PatternConfig
import argparse

patterns = []
L = []
G = set()
f_out = 0
concept_count = 0
max_concept_count = 0

def process(A, g, concept, start_time):  # see thesis of Kaytoue
    #print 'process', A
    global concept_count
    global max_concept_count
    if concept_count == max_concept_count and max_concept_count > 0:
        return
    C_minus_A = set(concept.objects) - A
    less_than_g = set(filter(lambda x: x < g, C_minus_A))
    if len(less_than_g) == 0:
        # if len(concept.objects) >= concept.cfg.min_row:
        if len(concept.objects) >= concept.cfg.theta:
            f_out.write(str(concept) + '|' + str(len(concept.objects)) + '|' + str(concept.objects) + '\n')
            concept_count += 1
            if concept_count == 10 or concept_count == 100 or concept_count == 500:
                f_out.write('check point ' + str(concept_count) + ' ' + str(time() - start_time) + '\n')
        G_minus_C = G - set(concept.objects)
        more_than_g = set(filter(lambda l: l > g, G_minus_C))
        for f in more_than_g:
            z = set(concept.objects)
            z.add(f)
            new_candidate = concept.intersect(patterns[f])
            # if len(z) >= new_candidate.cfg.min_row:  # kalau len(extent) kurang dari min_row, ngitung closure nya nanti aja, setelah nemu len(extent) yang >= minrow
            if len(z) >= new_candidate.cfg.theta:
                x = closure(new_candidate, z)
            else:
                x = list(z)
            new_candidate.objects = x
            process(z, f, new_candidate, start_time)
    else:
        #print 'backtrack'
        pass
    return


def closure(input_pattern, no_check):
    #print 'closure', no_check
    time1 = time()
    objects = list(no_check)
    check = G - no_check
    for c in check:
        if input_pattern <= patterns[c]:
            objects.append(patterns[c].objects[0])
    '''
    for p in patterns:
        if input_pattern <= p:
            objects.append(p.objects[0])'''
    time2 = time()
    # print len(check), 'compute closure', (time2 - time1)
    return objects


def alignment(base_row, input_row, base_column_index, type):
    base_cells = base_row.split()
    input_cells = input_row.split()
    difference_a = float(input_cells[base_column_index]) - float(base_cells[base_column_index])
    if base_cells[base_column_index] == '0':
        difference_m = float(input_cells[base_column_index]) / 0.1
    else:
        difference_m = float(input_cells[base_column_index]) / float(base_cells[base_column_index])

    if difference_m == 0:
        difference_m = 0.1

    out = ''
    for cell in input_cells:
        if type == 'a':
            out += str(float(cell) - difference_a) + '\t'
        elif type == 'm':
            out += str(round(float(cell) / difference_m, 1)) + '\t'
    return out[:-1]


if __name__ == "__main__":
    __parser__ = argparse.ArgumentParser(description='CBO')
    __parser__.add_argument('context_path', metavar='context_path', type=str, help='path to the formal context')
    __parser__.add_argument('-c', '--min_columns', metavar='min_columns', type=int, help='minimum number of columns', default=1)
    __parser__.add_argument('-r', '--min_rows', metavar='min_rows', type=int, help='minimum number of rows', default=1)
    __parser__.add_argument('-n', '--number_of_concepts', metavar='number_of_concepts', type=int, help='max number of concepts', default=0)
    __args__ = __parser__.parse_args()

    start_time = time()
    # cfg = PatternConfig(min_col=__args__.min_columns, min_row=__args__.min_rows)
    cfg = PatternConfig(theta=__args__.min_columns)
    f_out = open('cbo_' + __args__.context_path + '_c' + str(__args__.min_columns) + '_r' + str(__args__.min_rows) + '.txt', 'w')

    with open(__args__.context_path, 'r') as f_in:
        for object_id, line in enumerate(f_in):
            #print 'line', object_id
            raw_entry = line.replace('\n', '').replace('\r', '')
            patterns.append(Pattern(instance=raw_entry, config=cfg, object=object_id))
            G.add(object_id)

    max_concept_count = __args__.number_of_concepts
    for i in patterns:
        print ('iter', i.objects[0])
        if len(i.partition) == 0:
            continue
        start_iter = time()
        if concept_count == max_concept_count and max_concept_count > 0:
            break
        candidate_concept = Pattern(dirty=False, config=cfg)
        candidate_concept.partition = list(i.partition)
        candidate_concept.objects = closure(i, set(i.objects))
        process(set(i.objects), i.objects[0], candidate_concept, start_time)
        print(time() - start_iter)

    end_time = time()
    print(end_time - start_time)
