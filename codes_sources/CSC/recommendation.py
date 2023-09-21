import sys
from operator import itemgetter
from time import time
from CRPattern import Pattern, PatternConfig
import argparse

patterns = list()
patterns.append(0)

if __name__ == "__main__":
    __parser__ = argparse.ArgumentParser(description='recommendation')
    __parser__.add_argument('-context', '--context_path', metavar='context_path', type=str, help='path to the formal context')
    __parser__.add_argument('-test', '--test', metavar='test', type=str, help='path to the true test value')
    __parser__.add_argument('-c', '--min_columns', metavar='min_columns', type=int, help='minimum number of columns', default=1)
    __parser__.add_argument('-r', '--min_rows', metavar='min_rows', type=int, help='minimum number of rows', default=1)
    __parser__.add_argument('-type', '--context_type', metavar='context_type', type=str, help='o=original, t=transpose', default='o')
    __args__ = __parser__.parse_args()
    start_time = time()
    cfg = PatternConfig(min_col=__args__.min_columns, min_row=__args__.min_rows)
    with open(__args__.context_path, 'r') as f_in:
        for object_id, line in enumerate(f_in):
            raw_entry = line.replace('\n', '').replace('\r', '')
            patterns.append(Pattern(instance=raw_entry, config=cfg, object=object_id+1))

    if __args__.context_type == 'o':
        to_be_tested = {}  # dictionary. key = user id, value = set of items to be tested
        with open(__args__.test, 'r') as f_true:
            for line in f_true:
                numbers = line.split('\t')
                user_id = int(numbers[0])
                if user_id in to_be_tested:
                    to_be_tested[user_id].add(int(numbers[1]) - 1)
                else:
                    to_be_tested[user_id] = {int(numbers[1]) - 1}

        file_out = open(__args__.context_path + '.prediction', 'w+')
        stored_extent = {}
        for target_user_id in to_be_tested:
            print (target_user_id)
            concepts = []
            empty_items = to_be_tested[target_user_id]
            for u in range(1, len(patterns)):
                if u == target_user_id:
                    continue
                intersection = patterns[target_user_id].intersect(patterns[u])
                if not stored_extent.get((str(u) + '-' + str(target_user_id))):
                    intersection.objects = {target_user_id, u}
                    stored_extent[str(target_user_id) + '-' + str(u)] = intersection.objects
                else:
                    intersection.objects = list(stored_extent.pop(str(u) + '-' + str(target_user_id)))
                concepts.append((intersection.partition_size_sum(), -1 * len(intersection.partition), intersection))

            concepts.sort(reverse=True)
            prediction_dict = {}
            for i in empty_items:
                prediction_dict[i] = [1000, -1000]
            for c in concepts:
                print ('concept', c[0], c[1], ':', c[2])
                filled_items = set()
                for similar_user in c[2].objects:
                    if similar_user == target_user_id:
                        continue
                    to_be_predicted = patterns[similar_user].set_of_filled_columns() & empty_items  # bisa disederhanakan
                    for tb in to_be_predicted:
                        print ('from u', similar_user, 'item', (tb + 1), patterns[similar_user].number_list[tb])
                        # file_out.write(str(target_user_id) + '\t' + str(tb + 1) + '\t' + str(patterns[similar_user].intervals[tb][0]) + '\t' + str(patterns[similar_user].intervals[tb][1]) + '\n')
                        prediction_dict[tb][0] = min(prediction_dict[tb][0], patterns[similar_user].number_list[tb])
                        prediction_dict[tb][1] = max(prediction_dict[tb][1], patterns[similar_user].number_list[tb])
                        filled_items.add(tb)
                #for fi in filled_items:
                    #file_out.write(str(target_user_id) + '\t' + str(fi + 1) + '\t' + str(prediction_dict[fi][0]) + '\t' + str(prediction_dict[fi][1]) + '\t' + str(c[0]) + '\t' + str(c[1]) + '\n')
                empty_items = empty_items - filled_items
        file_out.close()
    end_time = time()
    print(end_time - start_time)
