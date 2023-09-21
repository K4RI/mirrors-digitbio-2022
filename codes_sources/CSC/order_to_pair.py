import argparse


def convert(filepath, to_binary):
    file_out1 = open('converted_greater_than.txt', 'w+')
    file_out2 = open('converted_lesser_than.txt', 'w+')
    with open(filepath) as f:
        number_of_attr = 0
        for line in f:
            numbers = line.split()
            if number_of_attr == 0:
                number_of_attr = len(numbers)
            for attr1 in range(0, number_of_attr):
                for attr2 in range(attr1 + 1, number_of_attr):
                    if to_binary == 1:
                        if int(numbers[attr1]) > int(numbers[attr2]):
                            file_out1.write('1 ')
                            file_out2.write('0 ')
                        elif int(numbers[attr1]) < int(numbers[attr2]):
                            file_out1.write('0 ')
                            file_out2.write('1 ')
                        else:
                            file_out1.write('0 ')
                            file_out2.write('0 ')
                    else:
                        if int(numbers[attr1]) > int(numbers[attr2]):
                            file_out1.write(str(attr1) + '-' + str(attr2) + ' ')
                        elif int(numbers[attr1]) < int(numbers[attr2]):
                            file_out2.write(str(attr1) + '-' + str(attr2) + ' ')
                        else:
                            file_out1.write(str(attr1) + '-' + str(attr2) + ' ')
                            file_out2.write(str(attr1) + '-' + str(attr2) + ' ')
            file_out1.write('\n')
            file_out2.write('\n')
    file_out1.close()
    file_out2.close()
    return


def convert_for_7(filepath):  # beware of the equality
    file_out = open('binary.txt', 'w+')
    with open(filepath) as f:
        number_of_attr = 0
        for line in f:
            numbers = line.split()
            if number_of_attr == 0:
                number_of_attr = len(numbers)
            for attr1 in range(0, number_of_attr):
                for attr2 in range(attr1 + 1, number_of_attr):
                    if float(numbers[attr1]) > float(numbers[attr2]):
                        file_out.write('g ')
                    #elif int(numbers[attr1]) < int(numbers[attr2]):
                        #file_out.write('l ')
                    else:
                        file_out.write('l ')
            file_out.write('\n')
    file_out.close()
    return


if __name__ == '__main__':
    __parser__ = argparse.ArgumentParser(description='Convert a matrix (e.g. a sequence) a numeric a pair')
    __parser__.add_argument('context_path', metavar='context_path', type=str, help='path to the formal context')
    __parser__.add_argument('output_format', metavar='output_format', type=int, help='1 = yes = binary formal context.\n'
                                                                                     '0 = no = special formal context for ex1\n'
                                                                             '7 = output \'lesser\' or \'greater\', for ex_7')
    __args__ = __parser__.parse_args()
    if __args__.output_format == 7:
        convert_for_7(__args__.context_path)
    else:
        convert(__args__.context_path, __args__.output_format)