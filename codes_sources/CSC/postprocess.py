import sys


if __name__ == "__main__":
    path = sys.argv[1]
    string_out = ''

    # caesarian
    caesarian_class = list()
    with open('caesarian-class.txt', 'r') as f:
        for line in f:
            caesarian_class.append(line.replace('\n', '').replace('\r', ''))

    file_out = open('caesarian-postproc.txt', 'w+')
    with open(path, 'r') as f:
        for line in f:
            extent_intent = line.replace('\n', '').replace('\r', '').split(':')
            extent = extent_intent[2].replace('[', '').replace(']', '').split(',')
            for extent_member in extent:
                file_out.write(caesarian_class[int(extent_member)])
            file_out.write('\n')

    # dishonest
    '''dishonest_class = list()
    with open('dishonest-class.txt', 'r') as f:
        for line in f:
            if line.replace('\n', '').replace('\r', '') == 'trustworthy':
                dishonest_class.append('0')
            else:
                dishonest_class.append('1')

    file_out = open('dishonest-postproc.txt', 'w+')
    with open(path, 'r') as f:
        for line in f:
            extent_intent = line.replace('\n', '').replace('\r', '').split(':')
            extent = extent_intent[2].replace('[', '').replace(']', '').split(',')
            for extent_member in extent:
                file_out.write(dishonest_class[int(extent_member)])
            file_out.write('\n')'''

