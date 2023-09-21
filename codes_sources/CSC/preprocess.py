import sys


if __name__ == "__main__":
    path = sys.argv[1]
    string_out = ''

    # BicatYeast.txt
    file_out = open('BicatYeast-binary.txt', 'w+')
    with open(path, 'r') as f:
        for line in f:
            numbers = line.replace('\n', '').replace('\r', '').split()
            for number in numbers:
                if float(number) < 0:
                    string_out += '- '
                elif float(number) > 0:
                    string_out += '+ '
            file_out.write(string_out[:-1] + '\n')
            string_out = ''

    # bicpam - disjoint equivalence blocks
    '''file_out = open('dataConstant60x500-discr.txt', 'w+')
    max = 0
    min = 0
    min2 = 0
    max2 = 0
    with open(path, 'r') as f:
        for line in f:
            numbers = line.replace('\n', '').replace('\r', '').split()
            for number in numbers:
                if int(number) > max:
                    max = int(number)
                if int(number) < min:
                    min = int(number)
                converted = int(number) / 10  # <==  penting, menentukan banyaknya distinct values
                if converted < min2:
                    min2 = converted
                if converted > max2:
                    max2 = converted
                string_out += str(converted) + ' '
            file_out.write(string_out[:-1] + '\n')
            string_out = ''
    print 'before min = ', min, ' max = ', max
    print 'after  min = ', min2, ' max = ', max2'''

    # lymphoma.csv
    '''file_out = open('binary.txt', 'w+')
    with open(path, 'r') as f:
        for line in f:
            numbers = line.replace('\n', '').replace('\r', '').split(',')
            for number in numbers:
                #if int(number) == 999:
                    #string_out += '999 '
                if int(number)<0:
                    string_out += '- '
                else:
                    string_out += '+ '
            file_out.write(string_out[:-1] + '\n')
            string_out = ''
    '''

    # breast-cancer-wisconsin
    '''file_out = open('breast-cancer-binary.txt', 'w+')
    with open(path, 'r') as f:
        for line in f:
            numbers = line.replace('\n', '').replace('\r', '').split(',')
            for i in range(len(numbers)):
                if i == 0 or i == 10:
                    continue
                elif numbers[i] == '?':
                    string_out += '- '
                elif int(numbers[i]) < 6:
                    string_out += '- '
                else:
                    string_out += '+ '
            file_out.write(string_out[:-1] + '\n')
            string_out = ''
            '''

    # contraceptive
    '''file_out = open('contraceptive-binary.txt', 'w+')
    with open(path, 'r') as f:
        for line in f:
            numbers = line.replace('\n', '').replace('\r', '').split(',')
            for i in range(len(numbers)):
                if i == 0:
                    if int(numbers[i]) <= 30:
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 1 or i == 2 or i == 7:
                    if numbers[i] == '1':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 3:
                    if int(numbers[i]) <= 2:
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 4 or i == 5 or i == 8:
                    if numbers[i] == '0':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 6:
                    if int(numbers[i]) <= 2:
                        string_out += '- '
                    else:
                        string_out += '+ '
            file_out.write(string_out[:-1] + '\n')
            string_out = ''
            '''

    #mirTarBase
    '''file_out = open('mirTarBase-binary.txt', 'w+')
    mrna_dict = {}
    with open('mirTarBase_mRNA.txt', 'r') as f:
        line_number = 0
        for line in f:
            mrna_dict[line.replace('\n', '').replace('\r', '')] = line_number
            line_number += 1

    before = ''
    arr = ['-'] * 1730
    with open(path, 'r') as f:
        for line in f:
            cells = line.split()
            if cells[0] != before:
                arr_out = ' '.join(arr)
                file_out.write(arr_out + '\n')
                arr = ['-'] * 1730
            else:
                arr[mrna_dict[cells[1]]] = '+'
            before = cells[0]
        arr_out = ' '.join(arr)
        file_out.write(arr_out + '\n')
        arr = ['-'] * 1730'''

    # caesarian
    '''file_out = open('caesarian-binary.txt', 'w+')
    with open(path, 'r') as f:
        for line in f:
            numbers = line.replace('\n', '').replace('\r', '').split(',')
            for i in range(len(numbers)):
                if i == 0:  # age
                    if int(numbers[i]) <= 30:
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 1:  # delivery number
                    if int(numbers[i]) < 3:
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 2:  # delivery time
                    if numbers[i] == '1':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 3:  # blood pressure
                    if numbers[i] == '1':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 4:  # heart problem
                    if numbers[i] == '0':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 5:  # class
                    if numbers[i] == '0':
                        string_out += '- '
                    else:
                        string_out += '+ '
            file_out.write(string_out[:-1] + '\n')
            string_out = ''
'''

    # tripadvisor_review
    '''file_out = open('tripadvisor-binary.txt', 'w+')
    with open(path, 'r') as f:
        for line in f:
            numbers = line.replace('\n', '').replace('\r', '').split(',')
            for number in numbers:
                if float(number) <= 2:
                    string_out += '- '
                else:
                    string_out += '+ '
            file_out.write(string_out[:-1] + '\n')
            string_out = ''
            '''

    # student-performance-India
    '''file_out = open('student-binary.txt', 'w+')
    with open(path, 'r') as f:
        for line in f:
            entries = line.replace('\n', '').replace('\r', '').split(',')
            for i in range(len(entries)):
                if i == 0:  # gender
                    if entries[i] == 'M':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 1:  # caste
                    if entries[i] == 'G':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif 2 <= i <= 5:  # percentage
                    if entries[i] == 'Fail':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 6:
                    if entries[i] == 'Y':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 7:
                    if entries[i] == 'Married':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 8:
                    if entries[i] == 'T':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 9:
                    if entries[i] == 'Free':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 10:
                    if entries[i] == 'Medium' or entries[i] == 'Low':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 11 or i == 16 or i == 20:
                    if entries[i] == 'Small':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 12 or i == 13:
                    if entries[i] == 'Il' or entries[i] == 'Um' or entries[i] == '10' or entries[i] == '12':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 14 or i == 15:
                    if entries[i] == 'Retired':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 17 or i == 21:
                    if entries[i] == 'Poor':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 18:
                    if entries[i] == 'Govt':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 19:
                    if entries[i] == 'Eng':
                        string_out += '- '
                    else:
                        string_out += '+ '
            file_out.write(string_out[:-1] + '\n')
            string_out = ''
            '''

    # dishonesty
    '''file_out = open('dishonest-binary.txt', 'w+')
    with open(path, 'r') as f:
        for line in f:
            numbers = line.replace('\n', '').replace('\r', '').split()
            for i in range(len(numbers)):
                if i == 0:
                    if numbers[i] == 'CT_range_1' or numbers[i] == 'CT_range_2':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 1:
                    if numbers[i] == 'CU_range_1' or numbers[i] == 'CU_range_2':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 2:
                    if numbers[i] == 'LT_range_1' or numbers[i] == 'LT_range_2':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 3:
                    if numbers[i] == 'ECommerce':
                        string_out += '- '
                    else:
                        string_out += '+ '
                elif i == 4:
                    if numbers[i] == 'trustworthy':
                        string_out += '- '
                    else:
                        string_out += '+ '
            file_out.write(string_out[:-1] + '\n')
            string_out = ''
            '''
