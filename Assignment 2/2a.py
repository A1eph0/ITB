import csv
import numpy as np


# takes two sequences and calculates the distace between them
def distance_calculation(seq1, seq2, index, mat):
    # input seq1 and seq2 are strings, index is map, matrix is a matrix

    distance = 0
    length = len(seq1)
    is_open = False
    for i in range(length):
        a = seq1[i]
        b = seq2[i]
        if a == "-" and b == "-":
            continue                # goes directly to next iteration of loop
        if a == "-" or b == "-":
            if is_open:
                distance -= 1
            else:
                is_open = True
                distance -= 11
        else:
            is_open = False
            distance += int(mat[index[a]][index[b]])

    return distance


# calculating the redunant length due to overlaying gaps ("-" v/s "-")
def extra_len_calc(seq1, seq2):
    red_distance = 0
    length = len(seq1)
    for i in range(length):
        if seq1 == "-" and seq2 == "-":
            red_distance += 1
    return red_distance


def final_distance_calculation(seq1, seq2, index, mat):

    """ distance calculation using:
    Sonnhammer, E. L., & Hollich, V. (2005)
    [BMC Bioinformatics, 6(1), 108. doi:10.1186/1471-2105-6-108]
    """

    sigma0 = (-4)
    actual_length = len(seq1) - extra_len_calc(seq1, seq2)
    sigma_r = sigma0*actual_length

    sigma_1_2 = distance_calculation(seq1, seq2, index, mat)
    sigma_N = sigma_1_2 - (sigma_r * actual_length)

    sigma_1_1 = distance_calculation(seq1, seq1, index, mat)
    sigma_2_2 = distance_calculation(seq2, seq2, index, mat)
    sigma_U = float((float(sigma_1_1) + float(sigma_2_2))/float(2))
    sigma_UN = sigma_U - (sigma_r * actual_length)

    distance = (np.log(sigma_UN) - np.log(sigma_N))*100

    return distance


# opening files and defining reader/writer
in_file = open("Protein_alignment.txt", "r")
use_file = open("BLOSUM62.txt", "r")
reader = csv.reader(use_file, delimiter=" ")
out_file = open("Pdistance.txt", "w")
writer = csv.writer(out_file, delimiter='|', quotechar='"',
                    quoting=csv.QUOTE_NONNUMERIC)

use_matrix = list(reader)       # regenerating matrix for score calculation
use_matrix.pop(0)               # removing row with headers
amino_to_index = {}             # maps amino acid to index

# making  matrix without header in collumns and mapping indices to amino acid
for i in range(len(use_matrix)):
    amino_to_index[use_matrix[i][0]] = i
    use_matrix[i].pop(0)

# removing extra spaces
for i in range(len(use_matrix)):
    temp = []
    for j in range(len(use_matrix[i])):
        if use_matrix[i][j] == "":
            temp.insert(0, j)
    for j in temp:
        use_matrix[i].pop(j)
list_of_names = []  # list of all names for sequences
name_to_seq = {}    # dictionary mapping name to corresponding sequence
num = 0             # counter for storing number of sequences obtained as input

temp_seq = ''
temp_name = ''

while True:
    line = in_file.readline()                   # reading file line-by-line

    if not line:                                # break condition for EOF
        name_to_seq[temp_name] = temp_seq
        num += 1
        break

    if line[0] == '>':
        if temp_seq != '':
            name_to_seq[temp_name] = temp_seq
            num += 1
            temp_seq = ''
        temp_name = line[1:]
        temp_name = temp_name[:-1]
        list_of_names.append(temp_name)
    else:
        if line[-1] == "\n":
            line = line[:-1]
        temp_seq += line


# initiating a square matrix of dimension num*num
matrix = [list(range(1 + num * i, 1 + num * (i + 1)))
          for i in range(num)]

# filling in distance values
for i in range(num):
    for j in range(i, num):
        matrix[i][j] = final_distance_calculation(name_to_seq[list_of_names
                                                              [i]],
                                                  name_to_seq[list_of_names
                                                              [j]],
                                                  amino_to_index,
                                                  use_matrix)
        matrix[j][i] = matrix[i][j]

# adding header row
list_of_names.insert(0, "*")
matrix.insert(0, list_of_names)

# adding header collumn and giving output
writer.writerow(matrix[0])
for i in range(1, num+1):
    matrix[i].insert(0, list_of_names[i])
    writer.writerow(matrix[i])

# closing files
in_file.close()
use_file.close()
out_file.close()

# end of code
