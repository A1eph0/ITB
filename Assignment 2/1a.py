import csv


# takes two sequences and calculates the distace between them
def distance_calculation(seq1, seq2):
    # input seq1 and seq2 are strings

    if seq1 == seq2:
        return float(0)
    mismatch_count = 0
    length = len(seq1)
    blank_length = 0            # counts number of common blanks ('-')
    for i in range(length):
        if seq1[i] == seq2[i] and seq1[i] == "-":
            blank_length += 1
        if seq1[i] != seq2[i]:
            mismatch_count += 1
    true_length = length - blank_length    # actual length to be considered
    distance = float(float(mismatch_count)/float(true_length))
    return distance


# opening files and defining writer
in_file = open("Nucleotide_alignment.txt", "r")
out_file = open("Ndistance.txt", "w")
writer = csv.writer(out_file, delimiter='|', quotechar='"',
                    quoting=csv.QUOTE_NONNUMERIC)

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
        matrix[i][j] = distance_calculation(name_to_seq[list_of_names[i]],
                                            name_to_seq[list_of_names[j]])
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
out_file.close()

# end of code
