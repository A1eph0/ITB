import csv

# opening files and defining reader
in_file = open("Ndistance.txt", "r")
out_file = open("Final_N_out.txt", "w")
reader = csv.reader(in_file, delimiter="|")


matrix = list(reader)       # regenerating matrix

matrix.pop(0)               # removing row with headers

list_of_names = []          # list of names
name_to_row = {}            # maps the names to the corresponding row
distance_from_zero = {}     # maps names to corresponding hieght
num = len(matrix)           # stores number of sequences

# making lower triangular matrix using the original matrix
for i in range(num):
    list_of_names.append("\"" + matrix[i][0] + "\"")
    temp_list = matrix[i][1:i+2]
    name_to_row[list_of_names[i]] = temp_list
    distance_from_zero[list_of_names[i]] = float(0)

# generating Newick Tree using UPGMA method
while len(list_of_names) != 1:
    min_val = float(10000)      # initiated with large value (inf)
    epsilon = 1e-5              # to account for floating point errors

    # to store the indices of sequences to be clustered
    min_i = int()
    min_j = int()

    # finding minimum distance value in matrix
    for i in range(len(list_of_names)):
        for j in range(i):
            if min_val-float(name_to_row[list_of_names[i]][j]) > (epsilon):
                min_val = float(name_to_row[list_of_names[i]][j])
                min_i = i
                min_j = j

    # modifying the row corresponding to lower index
    for i in range(len(name_to_row[list_of_names[min_j]])):
        name_to_row[list_of_names[min_j]][i] = str((float(name_to_row
                                                          [list_of_names
                                                           [min_j]][i]) +
                                                    float(name_to_row
                                                          [list_of_names
                                                           [min_i]][i]))
                                                   / float(2))

    # modifying the row corresponding to values between lower and higher index
    for i in range(min_j+1, min_i):
        name_to_row[list_of_names[i]][min_j] = str((float(name_to_row
                                                          [list_of_names
                                                           [i]][min_j]) +
                                                    float(name_to_row
                                                          [list_of_names
                                                           [min_i]][i]))
                                                   / float(2))

    # modifying the row corresponding to values above the higher index
    for i in range(min_i+1, len(list_of_names)):
        name_to_row[list_of_names[i]][min_j] = str((float(name_to_row
                                                          [list_of_names
                                                           [i]][min_j]) +
                                                    float(name_to_row
                                                          [list_of_names
                                                           [i]][min_i]))
                                                   / float(2))
        name_to_row[list_of_names[i]].pop(min_i)

    # clustering with the format ("seq1": distance1, "seq2" : distance2)
    temp_var = ("(" + list_of_names[min_j] +
                " : " + str(float(min_val)/float(2) -
                            float(distance_from_zero[list_of_names[min_j]])) +
                ", " + list_of_names[min_i] +
                " : " + str(float(min_val)/float(2)
                            - float(distance_from_zero[list_of_names[min_i]]))
                + ")")

    temp_dist = float(min_val)/float(2)

    # deleting row of higher index and creating cluster
    del name_to_row[list_of_names[min_i]]
    list_of_names.pop(min_i)
    name_to_row[temp_var] = name_to_row[list_of_names[min_j]]
    distance_from_zero[temp_var] = temp_dist
    del name_to_row[list_of_names[min_j]]
    list_of_names[min_j] = temp_var

# giving output
out_file.write(list_of_names[0])

# closing files
in_file.close()
out_file.close()

# end of code
