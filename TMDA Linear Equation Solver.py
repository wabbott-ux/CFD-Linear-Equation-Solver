import numpy as np
from scipy.sparse import coo_array
from scipy.sparse.linalg import cgs
from scipy.sparse import lil_array
from numpy.random import rand
from random import randint


N = 10**3
M = 100


A_data = []
A_row = []
A_col = []

A_order = [A_data,A_row,A_col]

counter = 0



with open("Matrix A.csv", "w") as csvfile:

        for j in range(0,N):
                csvfile.write(str(randint(0,M))+",")
        csvfile.write("\n")
        for j in range(0,N):
                csvfile.write(str(randint(0,N-1))+",")
        csvfile.write("\n")
        for j in range(0,N):
                csvfile.write(str(randint(0,N-1))+",")
        csvfile.write("\n")
        
        
##for i in range(0,N):
##    A_row.append(int(i))
##    A_col.append(int(i))
##    counter += 1

with open("Matrix A.csv") as my_file:
    for line in my_file:
        Mat = line.strip("\n").split(",")
        for i in range(0,len(Mat)-1):
                A_order[counter].append(int(Mat[i]))
        counter += 1

    

A = coo_array((np.array(A_data),(np.array(A_row),np.array(A_col))),shape=(N,N))

A = A.tocsr()
b = rand(N)
x = cgs(A, b)


print("done")
        
