import numpy as np
from scipy.sparse import csc_array
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_array
from numpy.random import rand

A_data = []
A_row = []
A_col = []

with open("new_file.txt") as my_file:

    for line in my_file:
        M = line.split(",")
        for i in M:
            A_data.append(int(i))

        
        

##A = lil_array((1000000, 1000000))
##A.setdiag(rand(1000000))

##A = A.tocsr()
##b = rand(1000000)
##x = spsolve(A, b)

print("done")
