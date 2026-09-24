import numpy as np
from scipy.sparse import csc_array
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_array

def str_to_float(string):

    data = string.split(',')
    
    for i in range(0,len(data)):
        data[i] = float(data[i])
        
    return np.array(data)

data_string = input("Insert Data (Up to Down, then Left to Right) Separated by Commas e.g. '1,2,3': ")

row_ind_string = input("Insert Row Indices Separated by Commas e.g. '1,2,3': ")

col_ptr_string = input("Insert Row Indices Separated by Commas e.g. '1,2,3': ")
    

A = csc_array([[3, 2, 0], [1, -1, 0], [0, 5, 1]])
B = csc_array([[2, 0], [-1, 0], [2, 0]])
x = spsolve(A, B)

print(x)

        
