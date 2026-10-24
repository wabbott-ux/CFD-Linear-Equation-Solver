### Finite Volume Solver
### Will Abbott, 1013944706, Assignment 2

### This code takes a specified csv file and converts it into a sparse matrix A and corresponding vector b.
### The sparse matrix A is then solved based on the specified vector b to obtain a solution vector x.
### This code assumes csv file is to be given in a specific order:
### First row of csv file: values correspond to row indicies of non-zero entries of sparse matrix
### Second row of csv file: values correspond to column indicies of non-zero entries of sparse matrix
### Third row of csv file: values correspond to numerical value of non-zero entries of sparse matrix
### Fourth row of csv file: values correspond to numerical value of vector entries


import numpy as np
from scipy.sparse import coo_array
from scipy.sparse import csr_array
from scipy.sparse.linalg import bicg
from scipy.sparse.linalg import norm
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_array
from numpy.random import rand
from random import randint
import time

print("Program Start")
t_start = time.time()

def create_matrix(Matrix_A, Matrix_b):
    
    A_data = []
    A_row = []
    A_col = []
    b_data = []
    b_row = [] 


    with open(Matrix_A) as my_file:
        for line in my_file:
            Mat = line.strip("\n").split(",") #reading line of data from csv file
            A_data.append(float(Mat[0]))
            A_row.append(int(Mat[1]))
            A_col.append(int(Mat[2]))
            
    with open(Matrix_b) as my_file:
        for line in my_file:
            Mat = line.strip("\n").split(",") #reading line of data from csv file
            b_data.append(float(Mat[0]))
            b_row.append(int(Mat[1]))

    b_coords = (b_row,np.zeros(len(b_row),dtype="int"))
    A = coo_array((np.array(A_data),(np.array(A_row),np.array(A_col)))) #Create matrix A of size NxN
    b = coo_array((np.array(b_data),b_coords)) #Create vector b of size Nx1

    return [A,b,A_data]

finest = create_matrix("Matrix A Finest (320).csv","Matrix B Finest (320).csv")
A_finest = csr_array(finest[0])
b_finest = csr_array(finest[1])

fine = create_matrix("Matrix A Fine (160).csv","Matrix B Fine (160).csv")

A_fine = csr_array(fine[0])
b_fine = csr_array(fine[1])

coarse = create_matrix("Matrix A Coarse (80).csv","Matrix B Coarse (80).csv")
A_coarse = csr_array(coarse[0])
b_coarse = csr_array(coarse[1])



t_solve = time.time()
print(f"Setup finished in {round(t_solve-t_start,2)} seconds")

print(f"Starting Matrix Calculation")

x_finest = spsolve(A_finest,b_finest)
x_fine = spsolve(A_fine,b_fine)
x_coarse = spsolve(A_coarse,b_coarse)





t_finish = time.time()

print(f"Finished in {round(t_finish-t_solve,2)} seconds")



