### Finite Volume Solver
### Will Abbott, 1013944706, Assignment 2

### This code takes csv files and converts them into a sparse matrix A and corresponding vector b.
### The csv files provided represent the coefficient matrices of a system of linear equations.
### These linear equations are solved to find solutions for the equation provided in Assignment 2.
### The sparse matrix A is solved based on the specified vector b to obtain a solution vector x.
### This code assumes csv file is to be given in a specific order:
### First column of csv file: values at a particular point of the matrix
### Second column of csv file: row index of matrix value is located
### Third column of csv file: column index of matrix value is located

# Importing necessary functions from modules
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
import math as math
import matplotlib.pyplot as plt

print("Program Start")
t_start = time.time()

#Approach to solve matrices was changed from last report.
#Now a function is created which takes a matrix as an input, and solves it with its corresponding vector.

def create_matrix(Matrix_A):
    
    #Initializing lists to organize data received from csv files into values, rows, and columns
    A_data = []
    A_row = []
    A_col = []
    b_data = []
    b_row = [] 


    with open(Matrix_A+".csv") as my_file:
        for line in my_file:
            Mat = line.strip("\n").split(",") #reading line of data from csv file
            A_data.append(float(Mat[0]))
            A_row.append(int(Mat[1]))
            A_col.append(int(Mat[2]))
            
    with open(Matrix_A+" Vector.csv") as my_file:
        for line in my_file:
            Mat = line.strip("\n").split(",") #reading line of data from csv file
            b_data.append(float(Mat[0]))
            b_row.append(int(Mat[1]))

    b_coords = (b_row,np.zeros(len(b_row),dtype="int"))
    A = coo_array((np.array(A_data),(np.array(A_row),np.array(A_col)))) #Create matrix A of size NxN
    b = coo_array((np.array(b_data),b_coords)) #Create vector b of size Nx1

    return [A,b,A_data,A_row,A_col,b_data,b_coords]

#Creating matrices from specified csv files for each mesh size

finest = create_matrix("Matrix A Finest (320)")
A_finest = csr_array(finest[0])
b_finest = csr_array(finest[1])

fine = create_matrix("Matrix A Fine (160)")
A_fine = csr_array(fine[0])
b_fine = csr_array(fine[1])

coarse = create_matrix("Matrix A Coarse (80)")
A_coarse = csr_array(coarse[0])
b_coarse = csr_array(coarse[1])

non_uni_coarse = create_matrix("Matrix A Non Uniform Coarse (80)")
A_non_uni_coarse = csr_array(non_uni_coarse[0])
b_non_uni_coarse = csr_array(non_uni_coarse[1])

non_uni_fine = create_matrix("Matrix A Non Uniform Fine (160)")
A_non_uni_fine = csr_array(non_uni_fine[0])
b_non_uni_fine = csr_array(non_uni_fine[1])

non_uni_finest = create_matrix("Matrix A Non Uniform Finest (320)")
A_non_uni_finest = csr_array(non_uni_finest[0])
b_non_uni_finest = csr_array(non_uni_finest[1])

t_solve = time.time()
print(f"Setup finished in {round(t_solve-t_start,2)} seconds")

print(f"Starting Matrix Calculation")

#Calculate solution vector x, representing temperature phi in Assignment 2

x_finest = spsolve(A_finest,b_finest)
x_fine = spsolve(A_fine,b_fine)
x_coarse = spsolve(A_coarse,b_coarse)
x_non_uni_coarse = spsolve(A_non_uni_coarse,b_non_uni_coarse)
x_non_uni_fine = spsolve(A_non_uni_fine,b_non_uni_fine)
x_non_uni_finest = spsolve(A_non_uni_finest,b_non_uni_finest)

#Interpolate solution vector phi based on phi for the finest mesh size

def interpolate(x):

  return np.interp(np.linspace(0,1,num=len(x_finest)),
                                  np.linspace(0,1,num=len(x)),
                                  x)
  
#Interpolation performed for all mesh sizes

x_coarse_interpolated = interpolate(x_coarse)
x_fine_interpolated = interpolate(x_fine)
x_non_uni_coarse_interpolated = interpolate(x_non_uni_coarse)
x_non_uni_fine_interpolated = interpolate(x_non_uni_fine)

#Error norm calculated as specified in Assignment 2

def error_norm(F,C):

    error = 0

    for i in range(0,len(F)):
        error += (F[i]-C[i])**2

    return(error/len(F))**0.5


#Error norm calculated for all mesh sizes

error_norm_coarse = error_norm(x_finest,x_coarse_interpolated)
error_norm_fine = error_norm(x_finest,x_fine_interpolated)
error_norm_non_uni_coarse = error_norm(x_non_uni_finest,x_non_uni_coarse_interpolated)
error_norm_non_uni_fine = error_norm(x_non_uni_finest,x_non_uni_fine_interpolated)

#Order of convergence calculated as specified in Assignment 2 for both uniform and non uniform meshes

convergence = math.log10(error_norm_coarse/error_norm_fine)/math.log10((1/80)/(1/160))
non_uni_convergence = math.log10(error_norm_non_uni_coarse/error_norm_non_uni_fine)/math.log10((1/80)/(1/160))

t_finish = time.time()

print(f"Finished in {round(t_finish-t_solve,2)} seconds")

print("For Uniform Spacing Mesh:")
print(f"Error Norm of Coarse Mesh: {error_norm_coarse}")
print(f"Error Norm of Fine Mesh: {error_norm_fine}")
print(f"Order of Convergence: {convergence}")

print("For Non Uniform Spacing Mesh:")
print(f"Error Norm of Coarse Mesh: {error_norm_non_uni_coarse}")
print(f"Error Norm of Fine Mesh: {error_norm_non_uni_fine}")
print(f"Order of Convergence: {non_uni_convergence}")

#Creating contour plot for solution vector based on x,y points and corresponding node

def make_plot(X):
    xi = np.linspace(1,int(len(X)**0.5),num=(int(len(X)**0.5)))
    yi = xi
    zi = np.reshape(X,(int(len(X)**0.5),int(len(X)**0.5)))

    cntr1 = plt.contourf(xi, yi, zi, levels=100, cmap="RdBu_r")

    plt.colorbar(cntr1)

    plt.show()

                    
make_plot(x_finest) #Creating countour plot of finest mesh

