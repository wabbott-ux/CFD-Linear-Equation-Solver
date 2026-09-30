### Sparse Matrix Linear Equation Solver
### Will Abbott, 1013944706, Assignment 1 Question 1

### This code generates a csv file containing data used to construct sparse matrix A and vector B for the Linear Equation Solver.
### The matrix will have random non-zero entries across the main diagonal, as well as generated randomly throughout.

import numpy as np
from random import randint
import time

N = 10**6 #Order of matrix can be specified
M = 10**5 #Tuning factor for how many random non zero entries to introduce, and what values they can have.

t1 = time.time()                    

with open("Matrix A.csv", "w") as csvfile:

        for i in range(0,N):
                csvfile.write(str(i)+",") #Writes all row indicies from 0->N-1
                
        for i in range(0,M):
                csvfile.write(str(randint(1,N-1))+",") #Introduces non-zero entries at random row incdicies
                
        csvfile.write("\n")
     
        for i in range(0,N):
                csvfile.write(str(i)+",") #Writes all column indicies from 0->N-1
                
        for i in range(0,M):
                csvfile.write(str(randint(1,N-1))+",") #Introduces non-zero entries at random column incdicies
                
        csvfile.write("\n")
     

        for i in range(0,N+M):
                csvfile.write(str(randint(1,M))+",") #Writes random values of non-zero entries, for appropriate length (Order of matrix + randomness)
        csvfile.write("\n")

        for i in range(0,N):
                csvfile.write(str(randint(1,M))+",") #Writes random values of non-zero entries for vector b
        csvfile.write("\n")
                

t2 = time.time()

print(f"Finished in {round(t2-t1,2)} seconds")
