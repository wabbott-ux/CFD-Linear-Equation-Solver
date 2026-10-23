### Sparse Matrix Linear Equation Solver
### Will Abbott, 1013944706, Assignment 1 Question 1

### This code generates a csv file containing data used to construct sparse matrix A and vector B for the Linear Equation Solver.
### The matrix will have random non-zero entries across the main diagonal, as well as generated randomly throughout.

import numpy as np
import time

L_x = 2
L_y = 2
Nodes = L_x*L_y

delta_x = 1/L_x
delta_y = 1/L_y

gamma = 20
hf = 10
Phi_x0 = 10
Phi_xL = 100
Phi_ext = 300

a_p_top = 3*gamma*(delta_y/delta_x) + gamma*(delta_x/delta_y)- hf*delta_x

a_p_bot = gamma*(3*(delta_y/delta_x) + (delta_x/delta_y))

a_p = 2*gamma*((delta_y/delta_x)+(delta_x/delta_y))

a_ns = -gamma*(delta_x/delta_y)

a_ew = -gamma*(delta_y/delta_x)


b_top_left = 2*gamma*(delta_y/delta_x)*Phi_x0 - hf*delta_x*Phi_ext

b_top_right = 2*gamma*(delta_y/delta_x)*Phi_xL - hf*delta_x*Phi_ext

b_bot_right = 2*gamma*(delta_y/delta_x)*Phi_xL

b_bot_left = 2*gamma*(delta_y/delta_x)*Phi_x0



D = {
        1:[a_p_top,0,a_ew,a_ns,0,b_top_left],
        L_x:[a_p_top,0,0,a_ns,a_ew,b_top_right],
        L_x*(L_y-1)+1:[a_p_bot,a_ns,a_ew,0,0,b_bot_left],
        L_y-1*L_x-1:[a_p_bot,a_ns,0,0,a_ew,b_bot_right]
        
        }

t1 = time.time()                    

with open("Matrix A.csv", "w") as csvfile:
        for node in range(1,Nodes+1):
                if node == 1:
                        for i in range(0,D[node]):
                                if D[node] != 0:
                                        csvfile.write(str(D[node][i])+",")
                                        csvfile.write(str(i)+",")
                                        csvfile.write(str(i)+",")
                                        
                                
                                
                        
t2 = time.time()

print(f"Finished in {round(t2-t1,2)} seconds")
