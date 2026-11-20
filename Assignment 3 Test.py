### Finite Volume Solver
### Will Abbott, 1013944706, Assignment 3

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

#L_x = 1
#L_y = 1

N_x = 100

delta_x = 1/N_x


gamma = 5
phi_B = 100
phi_b = 0

u_x = 2


D = gamma/delta_x

#Central Differencing
a_p_cd = 2*D
a_e_cd = -(D - 0.5*u_x)
a_w_cd = -(D + 0.5*u_x)


a_p_west_cd = 3*D + 0.5*u_x
a_p_east_cd = 3*D - 0.5*u_x

b_west_cd = (2*D+u_x)*phi_B
b_east_cd = (2*D-u_x)*phi_b

bc_list = {}
vector_list = {}

for col in range(1,N_x+1):
        if col == 1: #Equation for top left corner node 
                bc_list.update({col:[a_p_west_cd,0,a_e_cd]})
                vector_list.update({col:b_west_cd})

        elif col == N_x: #Equation for any node along right boundary
                bc_list.update({col:[a_p_east_cd,a_w_cd,0]})
                vector_list.update({col:b_east_cd})
                
        else: #Equation for any node not governed by a boundary condition
                bc_list.update({col:[a_p_cd,a_w_cd,a_e_cd]})

A_data = []
A_row = []
A_col = []
b_data = []
b_row = [] 

for col in range(1,N_x+1):

        if bc_list[col][0] != 0: #Writes 'ap' coefficient value for current principal node for each x,y
                
                A_data.append(bc_list[col][0])
                A_row.append(col-1)
                A_col.append(col-1)

        if bc_list[col][1] != 0: #Writes 'aw' coefficient value west of principal node

                A_data.append(bc_list[col][1])
                A_row.append(col-1)
                A_col.append(col-2)

        if bc_list[col][2] != 0: #Writes 'ae' coefficient value east of principal node
                
                A_data.append(bc_list[col][2])
                A_row.append(col-1)
                A_col.append(col)

for col in range(1,N_x+1):
        
        if col in vector_list: #Writes constant for vector b for node x,y

                b_data.append(vector_list[col])
                b_row.append(col-1)
                        
b_coords = (b_row,np.zeros(len(b_row),dtype="int"))
A = coo_array((np.array(A_data),(np.array(A_row),np.array(A_col)))) #Create matrix A of size NxN
b = coo_array((np.array(b_data),b_coords)) #Create vector b of size Nx1

#Calculate solution vector x, representing temperature phi
x = spsolve(csr_array(A),csr_array(b))

t_finish = time.time()

print(f"Finished in {round(t_finish-t_start,2)} seconds")

#Creating contour plot for solution vector based on x,y points and corresponding node

xi = np.linspace(1,len(x),num=(len(x)))

cntr1 = plt.bar(xi,x)

#plt.colorbar(cntr1)

plt.show()


