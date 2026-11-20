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

N_x = 10
N_y = N_x

delta_x = 1/N_x
delta_y = 1/N_y

gamma = 5
phi_B = 100
phi_b = 0

u_x = 2
u_y = 2

D = gamma/delta_x

#Central Differencing
a_p_cd = 4*D
a_n_cd = -(D - 0.5*u_y)
a_e_cd = -(D - 0.5*u_x)
a_s_cd = -(D + 0.5*u_y)
a_w_cd = -(D + 0.5*u_x)

a_p_north_west_cd = 6*D + 0.5*u_x - 0.5*u_y
a_p_north_east_cd = 6*D - 0.5*u_x - 0.5*u_y
a_p_south_west_cd = 6*D + 0.5*u_x + 0.5*u_y
a_p_south_east_cd = 6*D - 0.5*u_x + 0.5*u_y

a_p_north_cd = 5*D - 0.5*u_y
a_p_west_cd = 5*D + 0.5*u_x
a_p_east_cd = 5*D - 0.5*u_x
a_p_south_cd = 5*D + 0.5*u_y

b_north_west_cd = (4*D - u_y + u_x)*phi_B
b_north_east_cd = (2*D - u_y)*phi_B + (2*D - u_x)*phi_b
b_south_west_cd = (2*D + u_y)*phi_B + (2*D + u_x)*phi_b
b_south_east_cd = (4*D + u_y - u_x)*phi_b

b_north_cd = (2*D - u_y)*phi_B
b_west_cd = (2*D + u_x)*phi_B
b_east_cd = (2*D - u_x)*phi_b
b_south_cd = (2*D + u_y)*phi_b
  
bc_list = {}
vector_list = {}

for row in range(1,N_y+1):
        for col in range(1,N_x+1):
                if row == 1 and col == 1: #Equation for top left corner node 
                        bc_list.update({(row,col):[a_p_north_west_cd,0,a_e_cd,a_s_cd,0]})
                        vector_list.update({(row,col):b_north_west_cd})
                        
                elif row == 1 and col == N_x: #Equation for top right corner node
                        bc_list.update({(row,col):[a_p_north_east_cd,0,0,a_s_cd,a_w_cd]})
                        vector_list.update({(row,col):b_north_east_cd})

                elif row == N_y and col == 1: #Equation for bottom left corner node
                        bc_list.update({(row,col):[a_p_south_west_cd,a_n_cd,a_e_cd,0,0]})
                        vector_list.update({(row,col):b_south_west_cd})

                elif row == N_y and col == N_x: #Equation for bottom right corner node
                       bc_list.update({(row,col):[a_p_south_east_cd,a_n_cd,0,0,a_w_cd]})
                       vector_list.update({(row,col):b_south_east_cd})

                elif row == 1: #Equation for any node along top boundary
                        bc_list.update({(row,col):[a_p_north_cd,0,a_e_cd,a_s_cd,a_w_cd]})
                        vector_list.update({(row,col):b_north_cd})

                elif col ==1: #Equation for any node along left boundary
                        bc_list.update({(row,col):[a_p_west_cd,a_n_cd,a_e_cd,a_s_cd,0]})
                        vector_list.update({(row,col):b_west_cd})

                elif col == N_x: #Equation for any node along right boundary
                        bc_list.update({(row,col):[a_p_east_cd ,a_n_cd,0,a_s_cd,a_w_cd]})
                        vector_list.update({(row,col):b_east_cd})

                elif row == N_y: #Equation for any node along bottom boundary
                        bc_list.update({(row,col):[a_p_south_cd,a_n_cd,a_e_cd,0,a_w_cd]})
                        
                else: #Equation for any node not governed by a boundary condition
                        bc_list.update({(row,col):[a_p_cd,a_n_cd,a_e_cd,a_s_cd,a_w_cd]})

A_data = []
A_row = []
A_col = []
b_data = []
b_row = [] 

for row in range(1,N_y+1): 
        for col in range(1,N_x+1):

                node = (row-1)*N_x + col -1 #Converts x,y coordinate of point to node value

                if bc_list[(row,col)][0] != 0: #Writes 'ap' coefficient value for current principal node for each x,y
                        
                        A_data.append(bc_list[(row,col)][0])
                        A_row.append(node)
                        A_col.append(node)

                if bc_list[(row,col)][1] != 0: #Writes 'an' coefficient value north of principal node

                        A_data.append(bc_list[(row,col)][1])
                        A_row.append(node)
                        A_col.append((row-2)*N_x + col-1)

                if bc_list[(row,col)][2] != 0: #Writes 'ae' coefficient value east of principal node
                        
                        A_data.append(bc_list[(row,col)][2])
                        A_row.append(node)
                        A_col.append((row-1)*N_x + col)

                if bc_list[(row,col)][3] != 0: #Writes 'as' coefficient value south of principal node

                        A_data.append(bc_list[(row,col)][3])
                        A_row.append(node)
                        A_col.append((row)*N_x + col-1)

                if bc_list[(row,col)][4] != 0: #Writes 'aw' coefficient value west of principal node

                        A_data.append(bc_list[(row,col)][4])
                        A_row.append(node)
                        A_col.append((row-1)*N_x + col-2)

for row in range(1,N_y+1):
        for col in range(1,N_x+1):
                
                if (row,col) in vector_list: #Writes constant for vector b for node x,y

                        b_data.append(vector_list[(row,col)])
                        b_row.append((row-1)*N_x + col-1)
                        
b_coords = (b_row,np.zeros(len(b_row),dtype="int"))
A = coo_array((np.array(A_data),(np.array(A_row),np.array(A_col)))) #Create matrix A of size NxN
b = coo_array((np.array(b_data),b_coords)) #Create vector b of size Nx1

#Calculate solution vector x, representing temperature phi
x = spsolve(csr_array(A),csr_array(b))

t_finish = time.time()

print(f"Finished in {round(t_finish-t_start,2)} seconds")

#Creating contour plot for solution vector based on x,y points and corresponding node

xi = np.linspace(1,int(len(x)**0.5),num=(int(len(x)**0.5)))
yi = xi[::-1]
zi = np.reshape(x,(int(len(x)**0.5),int(len(x)**0.5)))

cntr1 = plt.contourf(xi, yi, zi, levels=100, cmap="RdBu_r")

plt.colorbar(cntr1)

plt.show()



