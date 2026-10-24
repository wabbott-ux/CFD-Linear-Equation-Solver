### Sparse Matrix Linear Equation Solver
### Will Abbott, 1013944706, Assignment 1 Question 1

### This code generates a csv file containing data used to construct sparse matrix A and vector B for the Linear Equation Solver.
### The matrix will have random non-zero entries across the main diagonal, as well as generated randomly throughout.

import numpy as np
import time

L_x = 80
L_y = L_x
Nodes = L_x*L_y

delta_x = 1/L_x
delta_y = 1/L_y

gamma = 20
hf = 10
Phi_x0 = 10
Phi_xL = 100
Phi_ext = 300

a_p_top_corner = gamma*(3*(delta_y/delta_x) +(delta_x/delta_y))- hf*delta_x

a_p_top = gamma*(2*(delta_y/delta_x) +(delta_x/delta_y))- hf*delta_x

a_p_ew = gamma*(3*(delta_y/delta_x) + 2*(delta_x/delta_y))

a_p_bot_corner = gamma*(3*(delta_y/delta_x) + (delta_x/delta_y))

a_p_bot = gamma*(2*(delta_y/delta_x) + (delta_x/delta_y))

a_p = 2*gamma*((delta_y/delta_x)+(delta_x/delta_y))

a_ns = -gamma*(delta_x/delta_y)

a_ew = -gamma*(delta_y/delta_x)


b_top_left = 2*gamma*(delta_y/delta_x)*Phi_x0 - hf*delta_x*Phi_ext

b_top_right = 2*gamma*(delta_y/delta_x)*Phi_xL - hf*delta_x*Phi_ext

b_top = -hf*delta_x*Phi_ext

b_right = 2*gamma*(delta_y/delta_x)*Phi_xL

b_left = 2*gamma*(delta_y/delta_x)*Phi_x0

corners = [(1,1),(1,L_x),(L_y,1),(L_y,L_x)]


bc_list = {}
vector_list = {}

for row in range(1,L_y+1):
        for col in range(1,L_x+1):
                if row == 1 and col == 1:
                        bc_list.update({(row,col):[a_p_top_corner,0,a_ew,a_ns,0]})
                        vector_list.update({(row,col):b_top_left})
                        
                elif row == 1 and col == L_x:
                        bc_list.update({(row,col):[a_p_top_corner,0,0,a_ns,a_ew]})
                        vector_list.update({(row,col):b_top_right})

                elif row == L_y and col == 1:
                        bc_list.update({(row,col):[a_p_bot_corner,a_ns,a_ew,0,0]})
                        vector_list.update({(row,col):b_left})

                elif row == L_y and col == L_x:
                       bc_list.update({(row,col):[a_p_bot_corner,a_ns,0,0,a_ew]})
                       vector_list.update({(row,col):b_right})

                elif row == 1:
                        bc_list.update({(row,col):[a_p_top,0,a_ew,a_ns,a_ew]})
                        vector_list.update({(row,col):b_top})

                elif col ==1:
                        bc_list.update({(row,col):[a_p_ew,a_ns,a_ew,a_ns,0]})
                        vector_list.update({(row,col):b_left})

                elif col == L_x:
                        bc_list.update({(row,col):[a_p_ew,a_ns,0,a_ns,a_ew]})
                        vector_list.update({(row,col):b_right})

                elif row == L_y:
                        bc_list.update({(row,col):[a_p_bot,a_ns,a_ew,0,a_ew]})
                else:
                        bc_list.update({(row,col):[a_p,a_ns,a_ew,a_ns,a_ew]})
                        
        

t1 = time.time()                    

with open("Matrix A.csv", "w") as csvfile:

        
        for row in range(1,L_y+1):
                for col in range(1,L_x+1):

                        node = (row-1)*L_x + col -1

                        if bc_list[(row,col)][0] != 0:

                                csvfile.write(str(bc_list[(row,col)][0])+",")
                                csvfile.write(str(node)+",")
                                csvfile.write(str(node)+",")
                                csvfile.write("\n")

                        if bc_list[(row,col)][1] != 0:

                                csvfile.write(str(bc_list[(row,col)][1])+",")
                                csvfile.write(str(node)+",")
                                csvfile.write(str((row-2)*L_x + col-1)+",")
                                csvfile.write("\n")

                        if bc_list[(row,col)][2] != 0:
                                
                                csvfile.write(str(bc_list[(row,col)][2])+",")
                                csvfile.write(str(node)+",")
                                csvfile.write(str((row-1)*L_x + col)+",")
                                csvfile.write("\n")

                        if bc_list[(row,col)][3] != 0:

                                csvfile.write(str(bc_list[(row,col)][3])+",")
                                csvfile.write(str(node)+",")
                                csvfile.write(str((row)*L_x + col-1)+",")
                                csvfile.write("\n")

                        if bc_list[(row,col)][4] != 0:

                                csvfile.write(str(bc_list[(row,col)][4])+",")
                                csvfile.write(str(node)+",")
                                csvfile.write(str((row-1)*L_x + col-2)+",")
                                csvfile.write("\n")

with open("Matrix B.csv", "w") as csvfile:

        for row in range(1,L_y+1):
                for col in range(1,L_x+1):
                        
                        if (row,col) in vector_list:
                                
                                csvfile.write(str(vector_list[(row,col)])+",")
                                csvfile.write(str((row-1)*L_x + col-1)+",")
                                csvfile.write("\n")
                       

                                
                                
                        
t2 = time.time()

print(f"Finished in {round(t2-t1,2)} seconds")
