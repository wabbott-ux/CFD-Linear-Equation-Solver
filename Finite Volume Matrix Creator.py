### Finite Volume Method Coefficient Matrix Creator
### Will Abbott, 1013944706, Assignment 2

### This code generates a csv file containing data used to construct
### the coefficient matrices and solution vectors for the Finite Volume Method
### equation that was provided in Assignment 2.

# Importing necessary functions and modules
import numpy as np
import time

t1 = time.time()


#Approach to solve matrices was changed from last report.
#Now a function is created which takes a matrix name and size as an input.
#Function creates a csv file using equations obtained from using Finite Volume Method (FVM) from Assignment 2.
#Each row of the created csv file corresponds to an ap, an, as, or aw value from the FVM equations.
#Each value also has a corresponding row and column value which will determine it's position in the matrix.
#CSV file will be turned into a coefficient matrix and then solved.

def Matrix_A_Creator (Matrix_Name,L_x,L_y):

        #Initializing constants for FVM equations

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

        #Initializing lists for specific boundary condition equations

        bc_list = {} #Each node x,y is stored in a dictionary with associated coefficients from FVM equation (For matrix A)
        vector_list = {} #Each node x,y is stored in a dictionary with associated constant values from FVM equation (For vector b)

        for row in range(1,L_y+1): #Cycles through each point in the domain, with coordinate x = col and coordinate y = row
                for col in range(1,L_x+1):
                        
                        if row == 1 and col == 1: #Equation for top left corner node 
                                bc_list.update({(row,col):[a_p_top_corner,0,a_ew,a_ns,0]})
                                vector_list.update({(row,col):b_top_left})
                                
                        elif row == 1 and col == L_x: #Equation for top right corner node
                                bc_list.update({(row,col):[a_p_top_corner,0,0,a_ns,a_ew]})
                                vector_list.update({(row,col):b_top_right})

                        elif row == L_y and col == 1: #Equation for bottom left corner node
                                bc_list.update({(row,col):[a_p_bot_corner,a_ns,a_ew,0,0]})
                                vector_list.update({(row,col):b_left})

                        elif row == L_y and col == L_x: #Equation for bottom right corner node
                               bc_list.update({(row,col):[a_p_bot_corner,a_ns,0,0,a_ew]})
                               vector_list.update({(row,col):b_right})

                        elif row == 1: #Equation for any node along top boundary
                                bc_list.update({(row,col):[a_p_top,0,a_ew,a_ns,a_ew]})
                                vector_list.update({(row,col):b_top})

                        elif col ==1: #Equation for any node along left boundary
                                bc_list.update({(row,col):[a_p_ew,a_ns,a_ew,a_ns,0]})
                                vector_list.update({(row,col):b_left})

                        elif col == L_x: #Equation for any node along right boundary
                                bc_list.update({(row,col):[a_p_ew,a_ns,0,a_ns,a_ew]})
                                vector_list.update({(row,col):b_right})

                        elif row == L_y: #Equation for any node along bottom boundary
                                bc_list.update({(row,col):[a_p_bot,a_ns,a_ew,0,a_ew]})
                                
                        else: #Equation for any node not governed by a boundary condition
                                bc_list.update({(row,col):[a_p,a_ns,a_ew,a_ns,a_ew]})
                                            
        
        with open(Matrix_Name+".csv", "w") as csvfile:

                
                for row in range(1,L_y+1): 
                        for col in range(1,L_x+1):

                                node = (row-1)*L_x + col -1 #Converts x,y coordinate of point to node value

                                if bc_list[(row,col)][0] != 0: #Writes 'ap' coefficient value for current principal node for each x,y

                                        csvfile.write(str(bc_list[(row,col)][0])+",")
                                        csvfile.write(str(node)+",")
                                        csvfile.write(str(node)+",")
                                        csvfile.write("\n")

                                if bc_list[(row,col)][1] != 0: #Writes 'an' coefficient value north of principal node

                                        csvfile.write(str(bc_list[(row,col)][1])+",")
                                        csvfile.write(str(node)+",")
                                        csvfile.write(str((row-2)*L_x + col-1)+",")
                                        csvfile.write("\n")

                                if bc_list[(row,col)][2] != 0: #Writes 'ae' coefficient value east of principal node
                                        
                                        csvfile.write(str(bc_list[(row,col)][2])+",")
                                        csvfile.write(str(node)+",")
                                        csvfile.write(str((row-1)*L_x + col)+",")
                                        csvfile.write("\n")

                                if bc_list[(row,col)][3] != 0: #Writes 'as' coefficient value south of principal node

                                        csvfile.write(str(bc_list[(row,col)][3])+",")
                                        csvfile.write(str(node)+",")
                                        csvfile.write(str((row)*L_x + col-1)+",")
                                        csvfile.write("\n")

                                if bc_list[(row,col)][4] != 0: #Writes 'aw' coefficient value west of principal node

                                        csvfile.write(str(bc_list[(row,col)][4])+",")
                                        csvfile.write(str(node)+",")
                                        csvfile.write(str((row-1)*L_x + col-2)+",")
                                        csvfile.write("\n")

        with open(Matrix_Name+" Vector.csv", "w") as csvfile:

                for row in range(1,L_y+1):
                        for col in range(1,L_x+1):
                                
                                if (row,col) in vector_list: #Writes constant for vector b for node x,y
                                        
                                        csvfile.write(str(vector_list[(row,col)])+",")
                                        csvfile.write(str((row-1)*L_x + col-1)+",")
                                        csvfile.write("\n")
                                
#Creating matrices for different mesh sizes
                                        
Matrix_A_Creator("Matrix A Coarse (80)",80,80)
Matrix_A_Creator("Matrix A Fine (160)",160,160)
Matrix_A_Creator("Matrix A Finest (320)",320,320)


t2 = time.time()

print(f"Finished in {round(t2-t1,2)} seconds")
