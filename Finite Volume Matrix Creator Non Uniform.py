### Finite Volume Method Non Uniform Coefficient Matrix Creator
### Will Abbott, 1013944706, Assignment 2

### This code generates a csv file containing data used to construct
### the coefficient matrices and solution vectors for the Finite Volume Method
### equation that was provided in Assignment 2.
### A slightly different approach was used to created non uniform meshing.
### This different approach required a separate code file.

# Initializing modules and constants

import numpy as np

gamma = 20
hf = 10
Phi_x0 = 10
Phi_xL = 100
Phi_ext = 300

#Function is created which takes a matrix name and size as an input.
#Function creates a csv file using equations obtained from using Finite Volume Method (FVM) from Assignment 2.
#Each row of the created csv file corresponds to an ap, an, as, or aw value from the FVM equations.
#Each value also has a corresponding row and column value which will determine it's position in the matrix.
#CSV file will be turned into a coefficient matrix and then solved.
#Code had to be changed to allow delta x and delta y to be specified within each node calculation, rather than taken as constant.

def Matrix_A_Creator(Matrix_Name,L_x,L_y):

    bc_list = {} #Each node x,y is stored in a dictionary with associated coefficients from FVM equation (For matrix A)
    vector_list = {} #Each node x,y is stored in a dictionary with associated constant values from FVM equation (For vector b)

    r = 1.01 #Inflation factor

    dx = 0.5*(L_x)*(1-r)/(1-r**(L_x/2)) #Initial delta x and delta y based on inflation factor
    dy = 0.5*(L_y)*(1-r)/(1-r**(L_y/2))


    for row in range(1,L_y+1): #Cycles through each point in the domain, with coordinate x = col and coordinate y = row
            for col in range(1,L_x+1):  
                    
                    if row == 1 and col == 1: #Equation for top left corner node 
                            bc_list.update({(row,col):[gamma*(3*(dy/dx) +(dx/dy))- hf*dx,
                                                       0,
                                                       -gamma*(dy/dx),
                                                       -gamma*(dx/dy),
                                                       0]})
                            vector_list.update({(row,col):2*gamma*(dy/dx)*Phi_x0 - hf*dx*Phi_ext})
                            
                    elif row == 1 and col == L_x: #Equation for top right corner node
                            bc_list.update({(row,col):[gamma*(3*(dy/dx) +(dx/dy))- hf*dx,
                                                       0,0,
                                                       -gamma*(dx/dy),
                                                       -gamma*(dy/dx)]})
                            vector_list.update({(row,col):2*gamma*(dy/dx)*Phi_xL - hf*dx*Phi_ext})

                    elif row == L_y and col == 1: #Equation for bottom left corner node
                            bc_list.update({(row,col):[gamma*(3*(dy/dx) + (dx/dy)),
                                                       -gamma*(dx/dy),
                                                       -gamma*(dy/dx),
                                                       0,0]})
                            vector_list.update({(row,col):2*gamma*(dy/dx)*Phi_x0})

                    elif row == L_y and col == L_x: #Equation for bottom right corner node
                           bc_list.update({(row,col):[gamma*(3*(dy/dx) + (dx/dy)),
                                                      -gamma*(dx/dy),
                                                      0,0,
                                                      -gamma*(dy/dx)]})
                           vector_list.update({(row,col):2*gamma*(dy/dx)*Phi_xL})

                    elif row == 1: #Equation for any node along top boundary
                            bc_list.update({(row,col):[gamma*(2*(dy/dx) +(dx/dy))- hf*dx,
                                                       0,
                                                       -gamma*(dy/dx),
                                                       -gamma*(dx/dy),
                                                       -gamma*(dy/dx)]})
                            vector_list.update({(row,col):-hf*dx*Phi_ext})

                    elif col ==1: #Equation for any node along left boundary
                            bc_list.update({(row,col):[gamma*(3*(dy/dx) + 2*(dx/dy)),
                                                       -gamma*(dx/dy),
                                                       -gamma*(dy/dx),
                                                       -gamma*(dx/dy),
                                                       0]})
                            vector_list.update({(row,col):2*gamma*(dy/dx)*Phi_x0})

                    elif col == L_x: #Equation for any node along right boundary
                            bc_list.update({(row,col):[gamma*(3*(dy/dx) + 2*(dx/dy)),
                                                       -gamma*(dx/dy),
                                                       0,
                                                       -gamma*(dx/dy),
                                                       -gamma*(dy/dx)]})
                            vector_list.update({(row,col):2*gamma*(dy/dx)*Phi_xL})

                    elif row == L_y: #Equation for any node along bottom boundary
                            bc_list.update({(row,col):[gamma*(2*(dy/dx) + (dx/dy)),
                                                       -gamma*(dx/dy),
                                                       -gamma*(dy/dx),
                                                       0,
                                                       -gamma*(dy/dx)]})
                            
                    else: #Equation for any node not governed by a boundary condition
                            bc_list.update({(row,col):[2*gamma*((dy/dx)+(dx/dy)),
                                                       -gamma*(dx/dy),
                                                       -gamma*(dy/dx),
                                                       -gamma*(dx/dy),
                                                       -gamma*(dy/dx)]})

                    if col <= L_x/2: #Condition to increase or decrease delta x recursively, based on current position along x domain
                            dx = r*dx
                    else:
                            dx = dx/r
                            
            if row <= L_y/2: #Condition to increase or decrease delta y recursively, based on current position along y domain
                    dy = r*dy
            else:
                    dy = dy/r
        
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


#Creating matrices for each mesh size

Matrix_A_Creator("Matrix A Non Uniform Coarse (80)",80,80)
Matrix_A_Creator("Matrix A Non Uniform Fine (160)",160,160)
Matrix_A_Creator("Matrix A Non Uniform Finest (320)",320,320)

print("Done")
